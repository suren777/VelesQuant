
#ifndef ShortRate2FPDE_H
#define ShortRate2FPDE_H

#include <algorithm>

#include <cmath>
#include <map>
#include <memory>
#include <ql/errors.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <string>
#include <vector>
#include <velesquant/models/utility.h>
#include <velesquant/numerics/tri_diag_matrix.h>
#include <velesquant/pde_solvers/cyclic_reduction.h>
#include <velesquant/volatility/lm.h>

namespace velesquant {

/**
 * @brief Two-Factor Short Rate PDE Solver (G2++).
 *
 * This class implements a PDE solver for two-factor short rate models
 * (specifically G2++) using an Alternating Direction Implicit (ADI)
 * Douglas scheme for the 2D PDE.
 *
 * It supports pricing of:
 * - Zero Coupon Bonds
 * - Swaptions (European)
 *
 * The solver handles calibration of the term structure (alphas) to match
 * input discount factors, as well as model parameter calibration to swaptions.
 *
 * @tparam ModelType The underlying short rate model type (e.g.,
 * ShortRate2FModel).
 */
template <typename ModelType> class ShortRate2FPDE {
public:
  /**
   * @brief Construct a new ShortRate2FPDE solver.
   *
   * @param model Shared pointer to the underlying 2-factor short rate model.
   * @param time_step Time step size for the PDE solver (default: 0.020833333333
   * approx 1/48).
   *
   * @throws std::runtime_error If model is null.
   */
public:
  ShortRate2FPDE(std::shared_ptr<ModelType> model,
                 double time_step = 0.020833333333, int grid_points = 21)
      : model_(std::move(model)), time_step_(time_step), Mxy_(grid_points) {
    if (!model_)
      throw std::runtime_error("ShortRate2FPDE: Model cannot be null");

    kappa1_ = model_->getKappa1();
    kappa2_ = model_->getKappa2();
    lambda_ = model_->getLambda();
    timeSigma1s_ = model_->getTimeSigma1s();
    sigma1s_ = model_->getSigma1s();
    timeSigma2s_ = model_->getTimeSigma2s();
    sigma2s_ = model_->getSigma2s();
    timeAlphas_ = model_->getTimeAlphas();
    alphas_ = model_->getAlphas();

    // Initialize grids
    initializeWxy();
  }

  virtual ~ShortRate2FPDE() = default;

  /**
   * @brief Calibrates the model to a set of Swap Quotes and Discount Factors.
   *
   * @param timeDFs Times corresponding to the discount factors.
   * @param DFs Discount factors.
   * @param swapQuotes Vector of Swap definitions to calibrate volatility to.
   * @param optimizer_params Optional map of optimizer parameters (ftol, xtol,
   * gtol, maxfev, epsfcn).
   */
  void calibrator(std::vector<double> timeDFs, std::vector<double> DFs,
                  std::vector<velesquant::defSwap> swapQuotes,
                  std::map<std::string, double> optimizer_params = {}) {
    DFs_ = DFs;
    quoteSwap_ = swapQuotes;
    int nD = timeDFs.size();
    timeAlphas_.resize(nD);
    timeAlphas_ = timeDFs;
    alphas_.resize(nD, 0.015);

    int ns1 = sigma1s_.size();
    int ns2 = sigma2s_.size();
    int n_params = ns1 + ns2 + 3; // no. of model volatility term paremeters
    std::vector<double> x_params(
        n_params); // initial estimate of parameters vector
    for (int i = 0; i < ns1; i++)
      x_params[i] = sigma1s_[i]; // sigma1s initial value
    for (int i = 0; i < ns2; i++)
      x_params[ns1 + i] = sigma2s_[i]; // sigma2s initial value
    x_params[n_params - 1] = kappa1_;  // kappa1 initial value
    x_params[n_params - 2] = kappa2_;  // kappa2 initial value
    x_params[n_params - 3] = lambda_;  // lambda initial value
    int m_obs = quoteSwap_.size();     // no. of observations
    QL_ENSURE(m_obs >= n_params, "too much freedom in Calibration");

    std::vector<double> fvec_vec(m_obs);
    // Default parameters for ShortRate2FPDE
    double ftol = 1e-5;
    double xtol = 1e-5;
    double gtol = 1e-5;
    int maxfev = 200;
    double epsfcn = 1e-10;

    if (optimizer_params.count("ftol"))
      ftol = optimizer_params.at("ftol");
    if (optimizer_params.count("xtol"))
      xtol = optimizer_params.at("xtol");
    if (optimizer_params.count("gtol"))
      gtol = optimizer_params.at("gtol");
    if (optimizer_params.count("maxfev"))
      maxfev = (int)optimizer_params.at("maxfev");
    if (optimizer_params.count("epsfcn"))
      epsfcn = optimizer_params.at("epsfcn");

    std::vector<double> diag(n_params);
    int mode = 1;
    double factor = 1;
    int nprint = 0;
    int info = 0;
    int nfev = 0;
    std::vector<double> fjac(m_obs * n_params);
    int ldfjac = m_obs;
    std::vector<int> ipvt(n_params);
    std::vector<double> qtf(n_params);
    std::vector<double> wa1(n_params);
    std::vector<double> wa2(n_params);
    std::vector<double> wa3(n_params);
    std::vector<double> wa4(m_obs);

    velesquant::lmfcn fcn = [this](int m, int n, double *x, double *fvec,
                                   int *iflag) {
      this->objFcnCalibrator(m, n, x, fvec, iflag);
    };

    velesquant::lmdif(m_obs, n_params, x_params.data(), fvec_vec.data(), ftol,
                      xtol, gtol, maxfev, epsfcn, diag.data(), mode, factor,
                      nprint, &info, &nfev, fjac.data(), ldfjac, ipvt.data(),
                      qtf.data(), wa1.data(), wa2.data(), wa3.data(),
                      wa4.data(), fcn);

    QL_ENSURE(info >= 1 && info <= 4,
              "Model Calibration Fails (info=" + std::to_string(info) + ")");

    // Update local members
    for (int i = 0; i < ns1; i++)
      sigma1s_[i] = std::fabs(x_params[i]);
    for (int i = 0; i < ns2; i++)
      sigma2s_[i] = std::fabs(x_params[ns1 + i]);
    kappa1_ = x_params[n_params - 1];
    kappa2_ = x_params[n_params - 2];
    lambda_ = x_params[n_params - 3];

    // Update Model
    model_->setKappa1(kappa1_);
    model_->setKappa2(kappa2_);
    model_->setLambda(lambda_);
    model_->setSigma1Structure(timeSigma1s_, sigma1s_);
    model_->setSigma2Structure(timeSigma2s_, sigma2s_);
    model_->setAlphaStructure(timeAlphas_, alphas_);
  }

  /**
   * @brief Calculates the price of a Zero Coupon Bond.
   *
   * @param Maturity Time to maturity of the bond.
   * @return double Price of the Zero Coupon Bond.
   */
  double pricingZB(double Maturity) {
    int Nt = int(Maturity / time_step_ + 0.5);
    if (Nt < 2)
      Nt = 2; // Ensuring minimum steps
    buildGrid(Maturity, Nt);
    std::vector<std::vector<double>> payoff;
    payoff.resize(Mxy_);
    for (int x = 0; x < Mxy_; x++) {
      payoff[x].resize(Mxy_);
      for (int y = 0; y < Mxy_; y++)
        payoff[x][y] = 1.0;
    }
    std::vector<std::vector<double>> pv(payoff);
    for (int t = Nt - 2; t >= 0; t--) {
      oneStepBackwardADIDouglas(t, payoff, pv);
      payoff = pv;
    }
    double zeroBond = 0.0;
    for (int x = 0; x < Mxy_; x++)
      if (gridX_[x] >= 0.0) {
        for (int y = 0; y < Mxy_; y++)
          if (gridY_[y] >= 0.0) {
            zeroBond = payoff[x][y];
            break;
          }
        break;
      }
    QL_ENSURE(zeroBond <= 5.0, "pricingZB Functor Fails " << zeroBond);
    return zeroBond;
  }

  /**
   * @brief Calculates the price of a European Swaption.
   *
   * @param Expiry Expiry time of the option.
   * @param Tenor Tenor of the underlying swap.
   * @param Strike Strike rate of the swaption.
   * @param PayFrequency Payment frequency (default 0.5).
   * @return double Price of the Swaption.
   */
  double pricingSwaption(double Expiry, double Tenor, double Strike,
                         double PayFrequency = 0.5) {
    int Nt = int((Expiry + Tenor) / time_step_ + 0.5);
    if (Nt < 2)
      Nt = 2;
    buildGrid((Expiry + Tenor), Nt);
    int iExpiry = int(Expiry / time_step_ + 0.5); // Aligning with grid steps
    // Ensure iExpiry aligns with grid points if needed, or approximate

    int Ncoupon = int(Tenor / PayFrequency + 0.5);
    std::vector<std::vector<double>> payoff;
    payoff.resize(Mxy_);
    for (int x = 0; x < Mxy_; x++) {
      payoff[x].resize(Mxy_);
      for (int y = 0; y < Mxy_; y++)
        payoff[x][y] = -1.0;
    }
    std::vector<std::vector<double>> pv(payoff);
    for (int t = Nt - 2; t >= iExpiry; t--) {
      double couponTime = Expiry + Ncoupon * PayFrequency;
      if (couponTime > gridT_[t] && couponTime <= gridT_[t + 1]) {
        for (int x = 0; x < Mxy_; x++)
          for (int y = 0; y < Mxy_; y++)
            payoff[x][y] -= PayFrequency * Strike;
        Ncoupon--;
      }
      oneStepBackwardADIDouglas(t, payoff, pv);
      payoff = pv;
    }
    for (int x = 0; x < Mxy_; x++)
      for (int y = 0; y < Mxy_; y++)
        payoff[x][y] = std::max(0.0, payoff[x][y] + 1.0);
    for (int t = iExpiry - 1; t >= 0; t--) {
      oneStepBackwardADIDouglas(t, payoff, pv);
      payoff = pv;
    }
    double swaptionValue = 0.0;
    for (int x = 0; x < Mxy_; x++)
      if (gridX_[x] >= 0.0) {
        for (int y = 0; y < Mxy_; y++)
          if (gridY_[y] >= 0.0) {
            swaptionValue = payoff[x][y];
            break;
          }
        break;
      }
    QL_ENSURE(swaptionValue >= -0.0001,
              "pricingSwaption Functor Fails " << swaptionValue);
    return swaptionValue;
  }

  /**
   * @brief Calculates Discount Factors from the model.
   *
   * @param timePoints Vector of times to calculate DFs for.
   * @return std::vector<double> Vector of discount factors.
   */
  std::vector<double> calculateDFs(std::vector<double> &timePoints) {
    int n = timePoints.size();
    std::vector<int> iTs(n);
    for (int i = 0; i < n; i++)
      iTs[i] = int(timePoints[i] / time_step_ + 0.5);

    buildGrid(timePoints[n - 1], iTs[n - 1]);
    std::vector<double> DFs(n, 0.0);
    std::vector<std::vector<double>> inV;
    inV.resize(Mxy_);
    for (int x = 0; x < Mxy_; x++) {
      inV[x].resize(Mxy_);
      for (int y = 0; y < Mxy_; y++)
        inV[x][y] = 0.0;
    }
    for (int x = 1; x < Mxy_ - 1; x++)
      if (gridX_[x] >= 0.0) {
        for (int y = 1; y < Mxy_ - 1; y++)
          if (gridY_[y] >= 0.0) {
            inV[x][y] = 4.0 / ((gridX_[x + 1] - gridX_[x - 1]) *
                               (gridY_[y + 1] - gridY_[y - 1]));
            break;
          }
        break;
      }
    std::vector<std::vector<double>> outV(inV);
    for (int i = 0; i < n; i++) {
      int sT = 0;
      if (i > 0)
        sT = iTs[i - 1] - 1;
      // Make sure loop runs correctly
      if (sT < 0)
        sT = 0;

      for (int t = sT; t < iTs[i] - 1; t++) {
        oneStepForwardADIDouglas(t, inV, outV);
        inV = outV;
      }
      DFs[i] = trapezoidal2D(inV);
      QL_ENSURE(DFs[i] <= 1.001, "pricingDF Functor Fails " << DFs[i]);
    }
    return DFs;
  }

  // Getters
  double getKappa1() { return kappa1_; };
  double getKappa2() { return kappa2_; };
  double getLambda() { return lambda_; };
  std::vector<double> getTimeSigma1s() { return timeSigma1s_; };
  std::vector<double> getSigma1s() { return sigma1s_; };
  std::vector<double> getTimeSigma2s() { return timeSigma2s_; };
  std::vector<double> getSigma2s() { return sigma2s_; };
  std::vector<double> getTimeAlphas() { return timeAlphas_; };
  std::vector<double> getAlphas() { return alphas_; };

private:
  std::shared_ptr<ModelType> model_;
  double kappa1_, kappa2_, lambda_;
  double time_step_; // Configurable time step
  std::vector<double> timeSigma1s_, sigma1s_;
  std::vector<double> timeSigma2s_, sigma2s_;
  mutable std::vector<double> timeAlphas_, alphas_;

  std::vector<double> DFs_;
  std::vector<velesquant::defSwap> quoteSwap_; // Now velesquant::defSwap

  int Mxy_;
  std::vector<double> Wxy;

  void initializeWxy() {
    Wxy.resize(Mxy_);
    // Dynamic grid generation: sinh distribution [-4, 4]
    double Rmax = 4.0;
    double factor = 0.7; // Concentration factor similar to what was likely
                         // used. Or implied.
    // The previous hardcoded grid was [-4 ... 4].
    // Let's implement a standard sinh mesh centered at 0.

    // We want P[0] = -Rmax, P[M-1] = Rmax.
    // P[i] = c * sinh(A * (i - (M-1)/2))

    // Mapping u in [-1, 1] -> x = c * sinh(u * asinh(Rmax/c))
    // c controls concentration. Small c -> more concentration near 0.
    double c = 0.5;
    double asinh_val = std::asinh(Rmax / c);

    for (int i = 0; i < Mxy_; ++i) {
      double u = -1.0 + 2.0 * i / (Mxy_ - 1);
      Wxy[i] = c * std::sinh(u * asinh_val);
    }
    // Ensure exact bounds and 0 center if odd
    Wxy[0] = -Rmax;
    Wxy[Mxy_ - 1] = Rmax;
    if (Mxy_ % 2 != 0)
      Wxy[Mxy_ / 2] = 0.0;

    // Resize workspace vectors
    sol_l_.resize(Mxy_ - 2);
    sol_c_.resize(Mxy_ - 2);
    sol_u_.resize(Mxy_ - 2);
    sol_d_.resize(Mxy_ - 2);
    sol_v_.resize(Mxy_ - 2);
    sol_tempC_.resize(Mxy_ - 2);
    sol_tempD_.resize(Mxy_ - 2);
  }

  void termStructureCalibrator() {
    int n = timeAlphas_.size();
    std::vector<int> iTs(n);
    for (int i = 0; i < n; i++)
      iTs[i] = int(timeAlphas_[i] / time_step_ + 0.5);

    buildGrid(timeAlphas_[n - 1], iTs[n - 1]);
    std::vector<std::vector<double>> inV, outV, lastV;
    inV.resize(Mxy_);
    for (int x = 0; x < Mxy_; x++) {
      inV[x].resize(Mxy_);
      for (int y = 0; y < Mxy_; y++)
        inV[x][y] = 0.0;
    }
    outV = inV;
    for (int x = 1; x < Mxy_ - 1; x++)
      if (gridX_[x] >= 0.0) {
        for (int y = 1; y < Mxy_ - 1; y++)
          if (gridY_[y] >= 0.0) {
            inV[x][y] = 4.0 / ((gridX_[x + 1] - gridX_[x - 1]) *
                               (gridY_[y + 1] - gridY_[y - 1]));
            break;
          }
        break;
      }
    lastV = inV;
    for (int i = 0; i < n; i++) {
      int sT = 0;
      if (i > 0)
        sT = iTs[i - 1] - 1;
      if (sT < 0)
        sT = 0;

      double df = 0.0;
      int niter = 0;
      do {
        niter++;
        inV = lastV;
        if (alphas_[i] == 0.0)
          alphas_[i] = 0.0005;
        alphas_[i] *= 1.001;
        for (int t = sT; t < iTs[i] - 1; t++) {
          oneStepForwardADIDouglas(t, inV, outV);
          inV = outV;
        }
        double dfUP = trapezoidal2D(inV);
        inV = lastV;
        alphas_[i] /= 1.001;
        for (int t = sT; t < iTs[i] - 1; t++) {
          oneStepForwardADIDouglas(t, inV, outV);
          inV = outV;
        }
        df = trapezoidal2D(inV);
        alphas_[i] *= (1 + 0.001 * (DFs_[i] - df) / (dfUP - df));
        alphas_[i] = std::max(-0.0050, alphas_[i]);
      } while (std::fabs(1.0 - df / DFs_[i]) >= 1.0E-6 && niter < 7);
      inV = lastV;
      for (int t = sT; t < iTs[i] - 1; t++) {
        oneStepForwardADIDouglas(t, inV, outV);
        inV = outV;
      }
      lastV = inV;
    };
  }

  void objFcnCalibrator(int m, int n, double *x, double *fvec,
                        int * /*iflag*/) {
    int ns1 = sigma1s_.size();
    int ns2 = sigma2s_.size();
    for (int i = 0; i < ns1; i++)
      sigma1s_[i] = std::fabs(x[i]);
    for (int i = 0; i < ns2; i++)
      sigma2s_[i] = std::fabs(x[ns1 + i]);
    kappa1_ = x[n - 1];
    kappa2_ = x[n - 2];
    lambda_ = x[n - 3];
    termStructureCalibrator();
    for (int i = 0; i < m; i++) {
      double model =
          pricingSwaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                          quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
      double market = quoteSwap_[i].Value;
      fvec[i] = model - market;
    }
  }

  mutable std::vector<double> gridT_;
  mutable std::vector<double> gridX_, gridY_;
  mutable std::vector<double> sol_l_, sol_c_, sol_u_, sol_d_, sol_v_,
      sol_tempC_, sol_tempD_;

  void buildGrid(double Time, int Nt = 5000) {
    gridT_.resize(Nt);
    double delT = Time / (Nt - 1);
    gridT_[0] = 0.0;
    for (int t = 1; t < Nt - 1; t++)
      gridT_[t] = t * delT;
    gridT_[Nt - 1] = Time;
    gridX_.resize(Mxy_);
    double avgSigma = 0.0;
    int nP = timeSigma1s_.size();
    for (int n = 0; n < nP; n++) {
      if (n == 0)
        avgSigma += sigma1s_[n] * sigma1s_[n] * timeSigma1s_[n];
      else
        avgSigma +=
            sigma1s_[n] * sigma1s_[n] * (timeSigma1s_[n] - timeSigma1s_[n - 1]);
    }
    avgSigma = sqrt(avgSigma / timeSigma1s_[nP - 1]) *
               sqrt((1 - exp(-2 * kappa1_ * Time)) / (2 * kappa1_));
    for (int x = 0; x < Mxy_; x++)
      gridX_[x] = Wxy[x] * avgSigma;
    gridY_.resize(Mxy_);
    avgSigma = 0.0;
    nP = timeSigma2s_.size();
    for (int n = 0; n < nP; n++) {
      if (n == 0)
        avgSigma += sigma2s_[n] * sigma2s_[n] * timeSigma2s_[n];
      else
        avgSigma +=
            sigma2s_[n] * sigma2s_[n] * (timeSigma2s_[n] - timeSigma2s_[n - 1]);
    }
    avgSigma = sqrt(avgSigma / timeSigma2s_[nP - 1]) *
               sqrt((1 - exp(-2 * kappa2_ * Time)) / (2 * kappa2_));
    for (int x = 0; x < Mxy_; x++)
      gridY_[x] = Wxy[x] * avgSigma;
  }

  void oneStepBackwardADIDouglas(int t,
                                 const std::vector<std::vector<double>> &inM,
                                 std::vector<std::vector<double>> &outM) {
    oneStepBackwardExplicit(t, inM, outM);
    std::vector<std::vector<double>> midM(outM);
    oneStepBackwardDouglasX(t, inM, midM, outM);
    midM = outM;
    oneStepBackwardDouglasY(t, inM, midM, outM);
  }

  void oneStepBackwardExplicit(int t,
                               const std::vector<std::vector<double>> &inM,
                               std::vector<std::vector<double>> &outM) {
    double te = gridT_[t + 1];
    double ti = gridT_[t];
    double delT = te - ti;
    double tm = 0.5 * (te + ti);
    double sigma1 = whichValue(tm, timeSigma1s_, sigma1s_);
    double difu1 = sigma1 * sigma1;
    double sigma2 = whichValue(tm, timeSigma2s_, sigma2s_);
    double difu2 = sigma2 * sigma2;
    double alpha = whichValue(tm, timeAlphas_, alphas_);

    // j=0
    double X0 = gridX_[0];
    double X1 = gridX_[1];
    double X2 = gridX_[2];
    double conv1 = kappa1_ * X0;
    // j=0, k=0 lower boundary condition
    double Y0 = gridY_[0];
    double Y1 = gridY_[1];
    double Y2 = gridY_[2];
    double conv2 = lambda_ * X0 + kappa2_ * Y0;
    double rate = alpha + X0 + Y0;
    outM[0][0] =
        (1.0 / delT - rate) * inM[0][0] +
        (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][0] +
        (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][0] +
        (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][0] +
        (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[0][0] +
        (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[0][1] +
        (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[0][2];
    outM[0][0] = delT * outM[0][0];
    // j=0, k\=0 middle range
    double Y, Yu, Yl;
    for (int k = 1; k < Mxy_ - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = lambda_ * X0 + kappa2_ * Y;
      rate = alpha + X0 + Y - kappa1_ - kappa2_;

      outM[0][k] =
          (1.0 / delT - rate) * inM[0][k] +
          (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
              inM[0][k] +
          (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
          (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
          (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[0][k - 1] +
          (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
              inM[0][k] +
          (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[0][k + 1];
      outM[0][k] = delT * outM[0][k];
    }
    // j=0, k=M-1 upper boundary condition
    double Y3l = gridY_[Mxy_ - 3];
    double Y2l = gridY_[Mxy_ - 2];
    double Y1l = gridY_[Mxy_ - 1];
    conv2 = lambda_ * X0 + kappa2_ * Y1l;
    rate = alpha + X0 + Y1l;
    outM[0][Mxy_ - 1] =
        (1.0 / delT - rate) * inM[0][Mxy_ - 1] +
        (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][Mxy_ - 1] +
        (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][Mxy_ - 1] +
        (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) *
            inM[2][Mxy_ - 1] +
        (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
            inM[0][Mxy_ - 3] +
        (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
            inM[0][Mxy_ - 2] +
        (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
            inM[0][Mxy_ - 1];
    outM[0][Mxy_ - 1] = delT * outM[0][Mxy_ - 1];

    // j\=0
    double X, Xu, Xl;
    for (int j = 1; j < Mxy_ - 1; j++) {
      X = gridX_[j];
      Xu = gridX_[j + 1];
      Xl = gridX_[j - 1];
      conv1 = kappa1_ * X;
      // j\=0, k=0 lower boundary condition
      conv2 = lambda_ * X + kappa2_ * Y0;
      rate = alpha + X + Y0;
      outM[j][0] =
          (1.0 / delT - rate) * inM[j][0] +
          (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][0] +
          (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
              inM[j][0] +
          (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][0] +
          (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
              inM[j][0] +
          (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
          (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2];
      outM[j][0] = delT * outM[j][0];
      // j\=0, k\=0 middle range
      for (int k = 1; k < Mxy_ - 1; k++) {
        Y = gridY_[k];
        Yu = gridY_[k + 1];
        Yl = gridY_[k - 1];
        conv2 = lambda_ * X + kappa2_ * Y;
        rate = alpha + X + Y;
        outM[j][k] =
            (1.0 / delT - rate) * inM[j][k] +
            (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
            (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
                inM[j][k] +
            (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
            (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
            (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
                inM[j][k] +
            (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1];
        outM[j][k] = delT * outM[j][k];
      }
      // j\=0, k=M-1 upper boundary condition
      conv2 = lambda_ * X + kappa2_ * Y1l;
      rate = alpha + X + Y1l;
      outM[j][Mxy_ - 1] = (1.0 / delT - rate) * inM[j][Mxy_ - 1] +
                          (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) *
                              inM[j - 1][Mxy_ - 1] +
                          (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) /
                              (Xu - X) * inM[j][Mxy_ - 1] +
                          (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) *
                              inM[j + 1][Mxy_ - 1] +
                          (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) /
                              (Y1l - Y3l) * inM[j][Mxy_ - 3] +
                          (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) /
                              (Y1l - Y2l) * inM[j][Mxy_ - 2] +
                          (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) /
                              (Y1l - Y2l) / (Y1l - Y3l) * inM[j][Mxy_ - 1];
      outM[j][Mxy_ - 1] = delT * outM[j][Mxy_ - 1];
    }

    // j=M-1
    double X3l = gridX_[Mxy_ - 3];
    double X2l = gridX_[Mxy_ - 2];
    double X1l = gridX_[Mxy_ - 1];
    conv1 = kappa1_ * X1l;
    // j=M-1, k=0 lower boundary condition
    conv2 = lambda_ * X1l + kappa2_ * Y0;
    rate = alpha + X1l + Y0;
    outM[Mxy_ - 1][0] =
        (1.0 / delT - rate) * inM[Mxy_ - 1][0] +
        (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
            inM[Mxy_ - 3][0] +
        (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
            inM[Mxy_ - 2][0] +
        (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
            inM[Mxy_ - 1][0] +
        (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[Mxy_ - 1][0] +
        (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[Mxy_ - 1][1] +
        (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[Mxy_ - 1][2];
    outM[Mxy_ - 1][0] = delT * outM[Mxy_ - 1][0];
    // j=M-1, k\=0 middle range
    for (int k = 1; k < Mxy_ - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = lambda_ * X1l + kappa2_ * Y;
      rate = alpha + X1l + Y;
      outM[Mxy_ - 1][k] = (1.0 / delT - rate) * inM[Mxy_ - 1][k] +
                          (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) /
                              (X1l - X3l) * inM[Mxy_ - 3][k] +
                          (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) /
                              (X1l - X2l) * inM[Mxy_ - 2][k] +
                          (-conv1 * (2 * X1l - X2l - X3l) + difu1) /
                              (X1l - X2l) / (X1l - X3l) * inM[Mxy_ - 1][k] +
                          (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) *
                              inM[Mxy_ - 1][k - 1] +
                          (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) /
                              (Yu - Y) * inM[Mxy_ - 1][k] +
                          (-conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) *
                              inM[Mxy_ - 1][k + 1];
      outM[Mxy_ - 1][k] = delT * outM[Mxy_ - 1][k];
    }
    // j=M-1, k=M-1 upper boundary condition
    conv2 = lambda_ * X1l + kappa2_ * Y1l;
    rate = alpha + X1l + Y1l;
    outM[Mxy_ - 1][Mxy_ - 1] =
        (1.0 / delT - rate) * inM[Mxy_ - 1][Mxy_ - 1] +
        (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
            inM[Mxy_ - 3][Mxy_ - 1] +
        (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
            inM[Mxy_ - 2][Mxy_ - 1] +
        (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
            inM[Mxy_ - 1][Mxy_ - 1] +
        (-conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
            inM[Mxy_ - 1][Mxy_ - 3] +
        (conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
            inM[Mxy_ - 1][Mxy_ - 2] +
        (-conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
            inM[Mxy_ - 1][Mxy_ - 1];
    outM[Mxy_ - 1][Mxy_ - 1] = delT * outM[Mxy_ - 1][Mxy_ - 1];
  }

  void oneStepForwardDouglasX(int T,
                              const std::vector<std::vector<double>> &inM,
                              const std::vector<std::vector<double>> &midM,
                              std::vector<std::vector<double>> &outM) {
    double Te = gridT_[T + 1];
    double Ti = gridT_[T];
    double delT = Te - Ti;
    double Tm = 0.5 * (Te + Ti);
    double sigma1 = whichValue(Tm, timeSigma1s_, sigma1s_);
    double difu1 = 0.5 * sigma1 * sigma1;
    double alpha = whichValue(Tm, timeAlphas_, alphas_);
    double conv1, rate;

    // for all y (k=0,...,M-1)
    for (int k = 0; k < Mxy_; k++) {
      double Y = gridY_[k];
      // x in the middle range
      for (int j = 1; j < Mxy_ - 1; j++) {
        double X = gridX_[j];
        double Xu = gridX_[j + 1];
        double Xl = gridX_[j - 1];
        conv1 = 0.5 * kappa1_ * X;
        rate = 0.25 * (alpha + X + Y);
        sol_l_[j - 1] = (-conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl);
        sol_c_[j - 1] =
            1.0 / delT + rate +
            (conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X);
        sol_u_[j - 1] = (conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl);
        sol_d_[j - 1] =
            (-conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
            (rate +
             (conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X)) *
                inM[j][k] +
            (conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
            midM[j][k] / delT;
      }
      // j=0 lower boundary condition
      double X0 = gridX_[0];
      double X1 = gridX_[1];
      double X2 = gridX_[2];
      conv1 = 0.5 * kappa1_ * X0;
      rate =
          0.25 * (alpha + X0 + Y - kappa1_ - kappa2_); // Correct this one too

      double k0 =
          1.0 / delT + rate +
          (-conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0);
      double k1 = (conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1);
      double k2 = (-conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0);
      double d0 =
          (rate +
           (-conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0)) *
              inM[0][k] +
          (conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
          (-conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
          midM[0][k] / delT;
      sol_c_[0] = sol_c_[0] - sol_l_[0] * k1 / k0;
      sol_u_[0] = sol_u_[0] - sol_l_[0] * k2 / k0;
      sol_d_[0] = sol_d_[0] - sol_l_[0] * d0 / k0;
      sol_l_[0] = 0.0;
      // j=M-1 upper boundary condition
      double X3 = gridX_[Mxy_ - 3];
      double X2l = gridX_[Mxy_ - 2];
      double X1l = gridX_[Mxy_ - 1];
      conv1 = 0.5 * kappa1_ * X1l;
      rate = 0.25 * (alpha + X1l + Y);
      double kl3 = (conv1 * (X1l - X2l) - difu1) / (X2l - X3) / (X1l - X3);
      double kl2 = (-conv1 * (X1l - X3) + difu1) / (X2l - X3) / (X1l - X2l);
      double kl1 =
          1.0 / delT + rate +
          (conv1 * (2.0 * X1l - X2l - X3) - difu1) / (X1l - X2l) / (X1l - X3);
      double dl1 = (conv1 * (X1l - X2l) - difu1) / (X2l - X3) / (X1l - X3) *
                       inM[Mxy_ - 3][k] +
                   (-conv1 * (X1l - X3) + difu1) / (X2l - X3) / (X1l - X2l) *
                       inM[Mxy_ - 2][k] +
                   (rate + (conv1 * (2.0 * X1l - X2l - X3) - difu1) /
                               (X1l - X2l) / (X1l - X3)) *
                       inM[Mxy_ - 1][k] +
                   midM[Mxy_ - 1][k] / delT;
      sol_l_[Mxy_ - 3] = sol_l_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl3 / kl1;
      sol_c_[Mxy_ - 3] = sol_c_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl2 / kl1;
      sol_d_[Mxy_ - 3] = sol_d_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * dl1 / kl1;
      sol_u_[Mxy_ - 3] = 0.0;
      velesquant::TriDiagonalSolve(Mxy_ - 2, sol_l_, sol_c_, sol_u_, sol_d_,
                                   sol_v_, sol_tempC_, sol_tempD_);
      for (int j = 1; j < Mxy_ - 1; j++)
        outM[j][k] = sol_v_[j - 1];
      // update j=0 lower boundary
      outM[0][k] = (d0 - k1 * outM[1][k] - k2 * outM[2][k]) / k0;
      // update j=M-1 upper boundary
      outM[Mxy_ - 1][k] =
          (dl1 - kl3 * outM[Mxy_ - 3][k] - kl2 * outM[Mxy_ - 2][k]) / kl1;
    }
  }

  void oneStepForwardDouglasY(int T,
                              const std::vector<std::vector<double>> &inM,
                              const std::vector<std::vector<double>> &midM,
                              std::vector<std::vector<double>> &outM) {
    double Te = gridT_[T + 1];
    double Ti = gridT_[T];
    double delT = Te - Ti;
    double Tm = 0.5 * (Te + Ti);
    double sigma2 = whichValue(Tm, timeSigma2s_, sigma2s_);
    double difu2 = 0.5 * sigma2 * sigma2;
    double alpha = whichValue(Tm, timeAlphas_, alphas_);
    double conv2, rate;

    // for all x (j=0,...,M-1)
    for (int j = 0; j < Mxy_; j++) {
      double X = gridX_[j];
      // y in the middle range
      for (int k = 1; k < Mxy_ - 1; k++) {
        double Y = gridY_[k];
        double Yu = gridY_[k + 1];
        double Yl = gridY_[k - 1];
        conv2 = 0.5 * (lambda_ * X + kappa2_ * Y);
        rate = 0.25 * (alpha + X + Y);
        sol_l_[k - 1] = (-conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl);
        sol_c_[k - 1] =
            1.0 / delT + rate +
            (conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y);
        sol_u_[k - 1] = (conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl);
        sol_d_[k - 1] =
            (-conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
            (rate +
             (conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y)) *
                inM[j][k] +
            (conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1] +
            midM[j][k] / delT;
      }
      // k=0 lower boundary condition
      double Y0 = gridY_[0];
      double Y1 = gridY_[1];
      double Y2 = gridY_[2];
      conv2 = 0.5 * (lambda_ * X + kappa2_ * Y0);
      rate = 0.25 * (alpha + X + Y0);
      double k0 =
          1.0 / delT + rate +
          (-conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0);
      double k1 = (conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1);
      double k2 = (-conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0);
      double d0 =
          (rate +
           (-conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0)) *
              inM[j][0] +
          (conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
          (-conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2] +
          midM[j][0] / delT;
      sol_c_[0] = sol_c_[0] - sol_l_[0] * k1 / k0;
      sol_u_[0] = sol_u_[0] - sol_l_[0] * k2 / k0;
      sol_d_[0] = sol_d_[0] - sol_l_[0] * d0 / k0;
      sol_l_[0] = 0.0;
      // k=M-1 upper boundary condition
      double Y3 = gridY_[Mxy_ - 3];
      double Y2l = gridY_[Mxy_ - 2];
      double Y1l = gridY_[Mxy_ - 1];
      conv2 = 0.5 * (lambda_ * X + kappa2_ * Y1l);
      rate = 0.25 * (alpha + X + Y1l);
      double kl3 = (conv2 * (Y1l - Y2l) - difu2) / (Y2l - Y3) / (Y1l - Y3);
      double kl2 = (-conv2 * (Y1l - Y3) + difu2) / (Y2l - Y3) / (Y1l - Y2l);
      double kl1 =
          1.0 / delT + rate +
          (conv2 * (2.0 * Y1l - Y2l - Y3) - difu2) / (Y1l - Y2l) / (Y1l - Y3);
      double dl1 = (conv2 * (Y1l - Y2l) - difu2) / (Y2l - Y3) / (Y1l - Y3) *
                       inM[j][Mxy_ - 3] +
                   (-conv2 * (Y1l - Y3) + difu2) / (Y2l - Y3) / (Y1l - Y2l) *
                       inM[j][Mxy_ - 2] +
                   (rate + (conv2 * (2.0 * Y1l - Y2l - Y3) - difu2) /
                               (Y1l - Y2l) / (Y1l - Y3)) *
                       inM[j][Mxy_ - 1] +
                   midM[j][Mxy_ - 1] / delT;
      sol_l_[Mxy_ - 3] = sol_l_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl3 / kl1;
      sol_c_[Mxy_ - 3] = sol_c_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl2 / kl1;
      sol_d_[Mxy_ - 3] = sol_d_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * dl1 / kl1;
      sol_u_[Mxy_ - 3] = 0.0;
      velesquant::TriDiagonalSolve(Mxy_ - 2, sol_l_, sol_c_, sol_u_, sol_d_,
                                   sol_v_, sol_tempC_, sol_tempD_);
      for (int k = 1; k < Mxy_ - 1; k++)
        outM[j][k] = sol_v_[k - 1];
      // update k=0 lower boundary
      outM[j][0] = (d0 - k1 * outM[j][1] - k2 * outM[j][2]) / k0;
      // update k=M-1 upper boundary
      outM[j][Mxy_ - 1] =
          (dl1 - kl3 * outM[j][Mxy_ - 3] - kl2 * outM[j][Mxy_ - 2]) / kl1;
    }
  }

  double trapezoidal2D(std::vector<std::vector<double>> &inV) {
    double value = 0.0;
    for (int x = 1; x < Mxy_; x++)
      for (int y = 1; y < Mxy_; y++)
        value +=
            0.25 *
            (inV[x][y] + inV[x][y - 1] + inV[x - 1][y] + inV[x - 1][y - 1]) *
            (gridX_[x] - gridX_[x - 1]) * (gridY_[y] - gridY_[y - 1]);
    return value;
  }

  double whichValue(double t, const std::vector<double> &times,
                    const std::vector<double> &values) {
    int n = times.size();
    if (t < times[0])
      return values[0];
    if (t >= times[n - 1])
      return values[n - 1];
    for (int i = 1; i < n; i++) {
      if ((t >= times[i - 1]) && (t < times[i]))
        return values[i];
    }
    return values.back();
  }

  void oneStepForwardADIDouglas(int t,
                                const std::vector<std::vector<double>> &inM,
                                std::vector<std::vector<double>> &outM) {
    oneStepForwardExplicit(t, inM, outM);
    std::vector<std::vector<double>> midM(outM);
    oneStepForwardDouglasX(t, inM, midM, outM);
    midM = outM;
    oneStepForwardDouglasY(t, inM, midM, outM);
  }

  void oneStepForwardExplicit(int t,
                              const std::vector<std::vector<double>> &inM,
                              std::vector<std::vector<double>> &outM) {
    double te = gridT_[t + 1];
    double ti = gridT_[t];
    double delT = te - ti;
    double tm = 0.5 * (te + ti);
    double sigma1 = whichValue(tm, timeSigma1s_, sigma1s_);
    double difu1 = sigma1 * sigma1;
    double sigma2 = whichValue(tm, timeSigma2s_, sigma2s_);
    double difu2 = sigma2 * sigma2;
    double alpha = whichValue(tm, timeAlphas_, alphas_);

    // j=0
    double X0 = gridX_[0];
    double X1 = gridX_[1];
    double X2 = gridX_[2];
    double conv1 = -kappa1_ * X0;
    // j=0, k=0 lower boundary condition
    double Y0 = gridY_[0];
    double Y1 = gridY_[1];
    double Y2 = gridY_[2];
    double conv2 = -(lambda_ * X0 + kappa2_ * Y0);
    double rate = alpha + X0 + Y0 - kappa1_ - kappa2_;
    outM[0][0] =
        (1.0 / delT - rate) * inM[0][0] +
        (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][0] +
        (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][0] +
        (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][0] +
        (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[0][0] +
        (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[0][1] +
        (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[0][2];
    outM[0][0] = delT * outM[0][0];
    // j=0, k\=0 middle range
    double Y, Yu, Yl;
    for (int k = 1; k < Mxy_ - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = -(lambda_ * X0 + kappa2_ * Y);
      rate = alpha + X0 + Y - kappa1_ - kappa2_;
      outM[0][k] =
          (1.0 / delT - rate) * inM[0][k] +
          (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
              inM[0][k] +
          (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
          (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
          (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[0][k - 1] +
          (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
              inM[0][k] +
          (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[0][k + 1];
      outM[0][k] = delT * outM[0][k];
    }
    // j=0, k=M-1 upper boundary condition
    double Y3l = gridY_[Mxy_ - 3];
    double Y2l = gridY_[Mxy_ - 2];
    double Y1l = gridY_[Mxy_ - 1];
    conv2 = -(lambda_ * X0 + kappa2_ * Y1l);
    rate = alpha + X0 + Y1l - kappa1_ - kappa2_;
    outM[0][Mxy_ - 1] =
        (1.0 / delT - rate) * inM[0][Mxy_ - 1] +
        (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][Mxy_ - 1] +
        (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][Mxy_ - 1] +
        (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) *
            inM[2][Mxy_ - 1] +
        (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
            inM[0][Mxy_ - 3] +
        (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
            inM[0][Mxy_ - 2] +
        (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
            inM[0][Mxy_ - 1];
    outM[0][Mxy_ - 1] = delT * outM[0][Mxy_ - 1];

    // j\=0
    double X, Xu, Xl;
    for (int j = 1; j < Mxy_ - 1; j++) {
      X = gridX_[j];
      Xu = gridX_[j + 1];
      Xl = gridX_[j - 1];
      conv1 = -kappa1_ * X;
      // j\=0, k=0 lower boundary condition
      conv2 = -(lambda_ * X + kappa2_ * Y0);
      rate = alpha + X + Y0 - kappa1_ - kappa2_;
      outM[j][0] =
          (1.0 / delT - rate) * inM[j][0] +
          (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][0] +
          (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
              inM[j][0] +
          (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][0] +
          (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
              inM[j][0] +
          (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
          (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2];
      outM[j][0] = delT * outM[j][0];
      // j\=0, k\=0 middle range
      for (int k = 1; k < Mxy_ - 1; k++) {
        Y = gridY_[k];
        Yu = gridY_[k + 1];
        Yl = gridY_[k - 1];
        conv2 = -(lambda_ * X + kappa2_ * Y);
        rate = alpha + X + Y - kappa1_ - kappa2_;
        outM[j][k] =
            (1.0 / delT - rate) * inM[j][k] +
            (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
            (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
                inM[j][k] +
            (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
            (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
            (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
                inM[j][k] +
            (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1];
        outM[j][k] = delT * outM[j][k];
      }
      // j\=0, k=M-1 upper boundary condition
      conv2 = -(lambda_ * X + kappa2_ * Y1l);
      rate = alpha + X + Y1l - kappa1_ - kappa2_;
      outM[j][Mxy_ - 1] = (1.0 / delT - rate) * inM[j][Mxy_ - 1] +
                          (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) *
                              inM[j - 1][Mxy_ - 1] +
                          (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) /
                              (Xu - X) * inM[j][Mxy_ - 1] +
                          (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) *
                              inM[j + 1][Mxy_ - 1] +
                          (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) /
                              (Y1l - Y3l) * inM[j][Mxy_ - 3] +
                          (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) /
                              (Y1l - Y2l) * inM[j][Mxy_ - 2] +
                          (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) /
                              (Y1l - Y2l) / (Y1l - Y3l) * inM[j][Mxy_ - 1];
      outM[j][Mxy_ - 1] = delT * outM[j][Mxy_ - 1];
    }

    // j=M-1
    double X3l = gridX_[Mxy_ - 3];
    double X2l = gridX_[Mxy_ - 2];
    double X1l = gridX_[Mxy_ - 1];
    conv1 = -kappa1_ * X1l;
    // j=M-1, k=0 lower boundary condition
    conv2 = -(lambda_ * X1l + kappa2_ * Y0);
    rate = alpha + X1l + Y0 - kappa1_ - kappa2_;
    outM[Mxy_ - 1][0] =
        (1.0 / delT - rate) * inM[Mxy_ - 1][0] +
        (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
            inM[Mxy_ - 3][0] +
        (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
            inM[Mxy_ - 2][0] +
        (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
            inM[Mxy_ - 1][0] +
        (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[Mxy_ - 1][0] +
        (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[Mxy_ - 1][1] +
        (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[Mxy_ - 1][2];
    outM[Mxy_ - 1][0] = delT * outM[Mxy_ - 1][0];
    // j=M-1, k\=0 middle range
    for (int k = 1; k < Mxy_ - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = -(lambda_ * X1l + kappa2_ * Y);
      rate = alpha + X1l + Y - kappa1_ - kappa2_;
      outM[Mxy_ - 1][k] = (1.0 / delT - rate) * inM[Mxy_ - 1][k] +
                          (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) /
                              (X1l - X3l) * inM[Mxy_ - 3][k] +
                          (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) /
                              (X1l - X2l) * inM[Mxy_ - 2][k] +
                          (-conv1 * (2 * X1l - X2l - X3l) + difu1) /
                              (X1l - X2l) / (X1l - X3l) * inM[Mxy_ - 1][k] +
                          (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) *
                              inM[Mxy_ - 1][k - 1] +
                          (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) /
                              (Yu - Y) * inM[Mxy_ - 1][k] +
                          (-conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) *
                              inM[Mxy_ - 1][k + 1];
      outM[Mxy_ - 1][k] = delT * outM[Mxy_ - 1][k];
    }
    // j=M-1, k=M-1 upper boundary condition
    conv2 = -(lambda_ * X1l + kappa2_ * Y1l);
    rate = alpha + X1l + Y1l - kappa1_ - kappa2_;
    outM[Mxy_ - 1][Mxy_ - 1] =
        (1.0 / delT - rate) * inM[Mxy_ - 1][Mxy_ - 1] +
        (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
            inM[Mxy_ - 3][Mxy_ - 1] +
        (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
            inM[Mxy_ - 2][Mxy_ - 1] +
        (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
            inM[Mxy_ - 1][Mxy_ - 1] +
        (-conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
            inM[Mxy_ - 1][Mxy_ - 3] +
        (conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
            inM[Mxy_ - 1][Mxy_ - 2] +
        (-conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
            inM[Mxy_ - 1][Mxy_ - 1];
    outM[Mxy_ - 1][Mxy_ - 1] = delT * outM[Mxy_ - 1][Mxy_ - 1];
  }

  void oneStepBackwardDouglasX(int T,
                               const std::vector<std::vector<double>> &inM,
                               const std::vector<std::vector<double>> &midM,
                               std::vector<std::vector<double>> &outM) {
    // Reuse member vectors instead of allocating new ones
    // std::vector<double> l(Mxy_ - 2), c(Mxy_ - 2), u(Mxy_ - 2), d(Mxy_ - 2),
    // V(Mxy_ - 2);
    double Te = gridT_[T + 1];
    double Ti = gridT_[T];
    double delT = Te - Ti;
    double Tm = 0.5 * (Te + Ti);
    double sigma1 = whichValue(Tm, timeSigma1s_, sigma1s_);
    double difu1 = 0.5 * sigma1 * sigma1;
    double alpha = whichValue(Tm, timeAlphas_, alphas_);
    double conv1, rate;

    // for all y (k=0,...,M-1)
    for (int k = 0; k < Mxy_; k++) {
      double Y = gridY_[k];
      // x in the middle range
      for (int j = 1; j < Mxy_ - 1; j++) {
        double X = gridX_[j];
        double Xu = gridX_[j + 1];
        double Xl = gridX_[j - 1];
        conv1 = -0.5 * kappa1_ * X;
        rate = 0.25 * (alpha + X + Y - kappa1_ - kappa2_);
        sol_l_[j - 1] = (-conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl);
        sol_c_[j - 1] =
            1.0 / delT + rate +
            (conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X);
        sol_u_[j - 1] = (conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl);
        sol_d_[j - 1] =
            (-conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
            (rate +
             (conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X)) *
                inM[j][k] +
            (conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
            midM[j][k] / delT;
      }
      // j=0 lower boundary condition
      double X0 = gridX_[0];
      double X1 = gridX_[1];
      double X2 = gridX_[2];
      conv1 = -0.5 * kappa1_ * X0;
      rate = 0.25 * (alpha + X0 + Y - kappa1_ - kappa2_);
      double k0 =
          1.0 / delT + rate +
          (-conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0);
      double k1 = (conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1);
      double k2 = (-conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0);
      double d0 =
          (rate +
           (-conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0)) *
              inM[0][k] +
          (conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
          (-conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
          midM[0][k] / delT;
      sol_c_[0] = sol_c_[0] - sol_l_[0] * k1 / k0;
      sol_u_[0] = sol_u_[0] - sol_l_[0] * k2 / k0;
      sol_d_[0] = sol_d_[0] - sol_l_[0] * d0 / k0;
      sol_l_[0] = 0.0;
      // j=M-1 upper boundary condition
      double X3 = gridX_[Mxy_ - 3];
      X2 = gridX_[Mxy_ - 2];
      X1 = gridX_[Mxy_ - 1];
      conv1 = -0.5 * kappa1_ * X1;
      rate = 0.25 * (alpha + X1 + Y - kappa1_ - kappa2_);
      double kl3 = (conv1 * (X1 - X2) - difu1) / (X2 - X3) / (X1 - X3);
      double kl2 = (-conv1 * (X1 - X3) + difu1) / (X2 - X3) / (X1 - X2);
      double kl1 =
          1.0 / delT + rate +
          (conv1 * (2.0 * X1 - X2 - X3) - difu1) / (X1 - X2) / (X1 - X3);
      double dl1 = (conv1 * (X1 - X2) - difu1) / (X2 - X3) / (X1 - X3) *
                       inM[Mxy_ - 3][k] +
                   (-conv1 * (X1 - X3) + difu1) / (X2 - X3) / (X1 - X2) *
                       inM[Mxy_ - 2][k] +
                   (rate + (conv1 * (2.0 * X1 - X2 - X3) - difu1) / (X1 - X2) /
                               (X1 - X3)) *
                       inM[Mxy_ - 1][k] +
                   midM[Mxy_ - 1][k] / delT;
      sol_l_[Mxy_ - 3] = sol_l_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl3 / kl1;
      sol_c_[Mxy_ - 3] = sol_c_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl2 / kl1;
      sol_d_[Mxy_ - 3] = sol_d_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * dl1 / kl1;
      sol_u_[Mxy_ - 3] = 0.0;
      velesquant::TriDiagonalSolve(Mxy_ - 2, sol_l_, sol_c_, sol_u_, sol_d_,
                                   sol_v_, sol_tempC_, sol_tempD_);
      for (int j = 1; j < Mxy_ - 1; j++)
        outM[j][k] = sol_v_[j - 1];
      // update j=0 lower boundary
      outM[0][k] = (d0 - k1 * outM[1][k] - k2 * outM[2][k]) / k0;
      // update j=M-1 upper boundary
      outM[Mxy_ - 1][k] =
          (dl1 - kl3 * outM[Mxy_ - 3][k] - kl2 * outM[Mxy_ - 2][k]) / kl1;
    }
  }

  void oneStepBackwardDouglasY(int T,
                               const std::vector<std::vector<double>> &inM,
                               const std::vector<std::vector<double>> &midM,
                               std::vector<std::vector<double>> &outM) {
    // Reuse member vectors
    // std::vector<double> l(Mxy_ - 2), c(Mxy_ - 2), u(Mxy_ - 2), d(Mxy_ - 2),
    // V(Mxy_ - 2);
    double Te = gridT_[T + 1];
    double Ti = gridT_[T];
    double delT = Te - Ti;
    double Tm = 0.5 * (Te + Ti);
    double sigma2 = whichValue(Tm, timeSigma2s_, sigma2s_);
    double difu2 = 0.5 * sigma2 * sigma2;
    double alpha = whichValue(Tm, timeAlphas_, alphas_);
    double conv2, rate;

    // for all x (j=0,...,M-1)
    for (int j = 0; j < Mxy_; j++) {
      double X = gridX_[j];
      // y in the middle range
      for (int k = 1; k < Mxy_ - 1; k++) {
        double Y = gridY_[k];
        double Yu = gridY_[k + 1];
        double Yl = gridY_[k - 1];
        conv2 = -0.5 * (lambda_ * X + kappa2_ * Y);
        rate = 0.25 * (alpha + X + Y - kappa1_ - kappa2_);
        sol_l_[k - 1] = (-conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl);
        sol_c_[k - 1] =
            1.0 / delT + rate +
            (conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y);
        sol_u_[k - 1] = (conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl);
        sol_d_[k - 1] =
            (-conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
            (rate +
             (conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y)) *
                inM[j][k] +
            (conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1] +
            midM[j][k] / delT;
      }
      // k=0 lower boundary condition
      double Y0 = gridY_[0];
      double Y1 = gridY_[1];
      double Y2 = gridY_[2];
      conv2 = -0.5 * (lambda_ * X + kappa2_ * Y0);
      rate = 0.25 * (alpha + X + Y0 - kappa1_ - kappa2_);
      double k0 =
          1.0 / delT + rate +
          (-conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0);
      double k1 = (conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1);
      double k2 = (-conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0);
      double d0 =
          (rate +
           (-conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0)) *
              inM[j][0] +
          (conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
          (-conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2] +
          midM[j][0] / delT;
      sol_c_[0] = sol_c_[0] - sol_l_[0] * k1 / k0;
      sol_u_[0] = sol_u_[0] - sol_l_[0] * k2 / k0;
      sol_d_[0] = sol_d_[0] - sol_l_[0] * d0 / k0;
      sol_l_[0] = 0.0;
      // k=M-1 upper boundary condition
      double Y3 = gridY_[Mxy_ - 3];
      Y2 = gridY_[Mxy_ - 2];
      Y1 = gridY_[Mxy_ - 1];
      conv2 = -0.5 * (lambda_ * X + kappa2_ * Y1);
      rate = 0.25 * (alpha + X + Y1 - kappa1_ - kappa2_);
      double kl3 = (conv2 * (Y1 - Y2) - difu2) / (Y2 - Y3) / (Y1 - Y3);
      double kl2 = (-conv2 * (Y1 - Y3) + difu2) / (Y2 - Y3) / (Y1 - Y2);
      double kl1 = 1.0 / delT + rate +
                   (conv2 * (2 * Y1 - Y2 - Y3) - difu2) / (Y1 - Y2) / (Y1 - Y3);
      double dl1 = (conv2 * (Y1 - Y2) - difu2) / (Y2 - Y3) / (Y1 - Y3) *
                       inM[j][Mxy_ - 3] +
                   (-conv2 * (Y1 - Y3) + difu2) / (Y2 - Y3) / (Y1 - Y2) *
                       inM[j][Mxy_ - 2] +
                   (rate + (conv2 * (2 * Y1 - Y2 - Y3) - difu2) / (Y1 - Y2) /
                               (Y1 - Y3)) *
                       inM[j][Mxy_ - 1] +
                   midM[j][Mxy_ - 1] / delT;
      sol_l_[Mxy_ - 3] = sol_l_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl3 / kl1;
      sol_c_[Mxy_ - 3] = sol_c_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * kl2 / kl1;
      sol_d_[Mxy_ - 3] = sol_d_[Mxy_ - 3] - sol_u_[Mxy_ - 3] * dl1 / kl1;
      sol_u_[Mxy_ - 3] = 0.0;
      velesquant::TriDiagonalSolve(Mxy_ - 2, sol_l_, sol_c_, sol_u_, sol_d_,
                                   sol_v_, sol_tempC_, sol_tempD_);
      for (int k = 1; k < Mxy_ - 1; k++)
        outM[j][k] = sol_v_[k - 1];
      // update k=0 lower boundary
      outM[j][0] = (d0 - k1 * outM[j][1] - k2 * outM[j][2]) / k0;
      // update k=M-1 upper boundary
      outM[j][Mxy_ - 1] =
          (dl1 - kl3 * outM[j][Mxy_ - 3] - kl2 * outM[j][Mxy_ - 2]) / kl1;
    }
  }
};
} // namespace velesquant
#endif