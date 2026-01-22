//		ShortRate1FPDE.h

#ifndef ShortRate1FPDE_H
#define ShortRate1FPDE_H

#include <array>

#include <cmath>
#include <list>
#include <map>
#include <memory>

#include <ql/errors.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <string>
#include <vector>
#include <velesquant/models/concepts.h>
#include <velesquant/models/short_rate_1f_model.h>
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/cyclic_reduction.h>
#include <velesquant/volatility/lm.h>

namespace velesquant {

/**
 * @brief General 1-Factor Short Rate PDE Solver.
 *
 * This class implements a PDE solver for general one-factor short rate models
 * (e.g., Hull-White, Black-Karasinski) using a finite difference scheme.
 * It supports pricing of Zero Coupon Bonds, Coupon Bonds, Swaps, Swaptions,
 * Bermudan Swaptions, and other interest rate derivatives.
 *
 * The solver uses a time-dependent grid and handles calibration of the
 * term structure to match input discount factors.
 *
 * @tparam ModelType The underlying short rate model type (e.g.,
 * ShortRate1FModel).
 */
template <typename ModelType> class ShortRate1FPDE {
public:
  /**
   * @brief Construct a new ShortRate1FPDE solver.
   *
   * @param model Shared pointer to the underlying 1-factor short rate model.
   * @param grid_points Number of spatial grid points (default: 256).
   * @param time_step Time step size for the PDE solver (default: 0.0025).
   *
   * @throws std::runtime_error If model is null.
   */
public:
  ShortRate1FPDE(std::shared_ptr<ModelType> model, int grid_points = 256,
                 double time_step = 0.0025)
      : model_(std::move(model)), grid_points_(grid_points), delT_(time_step) {
    if (!model_)
      throw std::runtime_error("ShortRate1FPDE: Model cannot be null");

    // Initialize from model
    R0_ = model_->getR0();
    kappa_ = model_->getKappa();
    alpha_ = model_->getAlpha();
    beta_ = model_->getBeta();
    gamma_ = model_->getGamma();
    timeSigmas_ = model_->getTimeSigmas();
    sigmas_ = model_->getSigmas();
    timeDFs_ = model_->getTimeDFs();
    DFs_ = model_->getDFs();

    in_calibration_ = false;

    // Thetas Init
    timeThetas_.clear();
    timeThetas_.push_back(0.5);
    if (!timeDFs_.empty()) {
      double dt = timeDFs_.back() / 20.0;
      if (dt > 1e-6) {
        timeThetas_[0] = dt;
        for (double i = timeThetas_[0] + dt; i <= timeDFs_.back(); i += dt)
          timeThetas_.push_back(i);
      }
    }
    thetas_.assign(timeThetas_.size(), 0.001);

    buildGrid();

    if (!timeDFs_.empty()) {
      termStructureCalibrator();
    }
  }

  virtual ~ShortRate1FPDE() = default;

  // Pricing Methods
  /**
   * @brief Calculates the price of a Zero Coupon Bond.
   *
   * @param Maturity Time to maturity of the bond.
   * @return double Price of the Zero Coupon Bond.
   */
  double pricingZB(double Maturity) {
    std::vector<double> payoff(grid_points_, 1.0);
    discountBack(0, Maturity, payoff);
    return payoff[iR0_];
  }

  /**
   * @brief Calculates the price of a Zero Coupon Bond Option.
   *
   * @param Expiry Expiry time of the option.
   * @param Maturity Maturity of the underlying bond.
   * @param Strike Strike price.
   * @param type Option type (Call or Put).
   * @return double Price of the option.
   */
  double pricingZBO(double Expiry, double Maturity, double Strike,
                    velesquant::OptionType type) {
    std::vector<double> f(grid_points_, 1);
    discountBack(Expiry, Expiry + Maturity, f);
    for (int i = 0; i < grid_points_; i++)
      f[i] = (type == velesquant::OptionType::Call)
                 ? std::max(0.0, f[i] - Strike)
                 : std::max(0.0, -f[i] + Strike);
    if (Expiry > 0)
      discountBack(0, Expiry, f);
    return f[iR0_];
  }

  /**
   * @brief Calculates the price of a Coupon Bond.
   *
   * @param Expiry Expiry/Start time? (Often viewed as current time 0 for
   * pricing, but parameter name suggests otherwise). Actually interpreted as
   * start time of the bond or evaluation time? Implementation discounts back to
   * 0 if Expiry > 0.
   * @param Tenor Tenor of the bond.
   * @param Coupon Coupon rate.
   * @param PayFrequency Payment frequency of coupons.
   * @return double Price of the Coupon Bond.
   */
  double pricingCouponBond(double Expiry, double Tenor, double Coupon,
                           double PayFrequency) {
    std::vector<double> payoff(grid_points_, 1);
    pricingCouponBondt(Expiry, Expiry + Tenor, Coupon, PayFrequency, payoff);
    if (Expiry > 0)
      discountBack(0, Expiry, payoff);
    return payoff[iR0_];
  }

  /**
   * @brief Calculates the price of a Coupon Bond Option.
   *
   * @param Expiry Expiry time of the option.
   * @param Tenor Tenor of the underlying bond.
   * @param Coupon Coupon rate of the underlying bond.
   * @param Strike Strike price.
   * @param PayFrequency Payment frequency.
   * @param type Option type (Call or Put).
   * @return double Price of the option.
   */
  double pricingCBO(double Expiry, double Tenor, double Coupon, double Strike,
                    double PayFrequency, velesquant::OptionType type) {
    std::vector<double> f(grid_points_, 1);
    pricingCouponBondt(Expiry, Expiry + Tenor, Coupon, PayFrequency, f);
    for (int i = 0; i < grid_points_; i++)
      f[i] = (type == velesquant::OptionType::Call)
                 ? std::max(0.0, f[i] - Strike)
                 : std::max(0.0, Strike - f[i]);
    if (Expiry > 0)
      discountBack(0, Expiry, f);
    return f[iR0_];
  }

  /**
   * @brief Calculates the price of a generic Swap.
   *
   * @param Expiry Start time of the swap.
   * @param Tenor Length of the swap.
   * @param Strike Fixed rate (Strike) of the swap.
   * @param PayFrequency Frequency of payments.
   * @return double Value of the swap (payer/receiver depending on sign
   * convention, typically payer of fixed).
   */
  double pricingSwap(double Expiry, double Tenor, double Strike,
                     double PayFrequency) {
    std::vector<double> payoff(grid_points_, -1);
    pricingCouponBondt(Expiry, Expiry + Tenor, -Strike, PayFrequency, payoff);
    for (int r = 0; r < grid_points_; r++)
      payoff[r] += 1.0;
    if (Expiry > 0)
      discountBack(0, Expiry, payoff);
    return payoff[iR0_];
  }

  /**
   * @brief Calculates the price of a Callable Swap.
   *
   * @param Expiry Start time.
   * @param Tenor Swap tenor.
   * @param Exercises Vector of exercise times (excluding Expiry).
   * @param Coupon Fixed coupon rate.
   * @param Strike Strike price? (Usually 0 for par swaps, but effectively the
   * strike on the underlying value).
   * @param PayFrequency Payment frequency.
   * @param type Option type (Call means 'Right to Enter' usually).
   * @return double Price of the Callable Swap.
   */
  double pricingCallableSwap(double Expiry, double Tenor,
                             std::vector<double> Exercises, double Coupon,
                             double Strike, double PayFrequency,
                             velesquant::OptionType type) {
    std::vector<double> payoff(grid_points_, -1), call_value(grid_points_);
    Exercises.insert(Exercises.begin(), Expiry);
    Exercises.push_back(Expiry + Tenor);
    int Ne = Exercises.size();
    for (int e = Ne - 2; e >= 0; e--) {
      pricingCouponBondt(Exercises[e], Exercises[e + 1], -Coupon, PayFrequency,
                         payoff);
      if (e == Ne - 2) {
#pragma omp parallel for
        for (int r = 0; r < grid_points_; r++)
          call_value[r] = (type == velesquant::OptionType::Call)
                              ? std::max(0.0, (payoff[r] + 1.0) - Strike)
                              : std::max(0.0, -(payoff[r] + 1.0) + Strike);
      } else {
        discountBack(Exercises[e], Exercises[e + 1], call_value);
#pragma omp parallel for
        for (int r = 0; r < grid_points_; r++)
          call_value[r] = (type == velesquant::OptionType::Call)
                              ? std::max(0.0, (payoff[r] + 1.0) - Strike)
                              : std::max(0.0, -(payoff[r] + 1.0) + Strike);
      }
    }
    if (Expiry > 0)
      discountBack(0, Expiry, call_value);
    return call_value[iR0_];
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
    std::vector<double> payoff(grid_points_, -1.0);
    pricingCouponBondt(Expiry, Expiry + Tenor, -Strike, PayFrequency, payoff);
    for (int r = 0; r < grid_points_; r++)
      payoff[r] = std::max(0.0, payoff[r] + 1.0);
    if (Expiry > 0)
      discountBack(0, Expiry, payoff);
    return payoff[iR0_];
  }

  /**
   * @brief Calculates the price of a Bermudan Swaption.
   *
   * @param Expiry First exercise time (Expiry).
   * @param Tenor Tenor of the underlying swap.
   * @param Exercises Vector of subsequent exercise times.
   * @param Strike Strike rate.
   * @param PayFrequency Payment frequency.
   * @return double Price of the Bermudan Swaption.
   */
  double pricingBermudan(double Expiry, double Tenor,
                         std::vector<double> Exercises, double Strike,
                         double PayFrequency) {
    std::vector<double> payoff(grid_points_, -1.0), swaption(grid_points_);
    Exercises.insert(Exercises.begin(), Expiry);
    Exercises.push_back(Expiry + Tenor);
    int Ne = Exercises.size();
    for (int e = Ne - 2; e >= 0; e--) {
      pricingCouponBondt(Exercises[e], Exercises[e + 1], -Strike, PayFrequency,
                         payoff);
      if (e == Ne - 2) {
        for (int r = 0; r < grid_points_; r++)
          swaption[r] = std::max(0.0, payoff[r] + 1.0);
      } else {
        discountBack(Exercises[e], Exercises[e + 1], swaption);
        for (int r = 0; r < grid_points_; r++)
          swaption[r] = std::max(swaption[r], payoff[r] + 1.0);
      }
    }
    if (Expiry > 0)
      discountBack(0, Expiry, swaption);
    return swaption[iR0_];
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
      iTs[i] = int(timePoints[i] / delT_ + 1);
    std::vector<double> inV(grid_points_, 0.0), DFs(n, 0.0);
    for (int r = 1; r < grid_points_ - 1; r++)
      if (gridR_[r] >= R0_) {
        inV[r] = 2.0 / (gridR_[r + 1] - gridR_[r - 1]);
        break;
      }
    for (int i = 0; i < n; i++) {
      int sT = 0;
      if (i > 0)
        sT = iTs[i - 1] - 1;
      for (int t = sT; t < iTs[i] - 1; t++)
        oneStepForward(t, inV);
      DFs[i] = trapezoidal(inV);
      QL_ENSURE(DFs[i] <= 1.001, "pricingDF Functor Fails " << DFs[i]);
    }
    return DFs;
  }

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
    QL_REQUIRE(!swapQuotes.empty(), "Swap quotes cannot be empty");
    QL_REQUIRE(!timeDFs.empty(), "Time DFs cannot be empty");

    DFs_ = DFs;
    timeDFs_ = timeDFs;
    quoteSwap_ = swapQuotes;

    // Update model DFs
    model_->setDiscountFactors(timeDFs_,
                               DFs_); // Assuming model has this setter

    // Thetas Init
    timeThetas_.resize(1);
    timeThetas_[0] = 0.5;
    double dt = timeDFs_[timeDFs_.size() - 1] / 20;
    timeThetas_[0] = dt;
    for (double i = timeThetas_[0] + dt; i <= timeDFs_.back(); i += dt)
      timeThetas_.push_back(i);
    thetas_.assign(timeThetas_.size(), 0.001);

    // DF Interp
    DFinterp_.resize(timeThetas_.size());
    for (size_t i = 0; i < timeThetas_.size(); i++)
      DFinterp_[i] = getDFinterp(timeThetas_[i]);

    // Sigmas Init
    timeSigmas_.clear();
    timeSigmas_.push_back(quoteSwap_[0].Expiry);
    double pos = timeSigmas_[0];
    for (size_t i = 1; i < quoteSwap_.size(); i++) {
      if (quoteSwap_[i].Expiry > pos) {
        timeSigmas_.push_back(quoteSwap_[i].Expiry);
        pos = quoteSwap_[i].Expiry;
      }
    }
    sigmas_.assign(timeSigmas_.size(), 0.01);

    // Calibration Loop
    int ns = sigmas_.size();
    int n_params = ns + 1;
    std::vector<double> x_params(n_params);
    for (int i = 0; i < ns; i++)
      x_params[i] = sigmas_[i];

    // Kappa Calibration
    calibrateKappa();

    x_params[n_params - 1] = alpha_ / 1000; // Beta initial guess related
    int m_obs = quoteSwap_.size();
    QL_ENSURE(m_obs >= n_params, "too much freedom in Calibration");

    for (int i = 0; i < m_obs; i++)
      quoteSwap_[i].Value = swaptionATM(
          quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].VolATM);

    std::vector<double> fvec_vec(m_obs);
    int info = 0;
    std::vector<double> diag(n_params), fjac(m_obs * n_params), qtf(n_params),
        wa1(n_params), wa2(n_params), wa3(n_params), wa4(m_obs);
    std::vector<int> ipvt(n_params);
    int nfev = 0;

    lmfcn fcn = [this](int m, int n, double *x, double *fvec, int *iflag) {
      this->objFcnCalibrator(m, n, x, fvec, iflag);
    };

    // Optimizer parameters with defaults
    double ftol = 1e-6;
    double xtol = 1e-6;
    double gtol = 1e-6;
    int maxfev = 800;
    double epsfcn = 1e-7;

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

    in_calibration_ = true;
    lmdif(m_obs, n_params, x_params.data(), fvec_vec.data(), ftol, xtol, gtol,
          maxfev, epsfcn, diag.data(), 1, 1.0, 0, &info, &nfev, fjac.data(),
          m_obs, ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(),
          wa4.data(), fcn);
    in_calibration_ = false;

    QL_ENSURE(info >= 1 && info <= 4,
              "Model Calibration Fails (info=" + std::to_string(info) + ")");

    for (int i = 0; i < ns; i++)
      sigmas_[i] = std::fabs(x_params[i]);
    beta_ = (alpha_ - 1e-32) * std::sin(x_params[n_params - 1]) / 500;

    // Update Model
    model_->setVolatilityStructure(timeSigmas_, sigmas_);
    model_->setKappa(kappa_);
    model_->setBeta(beta_);
  }

  // Getters
  double getR0() { return R0_; }
  double getKappa() { return kappa_; }
  double getAlpha() { return alpha_; }
  double getBeta() { return beta_; }
  double getGamma() { return gamma_; }
  std::vector<double> getTimeSigmas() { return timeSigmas_; }
  std::vector<double> getSigmas() { return sigmas_; }
  std::vector<double> getTimeThetas() { return timeThetas_; }
  std::vector<double> getThetas() { return thetas_; }
  std::vector<double> getGrid() { return gridR_; }

  double getSwapRate(double Expiry, double Tenor, double PayFrequency) {
    double SwapLow = pricingSwap(Expiry, Tenor, 0.01, PayFrequency);
    double SwapHigh = pricingSwap(Expiry, Tenor, 0.02, PayFrequency);
    double SwapRate = 0.01 + SwapLow / (SwapLow - SwapHigh) * 0.01;
    return SwapRate;
  }

  double getImpVolATM(double Expiry, double Tenor, double PayFrequency) {
    double SwapLow = pricingSwap(Expiry, Tenor, 0.00, PayFrequency);
    double SwapHigh = pricingSwap(Expiry, Tenor, 0.01, PayFrequency);
    double Level = (SwapLow - SwapHigh) / 0.01;
    double SwapRate = SwapLow / Level;
    double swaption_atm_val =
        pricingSwaption(Expiry, Tenor, SwapRate, PayFrequency);

    double Lo = 0.001;
    double Hi = 4.999;
    double Vol = 0.5 * (Lo + Hi);
    double swaptionVal =
        Level * SwapRate *
        (2.0 * velesquant::cdf_normal(0.5 * Vol * std::sqrt(Expiry)) - 1.0);
    int Niter = 0;
    do {
      Niter++;
      if (swaptionVal < swaption_atm_val)
        Lo = Vol;
      else
        Hi = Vol;
      Vol = 0.5 * (Lo + Hi);
      swaptionVal =
          Level * SwapRate *
          (2.0 * velesquant::cdf_normal(0.5 * Vol * std::sqrt(Expiry)) - 1.0);
    } while (std::fabs(1.0 - swaptionVal / swaption_atm_val) >= 1.0E-6 &&
             Niter < 13);
    return Vol;
  }

  std::vector<double> simulationPDE(std::vector<double> times) const {

    return {};
  }

  std::vector<double> SwaptionDiagnostic(double Expiry, double Tenor,
                                         double Strike, double PayFrequency) {

    return {};
  }

private:
  std::shared_ptr<ModelType> model_;
  double R0_, kappa_, alpha_, beta_, gamma_;
  int grid_points_;
  double delT_;
  double vleft_, vright_;
  bool in_calibration_;
  int iR0_;
  mutable std::vector<double> timeSigmas_, sigmas_;
  mutable std::vector<double> timeThetas_, thetas_, DFinterp_;
  mutable std::vector<double> a_, b_, c_, e_, f_, g_;
  mutable std::vector<double> solveL_, solveC_, solveU_, solveD_, solveTempC_,
      solveTempD_;
  std::vector<double> DFs_, timeDFs_;
  std::vector<velesquant::defSwap> quoteSwap_;
  mutable std::vector<double> gridR_;

  // Grid Methods
  std::vector<double> grMesh(int Mv, double Rmax, double factor) {
    int N = Mv;
    std::vector<double> Mesh(Mv + 1);
    double d3 = Rmax / (N * factor);
    double start = velesquant::asinh((-Rmax - R0_) / d3);
    double dx = 1.0 / N * (velesquant::asinh((Rmax - R0_) / d3) - start);

    for (int i = 0; i <= N; i++)
      Mesh[i] = R0_ + d3 * std::sinh(start + i * dx);
    vleft_ = Mesh[0] - (Mesh[1] - Mesh[0]);
    vright_ = Mesh[Mv] + (Mesh[Mv] - Mesh[Mv - 1]);
    for (int i = 0; i < N; i++)
      if ((Mesh[i] < R0_) && (R0_ < Mesh[i + 1])) {
        Mesh[i + 1] = R0_;
        break;
      }
    return Mesh;
  }

  void buildGrid(double Rmax = 1, double factor = 0.7) {
    gridR_ = grMesh(grid_points_ - 1, Rmax, factor);
    a_.resize(grid_points_);
    b_.resize(grid_points_);
    c_.resize(grid_points_);
    e_.resize(grid_points_);
    e_.resize(grid_points_);
    f_.resize(grid_points_);
    g_.resize(grid_points_);
    solveL_.resize(grid_points_);
    solveC_.resize(grid_points_);
    solveU_.resize(grid_points_);
    solveD_.resize(grid_points_);
    solveTempC_.resize(grid_points_);
    solveTempD_.resize(grid_points_);

    for (int i = 0; i < grid_points_; i++) {
      double R = gridR_[i];
      double Ru = (i == grid_points_ - 1) ? vright_ : gridR_[i + 1];
      double Rl = (i == 0) ? vleft_ : gridR_[i - 1];
      a_[i] = 2.0 / (R - Rl) / (Ru - Rl);
      b_[i] = -2.0 / (R - Rl) / (Ru - R);
      c_[i] = 2.0 / (Ru - R) / (Ru - Rl);
      e_[i] = -(Ru - R) / (R - Rl) / (Ru - Rl);
      f_[i] = (Ru - 2 * R + Rl) / (R - Rl) / (Ru - R);
      g_[i] = (R - Rl) / (Ru - Rl) / (Ru - R);
    }
    for (int r = 0; r < grid_points_; r++) {
      if (gridR_[r] >= R0_) {
        iR0_ = r;
        break;
      }
    }
  }

  // PDE Solver Methods
  void discountBack(double t0, double Tn, std::vector<double> &f) {
    // const double delT = 10.0 / (4001 - 1.0); // using member delT_
    int Nt = int(Tn / delT_ + 1);
    int N0 = int(t0 / delT_ + 1);
    for (int t = Nt - 2; t >= N0; t--)
      oneStepBackward(t, f);
  }

  void oneStepBackward(const int t, std::vector<double> &inV) {
    double tm = (t + .5) * delT_;
    double sigma = whichSigma(tm);
    double sigma2 = sigma * sigma;
    double theta = whichTheta(tm);
    double ndt = 1.0 / delT_;

    for (int r = 0; r < grid_points_; r++) {
      double f1 = ((r == 0) || (r == grid_points_ - 1))
                      ? 0
                      : (theta - kappa_ * gridR_[r]);
      double f2 = sigma2 * std::pow(alpha_ + beta_ * gridR_[r], 2 * gamma_);
      double f3 = -gridR_[r];
      solveL_[r] = 0.5 * (f1 * e_[r] + f2 / 2 * a_[r]);
      solveC_[r] = -ndt + 0.5 * (f1 * f_[r] + f2 / 2 * b_[r] + f3);
      solveU_[r] = 0.5 * (f1 * g_[r] + f2 / 2 * c_[r]);
      if (r == 0) {
        solveU_[0] += solveL_[0];
        solveD_[0] = -((2 * ndt + solveC_[r]) * inV[0] + solveU_[r] * inV[1]);
      } else if (r == grid_points_ - 1) {
        solveL_[grid_points_ - 1] += solveU_[grid_points_ - 1];
        solveD_[grid_points_ - 1] =
            -(solveL_[grid_points_ - 1] * inV[grid_points_ - 2] +
              (2.0 * ndt + solveC_[r]) * inV[grid_points_ - 1]);
      } else
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (2.0 * ndt + solveC_[r]) * inV[r] +
              solveU_[r] * inV[r + 1]);
    }
    solveL_[0] = 0.0;
    solveU_[grid_points_ - 1] = 0.0;
    velesquant::TriDiagonalSolve(grid_points_, solveL_, solveC_, solveU_,
                                 solveD_, inV, solveTempC_, solveTempD_);
  }

  void oneStepForward(const int T, std::vector<double> &inV) {
    double Tm = (T + .5) * delT_;
    double sigma = whichSigma(Tm);
    double theta = whichTheta(Tm);
    double sigma2 = sigma * sigma;
    double p1 = 2 * sigma2 * beta_ * gamma_;
    double p2 = 0.5 * p1 * beta_ * (2 * gamma_ - 1);
    double ndt = 1.0 / delT_;

    for (int r = 0; r < grid_points_; r++) {
      double f1 = ((r == 0) || (r == grid_points_ - 1)) ? 0 : F1(r, p1, theta);
      double f2 = -F2(r, sigma2);
      double f3 = -F3(r, p2);
      solveL_[r] = 0.5 * (f1 * e_[r] + f2 / 2 * a_[r]);
      solveC_[r] = ndt + 0.5 * (f1 * f_[r] + f2 / 2 * b_[r] + f3);
      solveU_[r] = 0.5 * (f1 * g_[r] + f2 / 2 * c_[r]);
      if (r == 0) {
        solveU_[r] += solveL_[r];
        solveD_[r] =
            -((-2 * ndt + solveC_[r]) * inV[r] + solveU_[r] * inV[r + 1]);
      } else if (r == grid_points_ - 1) {
        solveL_[r] += solveU_[r];
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (-2.0 * ndt + solveC_[r]) * inV[r]);
      } else
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (-2.0 * ndt + solveC_[r]) * inV[r] +
              solveU_[r] * inV[r + 1]);
    }
    solveL_[0] = 0.0;
    solveU_[grid_points_ - 1] = 0.0;
    velesquant::TriDiagonalSolve(grid_points_, solveL_, solveC_, solveU_,
                                 solveD_, inV, solveTempC_, solveTempD_);
  }

  void pricingCouponBondt(double t0, double Tn, double Coupon,
                          double PayFrequency, std::vector<double> &f) {
    // const double delT = 10.0 / (4001 - 1.0); // using member delT_
    int Nt = int(Tn / delT_ + 1);
    int iExpiry = int(t0 / delT_ + 1);
    int Ncoupon = int((Tn - t0) / PayFrequency + 0.5);
    for (int t = Nt - 2; t >= iExpiry; t--) {
      double couponTime = t0 + Ncoupon * PayFrequency;
      if (couponTime > t * delT_ && couponTime <= (t + 1) * delT_) {
        for (int r = 0; r < grid_points_; r++)
          f[r] += PayFrequency * Coupon;
        Ncoupon--;
      }
      oneStepBackward(t, f);
    }
  }

  // Helpers
  double F1(int j, double p, double theta) {
    return theta - kappa_ * gridR_[j] -
           p * std::pow(alpha_ + beta_ * gridR_[j], 2 * gamma_ - 1);
  }
  double F2(int j, double sigma2) {
    return sigma2 * std::pow(alpha_ + beta_ * gridR_[j], 2 * gamma_);
  }
  double F3(int j, double p) {
    return p * std::pow(alpha_ + beta_ * gridR_[j], 2 * gamma_ - 2) + kappa_ -
           gridR_[j];
  }

  double trapezoidal(std::vector<double> &inV) {
    double value = 0.0;
    for (int r = 1; r < grid_points_; r++)
      value += 0.5 * (gridR_[r] - gridR_[r - 1]) * (inV[r] + inV[r - 1]);
    return value;
  }

  double getDFinterp(double t) {
    if (timeDFs_.empty())
      return 1.0;
    QuantLib::LogLinearInterpolation interp(timeDFs_.begin(), timeDFs_.end(),
                                            DFs_.begin());
    return interp(t, false);
  }

  double whichSigma(double t) const {
    if (timeSigmas_.empty())
      return 0.01;
    int n = sigmas_.size();
    if (t < timeSigmas_[0])
      return sigmas_[0];
    if (t >= timeSigmas_[n - 1])
      return sigmas_[n - 1];
    for (int i = 1; i < n; i++) {
      if ((t >= timeSigmas_[i - 1]) && (t < timeSigmas_[i]))
        return sigmas_[i];
    }
    return sigmas_.back();
  }

  double whichTheta(double t) const {
    if (timeThetas_.empty())
      return 0.0;
    int n = thetas_.size();
    if (t < timeThetas_[0])
      return thetas_[0];
    if (t >= timeThetas_[n - 1])
      return thetas_[n - 1];
    for (int i = 1; i < n; i++) {
      if (t < timeThetas_[i])
        return thetas_[i];
    }
    return thetas_.back();
  }

  double swaptionATM(double Expiry, double Tenor, double VolATM) {
    return (getDFinterp(Expiry) - getDFinterp(Expiry + Tenor)) *
           (2.0 * velesquant::cdf_normal(0.5 * VolATM * std::sqrt(Expiry)) -
            1.0);
  }

  void termStructureCalibrator() {
    // const double delT = 10.0 / (4001 - 1.0); // using member delT_
    int n = timeThetas_.size();
    std::vector<int> iTs(n);
    for (int i = 0; i < n; i++)
      iTs[i] = int(timeThetas_[i] / delT_ + 1);

    std::vector<double> inV(grid_points_, 0.0), lastV(grid_points_);
    for (int r = 1; r < grid_points_ - 1; r++)
      if (gridR_[r] >= R0_) {
        inV[r] = 2.0 / (gridR_[r + 1] - gridR_[r - 1]);
        break;
      }
    for (int r = 0; r < grid_points_; r++)
      lastV[r] = inV[r];

    for (int i = 0; i < n; i++) {
      int sT = 0;
      if (i > 0)
        sT = iTs[i - 1] - 1;
      double df = 0.0;
      int niter = 0;
      do {
        niter++;
        for (int r = 0; r < grid_points_; r++)
          inV[r] = lastV[r];

        if (thetas_[i] == 0.0)
          thetas_[i] = 0.0005;
        thetas_[i] *= 1.001;
        for (int t = sT; t < iTs[i] - 1; t++)
          oneStepForward(t, inV);
        double dfUP = trapezoidal(inV);

        for (int r = 0; r < grid_points_; r++)
          inV[r] = lastV[r];

        thetas_[i] /= 1.001;
        for (int t = sT; t < iTs[i] - 1; t++)
          oneStepForward(t, inV);
        df = trapezoidal(inV);

        if (std::abs(dfUP - df) > 1e-12)
          thetas_[i] *= (1 + 0.001 * (DFinterp_[i] - df) / (dfUP - df));
        thetas_[i] = std::max(-0.0050, thetas_[i]);
      } while (std::fabs(1.0 - df / DFinterp_[i]) >= 1.0E-6 && niter < 7);

      for (int r = 0; r < grid_points_; r++)
        inV[r] = lastV[r];
      for (int t = sT; t < iTs[i] - 1; t++)
        oneStepForward(t, inV);
      for (int r = 0; r < grid_points_; r++)
        lastV[r] = inV[r];
    }
  }

  void objFcnCalibrator(int m, int n, double *x, double *fvec,
                        int * /*iflag*/) {
    if (in_calibration_) {
      int ns = sigmas_.size();
      for (int i = 0; i < ns; i++)
        sigmas_[i] = std::fabs(x[i]);
      beta_ = (alpha_ - 1e-32) * std::sin(x[n - 1]) / 500;

      termStructureCalibrator();

      for (int i = 0; i < m; i++) {
        double model =
            pricingSwaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                            quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
        fvec[i] = 1 - quoteSwap_[i].Value / model;
      }
    }
  }

  void calibrateKappa() {
    kappa_ = 0.02;
    // ...
    int MAX_ITER = 1000;
    double X_TOL = 1e-5;
    double F_TOL = 1e-5;
    double a0 = 0.1;
    // ...

    if (quoteSwap_.empty())
      return;
    // ...

    // For now, assuming kappa_ is set by model or default
  }

}; // class

} // namespace velesquant

#endif