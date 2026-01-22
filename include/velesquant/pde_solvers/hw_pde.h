//		HWPDE.h

#ifndef HWPDE_H
#define HWPDE_H

#include <algorithm>

#include <cmath>
#include <list>
#include <map>
#include <memory>

#include <ql/errors.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <string>
#include <vector>
#include <velesquant/models/concepts.h>
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/cyclic_reduction.h>
#include <velesquant/volatility/lm.h>

namespace velesquant {

template <models::ShortRateModel ModelType> class HWPDE {
public:
  HWPDE(std::shared_ptr<ModelType> model, int grid_points = 512,
        double time_step = 0.001)
      : model_(std::move(model)), MR_(grid_points), delT_(time_step) {
    if (!model_)
      throw std::runtime_error("HWPDE: Model cannot be null");

    kappa_ = model_->getKappa();
    timeSigmas_ = model_->getTimeSigmas();
    sigmas_ = model_->getSigmas();
    timeDFs_ = model_->getTimeDFs();
    DFs_ = model_->getDFs();

    if (DFs_.size() > 1 && timeDFs_.size() > 1) {
      R0_ = -std::log(DFs_[1]) / timeDFs_[1];
    } else {
      R0_ = 0.03;
    }

    timeThetas_.push_back(0.5);
    if (!timeDFs_.empty()) {
      double dt = (timeDFs_.back() - timeThetas_[0]) / 20.0;
      if (dt > 1e-6) {
        for (double i = timeThetas_[0] + dt; i <= timeDFs_.back(); i += dt)
          timeThetas_.push_back(i);
      }
    }
    thetas_.resize(timeThetas_.size(), 0.001);

    buildGrid();

    if (!timeDFs_.empty()) {
      termStructureCalibrator();
    }

    sigmas0_ = sigmas_;
    kappa0_ = kappa_;
  }

  virtual ~HWPDE() = default;

  double getKappa() const { return kappa_; }
  double getR0() const { return R0_; }
  std::vector<double> getSigmas() const { return sigmas_; }
  std::vector<double> getThetas() const { return thetas_; }
  std::vector<double> getTimeSigmas() const { return timeSigmas_; }
  std::vector<double> getTimeThetas() const { return timeThetas_; }

  std::vector<double> getDFs(std::vector<double> &timePoints) {
    int N = timePoints.size();
    std::vector<int> iTs(N);
    std::vector<double> inV(MR_, 0.0), DFs(N, 0.0);
    for (int i = 0; i < N; i++)
      iTs[i] = int(timePoints[i] / delT_ + 1);

    inV[iR0_] = 2.0 / (gridR_[iR0_ + 1] - gridR_[iR0_ - 1]);

    for (int i = 0; i < N; i++) {
      int sT = 0;
      if (i > 0)
        sT = iTs[i - 1] - 1;
      for (int t = sT; t < iTs[i] - 1; t++)
        oneStepForward(t, inV);
      DFs[i] = trapezoidal(inV);
    }
    return DFs;
  }

  std::vector<double> simulationPDE(std::vector<double> times) const {
    int N = times.size();
    std::vector<double> path(N);
    int T = 0;
    double r = R0_;
    const double DT = 0.001;
    const double SDT = std::sqrt(DT);

    for (int i = 0; i <= int(times[N - 1] / DT + 0.5); i++) {
      if (T < N && i * DT <= times[T] && (i + 1) * DT > times[T]) {
        path[T] = r;
        T++;
      }
      double sigma = whichSigma(i * DT);
      double theta = whichTheta(i * DT);
      r += (theta - kappa_ * r) * DT + sigma * SDT * random_normal();
    }
    return path;
  }

  double pricingSwaption(double Expiry, double Tenor, double Strike,
                         double PayFrequency) {
    std::vector<double> payoff(MR_, -1.0);
    pricingCouponBondt(Expiry, Expiry + Tenor, -Strike, PayFrequency, payoff);
    for (int r = 0; r < MR_; r++)
      payoff[r] = std::max(0.0, payoff[r] + 1.0);
    if (Expiry > 0)
      discountBack(0, Expiry, payoff);
    return payoff[iR0_];
  }

  double pricingSwap(double Expiry, double Tenor, double Strike,
                     double PayFrequency) {
    std::vector<double> payoff(MR_, -1);
    pricingCouponBondt(Expiry, Expiry + Tenor, -Strike, PayFrequency, payoff);
    for (int r = 0; r < MR_; r++)
      payoff[r] += 1.0;
    if (Expiry > 0)
      discountBack(0, Expiry, payoff);
    return payoff[iR0_];
  }

  double pricingBermudan(double Expiry, double Tenor,
                         std::vector<double> Exercises, double Strike,
                         double PayFrequency) {
    std::vector<double> payoff(MR_, -1.0), swaption(MR_);
    Exercises.insert(Exercises.begin(), Expiry);
    Exercises.push_back(Expiry + Tenor);
    int Ne = Exercises.size();
    for (int e = Ne - 2; e >= 0; e--) {
      pricingCouponBondt(Exercises[e], Exercises[e + 1], -Strike, PayFrequency,
                         payoff);
      if (e == Ne - 2) {
        for (int r = 0; r < MR_; r++)
          swaption[r] = std::max(0.0, payoff[r] + 1.0);
      } else {
        discountBack(Exercises[e], Exercises[e + 1], swaption);
        for (int r = 0; r < MR_; r++)
          swaption[r] = std::max(swaption[r], payoff[r] + 1.0);
      }
    }
    if (Expiry > 0)
      discountBack(0, Expiry, swaption);
    return swaption[iR0_];
  }

  double pricingCallableSwap(double Expiry, double Tenor,
                             std::vector<double> Exercises, double Coupon,
                             double Strike, double PayFrequency,
                             OptionType type) {
    std::vector<double> payoff(MR_, -1), call_value(MR_);
    Exercises.push_back(Expiry + Tenor);
    int Ne = Exercises.size();
    int itype = (type == OptionType::Call) ? 1 : -1;
    for (int e = Ne - 2; e >= 0; e--) {
      pricingCouponBondt(Exercises[e], Exercises[e + 1], -Coupon, PayFrequency,
                         payoff);
      if (e == Ne - 2) {
        for (int r = 0; r < MR_; r++)
          call_value[r] = std::max(0.0, ((payoff[r] + 1.0) - Strike) * itype);
      } else {
        discountBack(Exercises[e], Exercises[e + 1], call_value);
        for (int r = 0; r < MR_; r++)
          call_value[r] =
              std::max(call_value[r], ((payoff[r] + 1.0) - Strike) * itype);
      }
    }
    if (Expiry > 0)
      discountBack(0, Exercises[0], call_value);
    return -itype * call_value[iR0_];
  }

  double pricingZBO(double Expiry, double Maturity, double Strike,
                    OptionType type) {
    std::vector<double> f(MR_, 1);
    discountBack(Expiry, Expiry + Maturity, f);
    for (int i = 0; i < MR_; i++)
      f[i] = (type == OptionType::Call) ? std::max(0.0, f[i] - Strike)
                                        : std::max(0.0, -f[i] + Strike);
    if (Expiry > 0)
      discountBack(0, Expiry, f);
    return f[iR0_];
  }

  double pricingZB(double Maturity) {
    std::vector<double> payoff(MR_, 1.0);
    discountBack(0, Maturity, payoff);
    return payoff[iR0_];
  }

  double pricingCBO(double Expiry, double Tenor, double Coupon, double Strike,
                    double PayFrequency, OptionType type) {
    std::vector<double> f(MR_, 1);
    pricingCouponBondt(Expiry, Expiry + Tenor, Coupon, PayFrequency, f);
    for (int i = 0; i < MR_; i++)
      f[i] = (type == OptionType::Call) ? std::max(0.0, f[i] - Strike)
                                        : std::max(0.0, Strike - f[i]);
    if (Expiry > 0)
      discountBack(0, Expiry, f);
    return f[iR0_];
  }

  double pricingCouponBond(double Expiry, double Tenor, double Coupon,
                           double PayFrequency) {
    std::vector<double> payoff(MR_, 1);
    pricingCouponBondt(Expiry, Expiry + Tenor, Coupon, PayFrequency, payoff);
    if (Expiry > 0)
      discountBack(0, Expiry, payoff);
    return payoff[iR0_];
  }

  void calibrator(std::vector<double> timeDFs, std::vector<double> DFs,
                  std::vector<defSwap> swapQuotes,
                  std::map<std::string, double> optimizer_params = {}) {
    QL_REQUIRE(!swapQuotes.empty(), "Swap quotes cannot be empty");
    QL_REQUIRE(!timeDFs.empty(), "Time DFs cannot be empty");

    timeDFs_ = timeDFs;
    DFs_ = DFs;
    quoteSwap_ = swapQuotes;

    std::vector<double> aux;
    for (const auto &sq : swapQuotes) {
      aux.push_back(sq.Expiry);
      aux.push_back(sq.Tenor);
      aux.push_back(sq.Expiry + sq.Tenor);
    }
    std::sort(aux.begin(), aux.end());

    timeThetas_.clear();
    timeThetas_.push_back(aux[0]);
    int ii = 1;
    int kk = 1;
    do {
      if (timeThetas_[ii - 1] < aux[kk]) {
        timeThetas_.push_back(aux[kk]);
        ii++;
      }
      kk++;
    } while ((timeThetas_[ii - 1] <=
              swapQuotes.back().Expiry + swapQuotes.back().Tenor) &&
             (aux.size() > kk));

    thetas_.assign(timeThetas_.size(), 0.005);

    calibrateKappa();

    timeSigmas_.clear();
    timeSigmas_.push_back(quoteSwap_[0].Expiry);
    double pos = timeSigmas_[0];
    for (size_t i = 1; i < quoteSwap_.size(); i++) {
      if (quoteSwap_[i].Expiry > pos) {
        timeSigmas_.push_back(quoteSwap_[i].Expiry);
        pos = quoteSwap_[i].Expiry;
      }
    }
    sigmas_ = sigmas0_;
    sigmas_.resize(timeSigmas_.size());

    int n = sigmas_.size() + 1;
    std::vector<double> x(n), lb(n), ub(n);
    for (int i = 0; i < n - 1; i++) {
      x[i] = sigmas_[i];
      lb[i] = 0;
      ub[i] = 1;
    }
    x[n - 1] = kappa0_;
    lb[n - 1] = 0;
    ub[n - 1] = 1.0;

    int m = quoteSwap_.size();
    marketSwaption_.resize(m);
    for (int i = 0; i < m; i++)
      marketSwaption_[i] = swaptionATM(
          quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].VolATM);

    std::vector<double> fvec(m);
    int info = 0;
    std::vector<double> diag(n), fjac(m * n), qtf(n), wa1(n), wa2(n), wa3(n),
        wa4(m);
    std::vector<int> ipvt(n);
    int nfev = 0;

    lmfcn fcn = [this, &lb, &ub](int m, int n, double *x, double *fvec,
                                 int *iflag) {
      this->objFcnCalibration(m, n, x, fvec, iflag, lb.data(), ub.data());
    };

    // Optimizer parameters with defaults
    double ftol = 1e-10;
    double xtol = 1e-10;
    double gtol = 1e-10;
    int maxfev = 10000;
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

    lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
          diag.data(), 1, 1.0, 0, &info, &nfev, fjac.data(), m, ipvt.data(),
          qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(), fcn);

    QL_ENSURE(info >= 1 && info <= 4,
              "Hull-White Model Calibration Fails (info=" +
                  std::to_string(info) + ")");

    for (int i = 0; i < n - 1; i++)
      sigmas_[i] = x[i];
    kappa_ = x[n - 1];

    model_->setVolatilityStructure(timeSigmas_, sigmas_);
    model_->setKappa(kappa_);
  }

  double getSwapRate(double Expiry, double Tenor, double PayFrequency) {
    double SwapLow = pricingSwap(Expiry, Tenor, 0.01, PayFrequency);
    double SwapHigh = pricingSwap(Expiry, Tenor, 0.02, PayFrequency);
    return 0.01 + SwapLow / (SwapLow - SwapHigh) * 0.01;
  }

  double getImpVolATM(double Expiry, double Tenor, double PayFrequency) {
    double SwapLow = pricingSwap(Expiry, Tenor, 0.00, PayFrequency);
    double SwapHigh = pricingSwap(Expiry, Tenor, 0.01, PayFrequency);
    double Level = (SwapLow - SwapHigh) / 0.01;
    double SwapRate = SwapLow / Level;
    double swaptionATM_val =
        pricingSwaption(Expiry, Tenor, SwapRate, PayFrequency);
    double Lo = 0.001;
    double Hi = 4.999;
    double Vol = 0.5 * (Lo + Hi);
    double swaptionVal =
        Level * SwapRate *
        (2.0 * velesquant::cdf_normal(0.5 * Vol * sqrt(Expiry)) - 1.0);
    int Niter = 0;
    do {
      Niter++;
      if (swaptionVal < swaptionATM_val)
        Lo = Vol;
      if (swaptionVal > swaptionATM_val)
        Hi = Vol;
      Vol = 0.5 * (Lo + Hi);
      swaptionVal =
          Level * SwapRate *
          (2.0 * velesquant::cdf_normal(0.5 * Vol * sqrt(Expiry)) - 1.0);
    } while (std::fabs(1.0 - swaptionVal / swaptionATM_val) >= 1.0E-6 &&
             Niter < 40);
    return Vol;
  }

private:
  std::shared_ptr<ModelType> model_;
  double R0_, kappa_, kappa0_;
  int MR_;
  double delT_;
  double vleft_, vright_;
  mutable std::vector<double> timeSigmas_, sigmas_, sigmas0_;
  mutable std::vector<double> timeThetas_, thetas_;
  mutable std::vector<double> a_, b_, c_, e_, f_, g_;
  mutable std::vector<double> solveL_, solveC_, solveU_, solveD_, solveTempC_,
      solveTempD_;
  std::vector<double> timeDFs_, DFs_;
  std::vector<defSwap> quoteSwap_;

  mutable std::vector<double> gridR_;
  int iR0_;
  mutable std::vector<double> marketSwaption_;

  // Grid Construction
  std::vector<double> grMesh(int Mv, double Rmax, double factor) {
    int N = Mv;
    std::vector<double> Mesh(Mv + 1);
    double d3 = Rmax / (N * factor);
    double start = asinh((-Rmax - R0_) / d3);
    double dx = 1.0 / N * (asinh((Rmax - R0_) / d3) - start);

    for (int i = 0; i <= N; i++)
      Mesh[i] = R0_ + d3 * sinh(start + i * dx);
    vleft_ = Mesh[0] - (Mesh[1] - Mesh[0]);
    vright_ = Mesh[Mv] + (Mesh[Mv] - Mesh[Mv - 1]);
    for (int i = 0; i < N; i++)
      if ((Mesh[i] < R0_) && (R0_ < Mesh[i + 1])) {
        Mesh[i + 1] = R0_;
        break;
      }
    return Mesh;
  }

  void buildGrid(double Rmax = 0.5, double factor = 0.7) {
    gridR_ = grMesh(MR_ - 1, Rmax, factor);
    a_.resize(MR_);
    b_.resize(MR_);
    c_.resize(MR_);
    e_.resize(MR_);
    f_.resize(MR_);
    g_.resize(MR_);
    solveL_.resize(MR_);
    solveC_.resize(MR_);
    solveU_.resize(MR_);
    solveD_.resize(MR_);
    solveTempC_.resize(MR_);
    solveTempD_.resize(MR_);

    for (int i = 0; i < MR_; i++) {
      double R = gridR_[i];
      double Ru = (i == MR_ - 1) ? vright_ : gridR_[i + 1];
      double Rl = (i == 0) ? vleft_ : gridR_[i - 1];
      a_[i] = 2.0 / (R - Rl) / (Ru - Rl);
      b_[i] = -2.0 / (R - Rl) / (Ru - R);
      c_[i] = 2.0 / (Ru - R) / (Ru - Rl);
      e_[i] = -(Ru - R) / (R - Rl) / (Ru - Rl);
      f_[i] = (Ru - 2 * R + Rl) / (R - Rl) / (Ru - R);
      g_[i] = (R - Rl) / (Ru - Rl) / (Ru - R);
    }
    for (int r = 0; r < MR_; r++) {
      if (gridR_[r] >= R0_) {
        iR0_ = r;
        break;
      }
    }
  }

  void oneStepBackward(const int t, std::vector<double> &inV) {
    double tm = (t + .5) * delT_;
    double sigma = whichSigma(tm);
    double f2 = sigma * sigma / 4;
    double theta = whichTheta(tm);
    double ndt = 1.0 / delT_;

    for (int r = 0; r < MR_; r++) {
      double f1 =
          ((r == 0) || (r == MR_ - 1)) ? 0 : 0.5 * (theta - kappa_ * gridR_[r]);
      solveL_[r] = f1 * e_[r] + f2 * a_[r];
      solveC_[r] = -ndt + f1 * f_[r] + f2 * b_[r] - gridR_[r] / 2;
      solveU_[r] = f1 * g_[r] + f2 * c_[r];
      if (r == 0) {
        solveU_[r] += solveL_[r];
        solveD_[r] =
            -((2 * ndt + solveC_[r]) * inV[r] + solveU_[r] * inV[r + 1]);
      } else if (r == MR_ - 1) {
        solveL_[r] += solveU_[r];
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (2.0 * ndt + solveC_[r]) * inV[r]);
      } else
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (2.0 * ndt + solveC_[r]) * inV[r] +
              solveU_[r] * inV[r + 1]);
    }
    solveL_[0] = 0.0;
    solveU_[MR_ - 1] = 0.0;
    velesquant::TriDiagonalSolve(MR_, solveL_, solveC_, solveU_, solveD_, inV,
                                 solveTempC_, solveTempD_);
  }

  void oneStepForward(const int T, std::vector<double> &inV) {
    double Tm = (T + .5) * delT_;
    double sigma = whichSigma(Tm);
    double f2 = -sigma * sigma / 4;
    double theta = whichTheta(Tm);
    double ndt = 1.0 / delT_;

    for (int r = 0; r < MR_; r++) {
      double f1 =
          ((r == 0) || (r == MR_ - 1)) ? 0 : 0.5 * (theta - kappa_ * gridR_[r]);
      double f3 = (gridR_[r] - kappa_) / 2;
      solveL_[r] = f1 * e_[r] + f2 * a_[r];
      solveC_[r] = ndt + f1 * f_[r] + f2 * b_[r] + f3;
      solveU_[r] = f1 * g_[r] + f2 * c_[r];
      if (r == 0) {
        solveU_[r] += solveL_[r];
        solveD_[r] =
            -((-2 * ndt + solveC_[r]) * inV[r] + solveU_[r] * inV[r + 1]);
      } else if (r == MR_ - 1) {
        solveL_[r] += solveU_[r];
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (-2.0 * ndt + solveC_[r]) * inV[r]);
      } else
        solveD_[r] =
            -(solveL_[r] * inV[r - 1] + (-2.0 * ndt + solveC_[r]) * inV[r] +
              solveU_[r] * inV[r + 1]);
    }
    solveL_[0] = 0.0;
    solveU_[MR_ - 1] = 0.0;
    velesquant::TriDiagonalSolve(MR_, solveL_, solveC_, solveU_, solveD_, inV,
                                 solveTempC_, solveTempD_);
  }

  double whichSigma(double t) const {
    int n = sigmas_.size();
    if (n == 0)
      return 0.003;
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
    int n = thetas_.size();
    if (n == 0)
      return 0.0;
    if (t < timeThetas_[0])
      return thetas_[0];
    if (t >= timeThetas_[n - 1])
      return thetas_[n - 1];
    for (int i = 1; i < n; i++) {
      if ((t >= timeThetas_[i - 1]) && (t < timeThetas_[i]))
        return thetas_[i];
    }
    return thetas_.back();
  }

  double getDFinterp(double t) {
    if (timeDFs_.empty())
      return 1.0;
    QuantLib::LogLinearInterpolation interp(timeDFs_.begin(), timeDFs_.end(),
                                            DFs_.begin());
    return interp(t, false);
  }

  double trapezoidal(std::vector<double> &inV) {
    double value = 0.0;
    for (int r = 1; r < MR_; r++)
      value += 0.5 * (gridR_[r] - gridR_[r - 1]) * (inV[r] + inV[r - 1]);
    return value;
  }

  void discountBack(double t0, double Tn, std::vector<double> &f) {
    int Nt = int(Tn / delT_ + 1);
    int N0 = int(t0 / delT_ + 1);
    for (int t = Nt - 2; t >= N0; t--)
      oneStepBackward(t, f);
  }

  void pricingCouponBondt(double t0, double Tn, double Coupon,
                          double PayFrequency, std::vector<double> &f) {
    int Nt = int(Tn / delT_ + 1);
    int iExpiry = int(t0 / delT_ + 1);
    int Ncoupon = int((Tn - t0) / PayFrequency + 0.5);
    for (int t = Nt - 2; t >= iExpiry; t--) {
      double couponTime = t0 + Ncoupon * PayFrequency;
      if (couponTime > t * delT_ && couponTime <= (t + 1) * delT_) {
        for (int r = 0; r < MR_; r++)
          f[r] += PayFrequency * Coupon;
        Ncoupon--;
      }
      oneStepBackward(t, f);
    }
  }

  void objFcnCalibration(int m, int n, double *x, double *fvec, int * /*iflag*/,
                         double *lb, double *ub) {
    double penalty = pen_fun(x, lb, ub, n);
    for (int i = 0; i < n - 1; i++)
      sigmas_[i] = x[i];
    kappa_ = x[n - 1];
    termStructureCalibrator();
    for (int j = 0; j < m; j++) {
      double modelSwaption =
          pricingSwaption(quoteSwap_[j].Expiry, quoteSwap_[j].Tenor,
                          quoteSwap_[j].SwapRate, quoteSwap_[j].Frequency);
      fvec[j] = modelSwaption - marketSwaption_[j] + penalty;
    }
  }

  void termStructureCalibrator() {
    int N = timeThetas_.size();
    std::vector<int> iTs(N);
    for (int i = 0; i < N; i++)
      iTs[i] = int(timeThetas_[i] / delT_ + 1);
    std::vector<double> DFinterp(N);
    for (int i = 0; i < N; i++)
      DFinterp[i] = getDFinterp(timeThetas_[i]);

    std::vector<double> inV(MR_, 0.0), lastV(MR_, 0);
    for (int r = 0; r < MR_; r++)
      if (gridR_[r] >= R0_) {
        inV[r] = 2.0 / (gridR_[r + 1] - gridR_[r - 1]);
        break;
      }
    for (int r = 0; r < MR_; r++)
      lastV[r] = inV[r];

    for (int i_p = 0; i_p < N; i_p++) {
      int sT = 0;
      if (i_p > 0)
        sT = iTs[i_p - 1] - 1;
      double df = 0.0;
      int niter = 0;
      do {
        niter++;
        for (int r = 0; r < MR_; r++)
          inV[r] = lastV[r];

        if (thetas_[i_p] == 0.0)
          thetas_[i_p] = 0.0005;
        thetas_[i_p] *= 1.001;
        for (int t = sT; t < iTs[i_p] - 1; t++)
          oneStepForward(t, inV);
        double dfUP = trapezoidal(inV);

        for (int r = 0; r < MR_; r++)
          inV[r] = lastV[r];

        thetas_[i_p] /= 1.001;
        for (int t = sT; t < iTs[i_p] - 1; t++)
          oneStepForward(t, inV);
        df = trapezoidal(inV);

        if (std::abs(dfUP - df) > 1e-12)
          thetas_[i_p] *= (1 + 0.001 * (DFinterp[i_p] - df) / (dfUP - df));

        thetas_[i_p] = std::max(-0.0050, thetas_[i_p]);
      } while (std::fabs(1.0 - df / DFinterp[i_p]) >= 1.0E-6 && niter < 10);

      for (int r = 0; r < MR_; r++)
        inV[r] = lastV[r];
      for (int t = sT; t < iTs[i_p] - 1; t++)
        oneStepForward(t, inV);
      for (int r = 0; r < MR_; r++)
        lastV[r] = inV[r];
    }
  }

  double pen_fun(double *x, double *lb, double *ub, int n) {
    double penalty = 0.0;
    double lambda = 1e8;
    for (int i = 0; i < n; i++)
      penalty +=
          ((x[i] <= lb[i]) || (x[i] >= ub[i])) ? lambda * std::pow(x[i], 2) : 0;
    return penalty;
  }

  double swaptionATM(double Expiry, double Tenor, double VolATM) {
    return (getDFinterp(Expiry) - getDFinterp(Expiry + Tenor)) *
           (2.0 * velesquant::cdf_normal(0.5 * VolATM * std::sqrt(Expiry)) -
            1.0);
  }

  double Bratio(double a, double Mi, double Tk, double Tj) {
    if (a == 0.0)
      return 100000;
    double t1 = (Tj - Mi);
    double t2 = (Tk - Mi);
    return (1.0 - std::exp(-a * t1)) / (1.0 - std::exp(-a * t2));
  }

  double objKappa(double a, const std::vector<double> &ivr,
                  const std::vector<double> &br, const std::vector<int> &ind) {
    int N = ivr.size();
    double sum = 0.0;
    for (int i = 1; i < N; i++) {
      double ex = quoteSwap_[ind[i]].Expiry;
      double t1 = quoteSwap_[ind[i]].Tenor + ex;
      double t2 = quoteSwap_[ind[i] + 1].Tenor + ex;
      double vs = std::fabs(br[i] * Bratio(a, ex, t1, t2)) - ivr[i];
      sum += vs * vs;
    }
    return sum;
  }

  void calibrateKappa() {
    int MAX_ITER = 1000;
    double X_TOL = 1e-5;
    double F_TOL = 1e-5;
    double a0 = kappa0_;
    double a1 = a0;
    double dF, F;
    double min = 1e-6;
    double max = 1;
    int N = quoteSwap_.size();
    if (N == 0)
      return;

    std::vector<double> iv_ratio, bond_ratio;
    std::vector<int> index;
    for (int i = 0; i < N - 1; i++) {
      if (quoteSwap_[i].Expiry == quoteSwap_[i + 1].Expiry) {
        double ex = quoteSwap_[i].Expiry;
        double t1 = quoteSwap_[i].Tenor + ex;
        double t2 = quoteSwap_[i + 1].Tenor + ex;
        index.push_back(i);
        iv_ratio.push_back(quoteSwap_[i + 1].VolATM / quoteSwap_[i].VolATM);
        bond_ratio.push_back((getDFinterp(ex) - getDFinterp(t1)) /
                             (getDFinterp(ex) - getDFinterp(t2)));
      }
    }
    if (index.empty())
      return;

    double delta;
    for (int i = 0; i < MAX_ITER; i++) {
      F = objKappa(a0, iv_ratio, bond_ratio, index);
      double Fu = objKappa(a0 + X_TOL, iv_ratio, bond_ratio, index);
      double Fd = objKappa(a0 - X_TOL, iv_ratio, bond_ratio, index);
      dF = 0.5 * (Fu - Fd) / X_TOL;
      if ((std::fabs(F) < F_TOL) || (std::fabs(dF) < F_TOL))
        break;

      a1 = a0;
      delta = F / dF;
      a1 -= delta;
      if (a1 <= min) {
        delta = 0.5 * (a0 - min);
        a1 = a0 - delta;
        if ((a1 == min) || (a1 == max))
          break;
      } else if (a1 >= max) {
        delta = 0.5 * (a0 - max);
        a1 = a0 - delta;
        if ((a1 == min) || (a1 == max))
          break;
      }
      if (delta > 0)
        max = a0;
      else
        min = a0;

      if (std::abs(a1 - a0) < X_TOL * std::fabs(a1))
        break;
      a0 = a1;
      if (std::fabs(max - min) < X_TOL)
        break;
    }
    kappa_ = a0;
  }
};

} // namespace velesquant

#endif