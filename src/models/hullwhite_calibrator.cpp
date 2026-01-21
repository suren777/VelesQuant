
#include <boost/math/distributions/normal.hpp>
#include <cmath>

#include <ql/qldefines.hpp>
#include <velesquant/constants.h>
#include <velesquant/models/hullwhite_calibrator.h>
#include <velesquant/volatility/lm.h>

namespace velesquant {
namespace models {

namespace {
double cdf_normal(double x) {
  boost::math::normal s;
  return boost::math::cdf(s, x);
}
} // namespace

HullWhiteCalibrator::HullWhiteCalibrator(
    std::shared_ptr<HullWhiteModel> model,
    std::shared_ptr<engines::HullWhiteAnalyticEngine<HullWhiteModel>> engine)
    : model_(std::move(model)), engine_(std::move(engine)) {}

void HullWhiteCalibrator::calibrate(const std::vector<defSwap> &swapQuotes,
                                    CalibrationTarget target) {
  quoteSwap_ = swapQuotes;

  // Initialize timeSigmas based on quotes
  std::vector<double> timeSigmas;
  if (!quoteSwap_.empty()) {
    timeSigmas.push_back(quoteSwap_[0].Expiry);
    double pos = timeSigmas[0];
    for (size_t i = 1; i < quoteSwap_.size(); i++) {
      if (quoteSwap_[i].Expiry > pos) {
        timeSigmas.push_back(quoteSwap_[i].Expiry);
        pos = quoteSwap_[i].Expiry;
      }
    }
  }

  // Resize sigmas to match structure
  std::vector<double> sigmas(timeSigmas.size(), 0.01);

  // Update model structure
  model_->setVolatilityStructure(timeSigmas, sigmas);

  // Optimization setup
  int n = static_cast<int>(sigmas.size()) + 1; // sigmas + kappa
  std::vector<double> x(n);
  std::vector<double> lb(n);
  std::vector<double> ub(n);

  for (int i = 0; i < n - 1; i++) {
    x[i] = sigmas[i];
    lb[i] = 1e-6;
    ub[i] = 1.0;
  }
  // Kappa
  x[n - 1] = model_->getKappa(); // Use current kappa as guess
  lb[n - 1] = 1e-6;
  ub[n - 1] = 1.0;

  int m = static_cast<int>(quoteSwap_.size());
  std::vector<double> fvec(m);

  // Minpack parameters
  double ftol = constants::calibration::TOLERANCE;
  double xtol = constants::calibration::TOLERANCE;
  double gtol = constants::calibration::TOLERANCE;
  int maxfev = constants::calibration::MAX_ITERATIONS;
  double epsfcn = constants::calibration::EPSILON;
  std::vector<double> diag(n);
  int mode = 1;
  double factor = 1;
  int nprint = 0;
  int info = 0;
  int nfev = 0;
  std::vector<double> fjac(m * n);
  int ldfjac = m;
  std::vector<int> ipvt(n);
  std::vector<double> qtf(n);
  std::vector<double> wa1(n);
  std::vector<double> wa2(n);
  std::vector<double> wa3(n);
  std::vector<double> wa4(m);

  // Prepare market data
  marketSwaption_.resize(m);
  if (target == CalibrationTarget::Price) {
    for (int i = 0; i < m; i++) {
      // Calculate ATM price if not provided directly or if Quote contains
      // VolATM Here we assume we want to match ATM price implied by VolATM
      // engine helper needed? No, use engine directly if we had a helper for
      // ATM price from Vol We can use Black formula. Actually `quoteSwap` has
      // SwapRate and VolATM. Price = Annuity * (2N(d1)-1) for ATM? Or just
      // match Vol? target == Price means we match Price. We calculate Market
      // Price from Market Vol.
      double dfT0 = model_->getDiscountFactor(quoteSwap_[i].Expiry);
      double dfTn =
          model_->getDiscountFactor(quoteSwap_[i].Expiry + quoteSwap_[i].Tenor);
      // This is simplified, usually exact swap schedule needed.
      // But existing code used get_swaptionATM logic:
      // (P(0,T0) - P(0,Tn)) * (2N - 1).
      // Let's replicate that logic.
      double priceScale = dfT0 - dfTn;
      // But wait, existing code calls 'swaptionATM' on *Model*? No, on itself
      // via `swaptionATM` helper. Helper used `getDF` (from model) and
      // `volATM`.

      marketSwaption_[i] =
          priceScale * (2.0 * cdf_normal(0.5 * quoteSwap_[i].VolATM *
                                         std::sqrt(quoteSwap_[i].Expiry)) -
                        1.0);
    }
  } else {
    for (int i = 0; i < m; i++)
      marketSwaption_[i] = quoteSwap_[i].VolATM;
  }

  double *lb_ptr = lb.data();
  double *ub_ptr = ub.data();
  bool isSwap = (target == CalibrationTarget::Price);

  lmfcn fcn = [this, lb_ptr, ub_ptr, isSwap](int m, int n, double *x,
                                             double *fvec, int *iflag) {
    if (isSwap)
      this->objFcnPrice(m, n, x, fvec, iflag, lb_ptr, ub_ptr);
    else
      this->objFcnIV(m, n, x, fvec, iflag, lb_ptr, ub_ptr);
  };

  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcn);

  // Update model with final results
  std::vector<double> finalSigmas(n - 1);
  for (int i = 0; i < n - 1; ++i)
    finalSigmas[i] = x[i];
  model_->setVolatilityStructure(timeSigmas, finalSigmas);
  model_->setKappa(x[n - 1]);
}

void HullWhiteCalibrator::objFcnPrice(int m, int n, double *x, double *fvec,
                                      [[maybe_unused]] int *iflag, double *lb,
                                      double *ub) {
  double penalty = pen_fun(x, lb, ub, n);

  std::vector<double> currentSigmas(n - 1);
  for (int i = 0; i < n - 1; i++) {
    currentSigmas[i] = x[i];
  }
  double currentKappa = x[n - 1];

  // Update model temporarily for pricing
  // Note: optimization modifies 'x', we must update model to reflect 'x'
  // before calling engine.
  model_->setVolatilityStructure(model_->getTimeSigmas(), currentSigmas);
  model_->setKappa(currentKappa);

  for (int i = 0; i < m; i++) {
    double modelSwaption =
        engine_->swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                          quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
    fvec[i] = modelSwaption - marketSwaption_[i] + penalty;
  }
}

void HullWhiteCalibrator::objFcnIV(int m, int n, double *x, double *fvec,
                                   [[maybe_unused]] int *iflag, double *lb,
                                   double *ub) {
  double penalty = pen_fun(x, lb, ub, n);

  std::vector<double> currentSigmas(n - 1);
  for (int i = 0; i < n - 1; i++) {
    currentSigmas[i] = x[i];
  }
  double currentKappa = x[n - 1];

  model_->setVolatilityStructure(model_->getTimeSigmas(), currentSigmas);
  model_->setKappa(currentKappa);

  for (int i = 0; i < m; i++) {
    double modelSwaption =
        engine_->swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                          quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
    // We need inverse black formula here
    // Helper swaptionIVblack logic:
    // Find vol such that BlackPrice(vol) == modelSwaption

    // ATM Pricing approximation inversion
    double dfT0 = model_->getDiscountFactor(quoteSwap_[i].Expiry);
    double dfTn =
        model_->getDiscountFactor(quoteSwap_[i].Expiry + quoteSwap_[i].Tenor);
    double annuityVal =
        dfT0 - dfTn; // Very rough approximation of annuity for ATM?
    // Actually (P0 - Pn) is the value of the "Swap" if coupon = swap rate.
    // ATM Swaption = Annuity * (2N(d1)-1).
    // And (P0 - Pn) = Annuity * SwapRate.
    // So ATM Swaption = (P0 - Pn) * (2N(...) - 1).

    // Using this relation:
    // modelSwaption / (P0 - Pn) = 2N(...) - 1
    // N(...) = 0.5 * (modelSwaption/(P0-Pn) + 1)
    // ... = normal_inverse(...)
    // vol = ...

    // Robustness check
    if (std::abs(annuityVal) < 1e-8) {
      fvec[i] = 100.0;
      continue;
    } // Error

    double ratio = modelSwaption / annuityVal;
    double arg = 0.5 * (ratio + 1.0);
    if (arg <= 0.0001 || arg >= 0.9999) {
      fvec[i] = 10.0; /* penalty */
      continue;
    }

    boost::math::normal s;
    double d1 = boost::math::quantile(s, arg);
    double impliedVol = d1 * 2 / std::sqrt(quoteSwap_[i].Expiry);

    fvec[i] = impliedVol - marketSwaption_[i] + penalty;
  }
}

double HullWhiteCalibrator::pen_fun(double *x, double *lb, double *ub, int n) {
  double penalty = 0.0;
  double lambda = constants::calibration::PENALTY_LAMBDA;
  for (int i = 0; i < n; i++)
    penalty += ((x[i] < lb[i]) ? lambda * std::pow(x[i], 2) : 0) +
               ((x[i] > ub[i]) ? lambda * std::pow(x[i], 2) : 0);
  return penalty;
}

// Stub for bootstrap if needed, or remove from header if not implementing yet.
void HullWhiteCalibrator::calibrateBootstrap(
    [[maybe_unused]] const std::vector<defSwap> &swapQuotes,
    [[maybe_unused]] CalibrationTarget target) {
  // Left empty for now, focus on main calibrator
}

double HullWhiteCalibrator::Bratio(double a, double M, double t1, double t2) {
  return (1.0 - std::exp(-a * (t2 - M))) / (1.0 - std::exp(-a * (t1 - M)));
}

// ... other helpers ...

} // namespace models
} // namespace velesquant
