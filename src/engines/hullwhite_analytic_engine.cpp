#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <velesquant/engines/hullwhite_analytic_engine.h>

namespace velesquant {
namespace engines {

namespace {
double cdf_normal(double x) {
  boost::math::normal s;
  return boost::math::cdf(s, x);
}
} // namespace

HullWhiteAnalyticEngine::HullWhiteAnalyticEngine(
    std::shared_ptr<models::HullWhiteModel> model)
    : model_(std::move(model)) {}

double HullWhiteAnalyticEngine::optionBond(double expiry, double maturity,
                                           double strike,
                                           OptionType type) const {
  double v0 = model_->variance(0, expiry);
  double kappa = model_->getKappa();
  double impVol = (1.0 - std::exp(-kappa * (maturity - expiry))) / kappa *
                  std::sqrt(v0 / expiry);

  return formulaBlack(model_->getDiscountFactor(expiry),
                      model_->getDiscountFactor(maturity), strike, impVol,
                      expiry, type);
}

double HullWhiteAnalyticEngine::swaption(double expiry, double tenor,
                                         double strike,
                                         double payFrequency) const {
  double CP = criticalPoint(expiry, tenor, strike, payFrequency);
  double v0 = model_->variance(0, expiry);
  double dfT0 = model_->getDiscountFactor(expiry);

  int nC = static_cast<int>(tenor / payFrequency + 0.5);
  double swaptionV = 0.0;
  double kappa = model_->getKappa();

  for (int i = 1; i <= nC; i++) {
    double Ti = expiry + i * payFrequency;
    double dfTi = model_->getDiscountFactor(Ti);
    double Gi = (1.0 - std::exp(-kappa * (Ti - expiry))) / kappa;
    double Ki = dfTi / dfT0 * std::exp(-CP * Gi - 0.5 * v0 * Gi * Gi);
    double impVol = (1.0 - std::exp(-kappa * (Ti - expiry))) / kappa *
                    std::sqrt(v0 / expiry);

    double Puti = formulaBlack(dfT0, dfTi, Ki, impVol, expiry, OptionType::Put);
    swaptionV += strike * payFrequency * Puti;

    if (i == nC)
      swaptionV += Puti;
  }
  return swaptionV;
}

double HullWhiteAnalyticEngine::criticalPoint(double expiry, double tenor,
                                              double strike,
                                              double payFrequency) const {
  double totalVar = model_->variance(0, expiry);
  double lowerBound = -0.5;
  double upperBound = 0.5;
  double shortRate = 0.5 * (lowerBound + upperBound);

  double equalityValue =
      equalityCP(shortRate, totalVar, expiry, tenor, strike, payFrequency);
  int Niter = 0;

  do {
    Niter++;
    if (equalityValue < 0.0)
      lowerBound = shortRate;
    if (equalityValue > 0.0)
      upperBound = shortRate;

    shortRate = 0.5 * (lowerBound + upperBound);
    equalityValue =
        equalityCP(shortRate, totalVar, expiry, tenor, strike, payFrequency);

  } while (std::abs(equalityValue) >= 1.0E-10 && Niter < 20 &&
           (upperBound - lowerBound) >= 10e-6);

  return shortRate;
}

double HullWhiteAnalyticEngine::equalityCP(double x, double v0, double expiry,
                                           double tenor, double strike,
                                           double payFrequency) const {
  double dfT0 = model_->getDiscountFactor(expiry);
  int nC = static_cast<int>(tenor / payFrequency + 0.5);
  double CP = 1.0;
  double kappa = model_->getKappa();

  for (int i = 1; i <= nC; i++) {
    double Ti = expiry + i * payFrequency;
    double Gi = (1.0 - std::exp(-kappa * (Ti - expiry))) / kappa;
    double Ki = model_->getDiscountFactor(Ti) / dfT0 *
                std::exp(-x * Gi - 0.5 * v0 * Gi * Gi);

    CP -= strike * payFrequency * Ki;

    if (i == nC)
      CP -= Ki;
  }
  return CP;
}

double HullWhiteAnalyticEngine::formulaBlack(double dfT0, double dfTN,
                                             double strike, double impVol,
                                             double T0, OptionType type) const {
  double sd = impVol * std::sqrt(T0);
  double d1 = std::log(dfTN / (dfT0 * strike)) / sd + 0.5 * sd;
  double d2 = d1 - sd;

  double price;
  if (type == OptionType::Call)
    price = dfTN * cdf_normal(d1) - (dfT0 * strike) * cdf_normal(d2);
  else
    price = (dfT0 * strike) * cdf_normal(-d2) - dfTN * cdf_normal(-d1);

  return price;
}

} // namespace engines
} // namespace velesquant
