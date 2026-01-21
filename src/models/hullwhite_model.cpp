#include <algorithm>
#include <cmath>
#include <velesquant/constants.h>
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/utility.h>

namespace velesquant {
namespace models {

HullWhiteModel::HullWhiteModel(double kappa,
                               const std::vector<double> &timeSigmas,
                               const std::vector<double> &sigmas,
                               const std::vector<double> &timeDFs,
                               const std::vector<double> &DFs)
    : kappa_(kappa), timeSigmas_(timeSigmas), sigmas_(sigmas),
      timeDFs_(timeDFs), DFs_(DFs) {}

double HullWhiteModel::getDiscountFactor(double t) const {
  return interpolation("CubicNaturalSpline", timeDFs_, DFs_, t);
}

double HullWhiteModel::getSigma(double t) const {
  if (t < timeSigmas_.front())
    return sigmas_.front();

  auto it = std::upper_bound(timeSigmas_.begin(), timeSigmas_.end(), t);

  if (it == timeSigmas_.end())
    return sigmas_.back();

  size_t i = std::distance(timeSigmas_.begin(), it);
  return sigmas_[i];
}

double HullWhiteModel::B(double t, double T) const {
  return (1.0 - std::exp(-kappa_ * (T - t))) / kappa_;
}

double HullWhiteModel::A(double t, double T) const {
  double P0t = getDiscountFactor(t);
  double P0T = getDiscountFactor(T);
  double f0t = -std::log(getDiscountFactor(t + 0.0001) / P0t) /
               0.0001; // Approximate f(0,t)
  // Actually, P(t,T) = A * exp(-B * r(t)).
  // Standard formula: A(t,T) = P(0,T)/P(0,t) * exp( B(t,T)*f(0,t) - ... ) is
  // for r(t) definition. Let's stick to what allows us to price options. For
  // bond options, we need P(0,T) and integrated variance.
  return P0T / P0t; // Simplified; actual affine term depends on calibration.
                    // Usually purely needed for simulating P(t,T) given r(t).
}

double HullWhiteModel::variance([[maybe_unused]] double t, double T) const {
  // Calculates Variance of integral_{t}^{T} sigma(s) * exp(-kappa(T-s)) ds ?
  // Or just variance of r(T) given r(t)?
  // The existing totalVariance(T0) in HW.cpp seems to be Var(r(T0) | F0).
  // Let's port that logic.

  // We assume current time is 0 for this calculation as per original
  // implementation logic structure If we want Var(r(T) | F0), here it is:
  double T0 = T; // using T as the point of variance calculation
  double var = 0.0;
  int nP = static_cast<int>(timeSigmas_.size());
  for (int n = 0; n < nP; n++) {
    double ti = (n > 0) ? timeSigmas_[n - 1] : 0.0;
    double te = timeSigmas_[n];

    double upper = (T0 > te) ? te : T0;

    if (upper > ti) {
      var += sigmas_[n] * sigmas_[n] *
             (std::exp(2 * kappa_ * upper) - std::exp(2 * kappa_ * ti)) /
             (2 * kappa_);
    }

    if (T0 <= te)
      break;
  }

  // Handle tail
  if (T0 > timeSigmas_.back()) {
    var += sigmas_.back() * sigmas_.back() *
           (std::exp(2 * kappa_ * T0) -
            std::exp(2 * kappa_ * timeSigmas_.back())) /
           (2 * kappa_);
  }

  return std::exp(-2 * kappa_ * T0) * var;
}

std::vector<double>
HullWhiteModel::simulation(std::vector<double> times) const {
  std::vector<double> r_path;
  // r(0) = f(0,0) approx
  double r0 = -std::log(getDiscountFactor(0.0001)) / 0.0001;

  double x = 0.0;
  double t_prev = 0.0;
  // kappa_ is member
  double kappa = kappa_;

  for (double t : times) {
    if (t <= 1e-6) {
      r_path.push_back(r0);
      t_prev = t;
      continue;
    }

    double dt = t - t_prev;
    if (dt < 1e-8) {
      // Duplicate time handling
      if (r_path.empty())
        r_path.push_back(r0);
      else
        r_path.push_back(r_path.back());
      continue;
    }

    double sigma = getSigma(t_prev);

    double mean = x * std::exp(-kappa * dt);
    double var =
        (sigma * sigma) / (2 * kappa) * (1.0 - std::exp(-2 * kappa * dt));
    double std_dev = std::sqrt(var);

    x = mean + std_dev * random_normal();

    double P_t = getDiscountFactor(t);
    double P_t_dt = getDiscountFactor(t + 0.0001);
    double f_0_t = -std::log(P_t_dt / P_t) / 0.0001;

    double alpha = f_0_t + (sigma * sigma) / (2 * kappa * kappa) *
                               std::pow(1.0 - std::exp(-kappa * t), 2.0);

    r_path.push_back(x + alpha);
    t_prev = t;
  }
  return r_path;
}

} // namespace models
} // namespace velesquant
