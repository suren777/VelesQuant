#include <cmath>
#include <memory>
#include <vector>

#include <velesquant/engines/hullwhite_analytic_engine.h>
#include <velesquant/models/hullwhite_calibrator.h>
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/hw.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>

namespace velesquant {

HullWhite::HullWhite(double kappa, std::vector<double> timeSigmas,
                     std::vector<double> sigmas, std::vector<double> timeDFs,
                     std::vector<double> DFs)
    : kappa_(kappa), timeDFs_(timeDFs), DFs_(DFs), timeSigmas_(timeSigmas),
      sigmas_(sigmas) {
  model_ = std::make_shared<models::HullWhiteModel>(kappa, timeSigmas, sigmas,
                                                    timeDFs, DFs);
  engine_ = std::make_shared<engines::HullWhiteAnalyticEngine>(model_);
  calibrator_ = std::make_shared<models::HullWhiteCalibrator>(model_, engine_);
}

double HullWhite::optionBond(double expiry, double maturity, double strike,
                             OptionType type) {
  return engine_->optionBond(expiry, maturity, strike, type);
}

double HullWhite::swaption(double expiry, double tenor, double strike,
                           double payFrequency) {
  return engine_->swaption(expiry, tenor, strike, payFrequency);
}

double HullWhite::ZC(double expiry) {
  return model_->getDiscountFactor(expiry);
}

void HullWhite::calibrator(const std::vector<defSwap> &swapQuotes,
                           CalibrationTarget target) {
  calibrator_->calibrate(swapQuotes, target);
  kappa_ = model_->getKappa();
  sigmas_ = model_->getSigmas();
  timeSigmas_ = model_->getTimeSigmas();
}

// Deprecated / Unimplemented for now
double HullWhite::BlackStrike([[maybe_unused]] double T0,
                              [[maybe_unused]] double TN,
                              [[maybe_unused]] double impVol) {
  return 0.0;
}
double HullWhite::BlackStrikePlain([[maybe_unused]] double T0,
                                   [[maybe_unused]] double TN) {
  return 0.0;
}
double HullWhite::swaptionIVblackPub([[maybe_unused]] double expiry,
                                     [[maybe_unused]] double tenor,
                                     [[maybe_unused]] double swap_price) {
  return 0.0;
}
double HullWhite::getSwapRate([[maybe_unused]] double expiry,
                              [[maybe_unused]] double tenor,
                              [[maybe_unused]] double payFrequency) {
  return 0.0;
}

std::vector<double> HullWhite::simulation(std::vector<double> times) const {
  std::vector<double> r_path;
  // r(0) = f(0,0) approx
  double r0 = -std::log(model_->getDiscountFactor(0.0001)) / 0.0001;

  double x = 0.0;
  double t_prev = 0.0;
  double kappa = model_->getKappa();

  for (double t : times) {
    if (t <= 1e-6) {
      r_path.push_back(r0);
      t_prev = t;
      continue;
    }

    double dt = t - t_prev;
    if (dt < 1e-8) {
      // Duplicate time handling: just push last state or r0
      if (r_path.empty())
        r_path.push_back(r0);
      else
        r_path.push_back(r_path.back());
      continue;
    }

    double sigma = model_->getSigma(t_prev);

    double mean = x * std::exp(-kappa * dt);
    double var =
        (sigma * sigma) / (2 * kappa) * (1.0 - std::exp(-2 * kappa * dt));
    double std_dev = std::sqrt(var);

    x = mean + std_dev * random_normal();

    double P_t = model_->getDiscountFactor(t);
    double P_t_dt = model_->getDiscountFactor(t + 0.0001);
    double f_0_t = -std::log(P_t_dt / P_t) / 0.0001;

    double alpha = f_0_t + (sigma * sigma) / (2 * kappa * kappa) *
                               std::pow(1.0 - std::exp(-kappa * t), 2.0);

    r_path.push_back(x + alpha);
    t_prev = t;
  }
  return r_path;
}

} // namespace velesquant