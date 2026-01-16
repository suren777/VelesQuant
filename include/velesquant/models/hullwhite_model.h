#ifndef VELESQUANT_MODELS_HULLWHITE_MODEL_H
#define VELESQUANT_MODELS_HULLWHITE_MODEL_H

#include <cmath>
#include <vector>
#include <velesquant/models/utility.h>
#include <velesquant/numerics/interpolation.h>

namespace velesquant {
namespace models {

/**
 * @class HullWhiteModel
 * @brief Represents the Hull-White 1-Factor Model parameters and core dynamics.
 *
 * This class holds the model parameters (Mean reversion, Volatility structure)
 * and the term structure of interest rates. It provides methods to access
 * these parameters and calculate model-specific quantities like A(t,T) and
 * B(t,T).
 */
class HullWhiteModel {
public:
  HullWhiteModel(double kappa, const std::vector<double> &timeSigmas,
                 const std::vector<double> &sigmas,
                 const std::vector<double> &timeDFs,
                 const std::vector<double> &DFs);

  virtual ~HullWhiteModel() = default;

  // Getters
  [[nodiscard]] double getKappa() const { return kappa_; }
  [[nodiscard]] const std::vector<double> &getTimeSigmas() const {
    return timeSigmas_;
  }
  [[nodiscard]] const std::vector<double> &getSigmas() const { return sigmas_; }
  [[nodiscard]] const std::vector<double> &getDFs() const { return DFs_; }
  [[nodiscard]] const std::vector<double> &getTimeDFs() const {
    return timeDFs_;
  }

  // Setters (useful for calibration)
  void setKappa(double kappa) { kappa_ = kappa; }
  void setSigmas(const std::vector<double> &sigmas) { sigmas_ = sigmas; }
  void setVolatilityStructure(const std::vector<double> &timeSigmas,
                              const std::vector<double> &sigmas) {
    timeSigmas_ = timeSigmas;
    sigmas_ = sigmas;
  }

  // Core Dynamics
  [[nodiscard]] double getDiscountFactor(double t) const;
  [[nodiscard]] double getSigma(double t) const;

  // Bond Price P(t,T) decomposition: P(t,T) = A(t,T) * exp(-B(t,T) * r(t))
  // Note: ZC usually refers to P(0, T) or P(t, T).
  // Here we provide the components A and B used in affine term structure
  // models.
  [[nodiscard]] double B(double t, double T) const;
  [[nodiscard]] double A(double t, double T) const;

  // Variance of the accumulated short rate, used in option pricing
  // Variance of the accumulated short rate, used in option pricing
  [[nodiscard]] double variance(double t, double T) const;

  // Simulation
  [[nodiscard]] std::vector<double> simulation(std::vector<double> times) const;

private:
  double kappa_;
  std::vector<double> timeSigmas_;
  std::vector<double> sigmas_;
  std::vector<double> timeDFs_;
  std::vector<double> DFs_;
};

} // namespace models
} // namespace velesquant

#endif // VELESQUANT_MODELS_HULLWHITE_MODEL_H
