#ifndef VELESQUANT_MODELS_SHORTRATE_1F_MODEL_H
#define VELESQUANT_MODELS_SHORTRATE_1F_MODEL_H

#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

namespace velesquant {
namespace models {

class ShortRate1FModel {
public:
  ShortRate1FModel(double R0, double kappa, double alpha, double beta,
                   double gamma, std::vector<double> timeSigmas,
                   std::vector<double> sigmas)
      : R0_(R0), kappa_(kappa), alpha_(alpha), beta_(beta), gamma_(gamma),
        timeSigmas_(std::move(timeSigmas)), sigmas_(std::move(sigmas)) {}

  virtual ~ShortRate1FModel() = default;

  // Getters
  double getR0() const { return R0_; }
  double getKappa() const { return kappa_; }
  double getAlpha() const { return alpha_; }
  double getBeta() const { return beta_; }
  double getGamma() const { return gamma_; }
  const std::vector<double> &getTimeSigmas() const { return timeSigmas_; }
  const std::vector<double> &getSigmas() const { return sigmas_; }
  const std::vector<double> &getTimeDFs() const { return timeDFs_; }
  const std::vector<double> &getDFs() const { return DFs_; }

  // Setters
  void setKappa(double kappa) { kappa_ = kappa; }
  void setAlpha(double alpha) { alpha_ = alpha; }
  void setBeta(double beta) { beta_ = beta; }
  void setGamma(double gamma) { gamma_ = gamma; }
  void setVolatilityStructure(const std::vector<double> &timeSigmas,
                              const std::vector<double> &sigmas) {
    timeSigmas_ = timeSigmas;
    sigmas_ = sigmas;
  }
  void setDiscountFactors(const std::vector<double> &timeDFs,
                          const std::vector<double> &DFs) {
    timeDFs_ = timeDFs;
    DFs_ = DFs;
  }
  void setR0(double R0) { R0_ = R0; }

  // ShortRateModel Concept Requirements

  double getDiscountFactor(double t) const {

    return std::exp(-R0_ * t); // Simplified
  }

  double variance(double t, double T) const {
    // Simplified variance placeholder
    return 0.0;
  }

private:
  double R0_;
  double kappa_;
  double alpha_;
  double beta_;
  double gamma_;
  std::vector<double> timeSigmas_;
  std::vector<double> sigmas_;
  std::vector<double> timeDFs_;
  std::vector<double> DFs_;
};

} // namespace models
} // namespace velesquant

#endif // VELESQUANT_MODELS_SHORTRATE_1F_MODEL_H
