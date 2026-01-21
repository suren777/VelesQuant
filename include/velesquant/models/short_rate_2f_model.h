
#ifndef SHORTRATE2FMODEL_H
#define SHORTRATE2FMODEL_H

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace velesquant {
namespace models {

class ShortRate2FModel {
public:
  ShortRate2FModel(double kappa1, double kappa2, double lambda,
                   std::vector<double> timeSigma1s, std::vector<double> sigma1s,
                   std::vector<double> timeSigma2s, std::vector<double> sigma2s,
                   std::vector<double> timeAlphas, std::vector<double> alphas)
      : kappa1_(kappa1), kappa2_(kappa2), lambda_(lambda),
        timeSigma1s_(std::move(timeSigma1s)), sigma1s_(std::move(sigma1s)),
        timeSigma2s_(std::move(timeSigma2s)), sigma2s_(std::move(sigma2s)),
        timeAlphas_(std::move(timeAlphas)), alphas_(std::move(alphas)) {}

  virtual ~ShortRate2FModel() = default;

  // Getters
  double getKappa1() const { return kappa1_; }
  double getKappa2() const { return kappa2_; }
  double getLambda() const { return lambda_; }
  const std::vector<double> &getTimeSigma1s() const { return timeSigma1s_; }
  const std::vector<double> &getSigma1s() const { return sigma1s_; }
  const std::vector<double> &getTimeSigma2s() const { return timeSigma2s_; }
  const std::vector<double> &getSigma2s() const { return sigma2s_; }
  const std::vector<double> &getTimeAlphas() const { return timeAlphas_; }
  const std::vector<double> &getAlphas() const { return alphas_; }

  // Setters for calibration
  void setKappa1(double k) { kappa1_ = k; }
  void setKappa2(double k) { kappa2_ = k; }
  void setLambda(double l) { lambda_ = l; }
  void setSigma1Structure(const std::vector<double> &times,
                          const std::vector<double> &values) {
    timeSigma1s_ = times;
    sigma1s_ = values;
  }
  void setSigma2Structure(const std::vector<double> &times,
                          const std::vector<double> &values) {
    timeSigma2s_ = times;
    sigma2s_ = values;
  }
  void setAlphaStructure(const std::vector<double> &times,
                         const std::vector<double> &values) {
    timeAlphas_ = times;
    alphas_ = values;
  }

private:
  double kappa1_;
  double kappa2_;
  double lambda_;
  std::vector<double> timeSigma1s_;
  std::vector<double> sigma1s_;
  std::vector<double> timeSigma2s_;
  std::vector<double> sigma2s_;
  std::vector<double> timeAlphas_;
  std::vector<double> alphas_;
};

} // namespace models
} // namespace velesquant

#endif // SHORTRATE2FMODEL_H
