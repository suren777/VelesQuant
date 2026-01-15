// schobzhu

//		schobzhu.h

#ifndef SCHOBZHU_H
#define SCHOBZHU_H

#include <complex>
#include <vector>
#include <velesquant/models/utility.h>

namespace velesquant {

typedef std::complex<double> Cdoub;
typedef std::vector<double> Vdoub;

/**
 * @class SchobelZhu
 * @brief Schobel-Zhu Stochastic Volatility Model.
 *
 * Supports pricing of options via Schobel-Zhu integration.
 * Supports calibration to market prices or volatilities.
 */
class SchobelZhu {
public:
  SchobelZhu(double spot, double var0, double kappa, double theta, double xi,
           double rho); // parameters constructor
  ~SchobelZhu() {};

  // Schobel&Zhu integrand_1 and integrand_2

  [[nodiscard]] double SchobelPrice(double maturity, double forward,
                                    double strike) const;
  [[nodiscard]] Vdoub simulation(const Vdoub &times,
                                         const Vdoub &forwards) const;

  [[nodiscard]] double getParameterVar0() const { return sigma0_; };
  [[nodiscard]] double getParameterKappa() const { return kappa_; };
  [[nodiscard]] double getParameterTheta() const { return theta_; };
  [[nodiscard]] double getParameterXi() const { return xi_; };
  [[nodiscard]] double getParameterRho() const { return rho_; };

  void setParameterVar0(double var0) { sigma0_ = var0; };
  void setParameterKappa(double kappa) { kappa_ = fabs(kappa); };
  void setParameterTheta(double theta) { theta_ = fabs(theta); };
  void setParameterXi(double xi) { xi_ = fabs(xi); };
  void setParameterRho(double rho) { rho_ = sin(rho); };

  /**
   * @brief Calibrates model parameters to market quotes.
   *
   * @param maturitys Vector of maturities.
   * @param forwards Vector of forward prices.
   * @param strikes Vector of strikes.
   * @param marketQuotes Vector of market quotes (Price or Volatility).
   * @param target Calibration target type.
   */
  void calibrator(const Vdoub &maturitys, const Vdoub &forwards,
                  const Vdoub &strikes, const Vdoub &marketQuotes,
                  CalibrationTarget target = CalibrationTarget::Price);

private:
  double s0_;     // Spot price
  double sigma0_; // Initial volatility
  double kappa_;  // Mean reversion rate for volatility
  double theta_;  // Long run volatility
  double xi_;     // Volatility of Volatility
  double rho_;    // Price-volatility correlation
  Vdoub maturitys_, forwards_, strikes_, marketQuotes_;
  double SchobelIntegrand(double k, double maturity, double forward,
                          double strike, OptionType type) const;
  CalibrationTarget target_;
  void objFcn(int m, int n, double *x, double *fvec, int *iflag);
};

} // namespace velesquant
#endif