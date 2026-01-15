//		Sabr.h

#ifndef SABR_H
#define SABR_H

#include <stdexcept>
#include <string>
#include <vector>
#include <velesquant/models/utility.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor

/**
 * @class Sabr
 * @brief SABR Model Implementation.
 *
 * Provides methods for calculating Implied Volatility (Black or
 * Normal/Bachelier), Pricing Options, Local Volatility, and Calibration.
 */
class Sabr {
public:
  Sabr() {};
  Sabr(double maturity, double forward, double beta = 0.85, double alpha = 0.5,
       double nu = 0.25, double rho = -0.75);
  Sabr(double maturity, double forward, double beta, double alpha, double nu,
       double rho, double shift);
  ~Sabr() {};

  /**
   * @brief Calculates Black Scholes Implied Volatility.
   * @param strike Strike price.
   */
  [[nodiscard]] double impliedVol(double strike) const;
  [[nodiscard]] double normalVol(double K) const;
  [[nodiscard]] double
  premiumBachelier(double strike,
                   OptionType callORput = OptionType::Call) const;
  [[nodiscard]] double
  premiumBlackScholes(double strike,
                      OptionType callORput = OptionType::Call) const;
  [[nodiscard]] double localVol(double spot) const;
  [[nodiscard]] double
  localVolCall(double spot) const; // based on the call premium
  [[nodiscard]] double localVolzabr(double spot) const;
  [[nodiscard]] double
  localVolIV(double spot) const; // based on the implied vol

  void qutlTable();
  double simulation(double corrRN);
  std::vector<double> simulations(std::vector<double> correlatedRNs);
  [[nodiscard]] std::vector<double> getTableQUTL1() const { return spots_; };
  [[nodiscard]] std::vector<double> getTableQUTL2() const { return qutls_; };
  [[nodiscard]] std::vector<double> getTableQUTL3() const { return cdfs_; };
  [[nodiscard]] double getVol(double strike) const;
  [[nodiscard]] double
  getPremium(double strike, OptionType callORput = OptionType::Call) const;
  void calibratorNormalVol(std::vector<double> strikes,
                           std::vector<double> marketQuotes);
  /**
   * @brief Calibrates SABR parameters to market quotes.
   *
   * @param strikes Vector of strikes.
   * @param marketQuotes Vector of market quotes (Premium or Implied Vol).
   * @param quoteType "premium" or "impliedVol".
   */
  void calibrator(std::vector<double> strikes, std::vector<double> marketQuotes,
                  CalibrationTarget quoteType = CalibrationTarget::Price);
  // void calibratorWithInitial(std::vector<double> strikes, std::vector<double>
  // marketQuotes, 	std::string quoteType, std::vector<double>
  // initialParams);
  void calibratorWithInitial(std::vector<double> strikes,
                             std::vector<double> marketQuotes,
                             CalibrationTarget quoteType);
  void calibratorWithInitialATM(std::vector<double> strikes,
                                std::vector<double> marketQuotes,
                                CalibrationTarget quoteType);

  void setParameterAlpha(double alpha) { alpha_ = alpha; };
  void setParameterNu(double nu) { nu_ = nu; };
  void setParameterRho(double rho) {
    if (rho * rho <= 1)
      rho_ = rho;
    else
      rho_ = sin(rho);
  };
  [[nodiscard]] double getParameterAlpha() const { return alpha_; };
  [[nodiscard]] double getParameterNu() const { return nu_; };
  [[nodiscard]] double getParameterRho() const { return rho_; };
  void setMaturity(double maturity) { maturity_ = maturity; };
  void setForward(double forward) { forward_ = forward; };
  void setBeta(double beta) { beta_ = beta; };
  void setATMvol(double atmVol) { atmVol_ = atmVol; };

  [[nodiscard]] double getMaturity() const { return maturity_; };
  [[nodiscard]] double getForward() const { return forward_; };
  [[nodiscard]] double getBeta() const { return beta_; };
  [[nodiscard]] double getShift() const { return shift_; };

private:
  bool calibrated_;
  double maturity_, forward_, beta_, shift_;
  double alpha_, nu_, rho_;
  double atmVol_;
  std::string type_;
  std::vector<double> strikes_;
  std::vector<double> marketQuotes_;

  double ATMvolPoly(double alpha);
  double ATMvolRoots(double lBound, double uBound, double tol);

  bool notQUTL_;
  std::vector<double> spots_;
  std::vector<double> qutls_;
  std::vector<double> cdfs_;
  void objFcn(int m, int n, double *x, double *fvec, int *iflag,
              CalibrationTarget quoteType = CalibrationTarget::Price);
  void objFcnATM(int m, int n, double *x, double *fvec, int *iflag,
                 CalibrationTarget quoteType);

  void parameters() {
    if (!((alpha_ >= 0.0) && (nu_ >= 0.0) && (rho_ * rho_ <= 1.0)))
      throw std::invalid_argument("SABR parameters are not valid");
    calibrated_ = true;
  };

  int amin(const std::vector<double> &q, int position);
  int amax(const std::vector<double> &q, int position);
  void asort(std::vector<double> &a, std::vector<double> &b,
             std::vector<double> &c);
  // Mesh initialisation
  std::vector<double> gsMesh(int);
  std::vector<double> grMesh(int);
  // Normal volatilities:
  double Nvolb0(double K) const; //	beta == 1
  double Nvolb(double K) const;  //	0 < beta <1
  double Nvolb1(double K) const; //	beta == 1
  // Obj function
  void objFcnNormalVolCalibration(int m, int n, double *x, double *fvec,
                                  int *iflag);
};

} // namespace velesquant
#endif