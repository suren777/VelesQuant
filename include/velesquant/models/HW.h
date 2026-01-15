//		HullWhite.h

#ifndef HW_H
#define HW_H

#include <vector>
#include <velesquant/models/utility.h>
#include <velesquant/numerics/interpolation.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor

/**
 * @class HullWhite
 * @brief Hull-White 1-Factor Model Implementation.
 *
 * Provides methods for pricing bonds, swaptions, and simulation.
 * Supports calibration to market swaptions (Price or Volatility).
 */
class HullWhite {
public:
  HullWhite(double kappa, std::vector<double> timeSigmas,
            std::vector<double> sigmas, std::vector<double> timeDFs,
            std::vector<double> DFs)
      : kappa_(kappa), timeDFs_(timeDFs), DFs_(DFs), timeSigmas_(timeSigmas),
        sigmas_(sigmas) {
    sigmas0_ = sigmas_;
    kappa0_ = kappa_;
  };

  ~HullWhite() {};

  /**
   * @brief Calculates the price of a Bond Option.
   *
   * @param Expiry Option expiry time.
   * @param Maturity Bond maturity time.
   * @param Strike Strike price.
   * @param type Option type (Call or Put).
   * @return Option price.
   */
  [[nodiscard]] double optionBond(double Expiry, double Maturity, double Strike,
                                  OptionType type);
  [[nodiscard]] double swaption(double Expiry, double Tenor, double Strike,
                                double PayFrequency = 0.5);
  [[nodiscard]] double BlackStrike(double T0, double TN, double impVol);
  [[nodiscard]] double BlackStrikePlain(double T0, double TN);
  [[nodiscard]] std::vector<double> simulation(std::vector<double> times) const;

  [[nodiscard]] double ZC(double Expiry);
  /**
   * @brief Calibrates model parameters (sigma) to market swaptions.
   *
   * @param swapQuotes Vector of market swaption quotes.
   * @param target Calibration target (Price or Volatility).
   */
  void calibrator(const std::vector<defSwap> &swapQuotes,
                  CalibrationTarget target);
  void calibratorBstrp(const std::vector<defSwap> &swapQuotes,
                       CalibrationTarget target);
  [[nodiscard]] double getKappa() { return kappa_; };
  [[nodiscard]] std::vector<double> getTimeSigmas() { return timeSigmas_; };
  [[nodiscard]] std::vector<double> getSigmas() { return sigmas_; };
  [[nodiscard]] double getSwapRate(double Expiry, double Tenor,
                                   double PayFrequency = 0.5);
  [[nodiscard]] double swaptionIVblackPub(double Expiry, double Tenor,
                                          double swap_price);
  [[nodiscard]] double get_swaptionATM(double Expiry, double Tenor,
                                       double VolATM) {
    return swaptionATM(Expiry, Tenor, VolATM);
  };
  [[nodiscard]] std::vector<double> getDFs() { return DFs_; };
  [[nodiscard]] std::vector<double> getDFsTimes() { return timeDFs_; };

private:
  double kappa_, kappa0_;
  std::vector<double> timeDFs_, DFs_; // time_[0] first point
  mutable std::vector<double> timeSigmas_, sigmas_, sigmas0_;
  int iter_;
  mutable std::vector<double> marketSwaption_;
  std::vector<defSwap> quoteSwap_;
  std::vector<defSwap> qSB_;

  double totalVariance(double T0);
  double equalityCP(double vX, double vO, double Expiry, double Tenor,
                    double Strike, double PayFrequency);
  double criticalPoint(double Expiry, double Tenor, double Strike,
                       double PayFrequency);

  double formulaBlack(double dfT0, double dfTN, double Strike, double impVol,
                      double T0, OptionType type);

  double getDF(double T) {
    return interpolation("CubicNaturalSpline", timeDFs_, DFs_, T);
  };
  double swaptionATM(double Expiry, double Tenor, double VolATM) {
    return (getDF(Expiry) - getDF(Expiry + Tenor)) *
           (2.0 * cdf_normal(0.5 * VolATM * sqrt(Expiry)) - 1.0);
  };

  double whichSigma(double t) const;
  double whichFwdRate(double t) const;

  double pen_fun(double *x, double *lb, double *ub, int n);
  void objFcn(int m, int n, double *x, double *fvec, int *iflag);
  void objFcnIV(int m, int n, double *x, double *fvec, int *iflag, double *lb,
                double *ub);
  //	void objFcnPrice(int m, int n, double* x, double* fvec, int* iflag);
  void objFcnPrice(int m, int n, double *x, double *fvec, int *iflag,
                   double *lb, double *ub);
  double objFcnBswap(double x, int m);
  double objFcnBIV(double x, int m);
  // not used
  double swap(double Expiry, double Tenor, double Strike,
              double PayFrequency = 0.5);
  double swaptionIVblack(double Expiry, double Tenor, double swaption_price);
  double Bratio(double a, double Mi, double Tj, double Tk);
  void calibrateKappa();
  mutable std::vector<double> iv_ratio, bond_ratio, index;
  double objKappa(double x, const std::vector<double> &ivr,
                  const std::vector<double> &br, const std::vector<int> &ind);
};

} // namespace velesquant
#endif