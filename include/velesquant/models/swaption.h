//		swaption.h

#ifndef SWAPTION_H
#define SWAPTION_H

#include <string>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class swaption {
public:
  swaption(double expiry, double tenor, double forward, double annuity,
           double beta = 0.85, double alpha = 0.5, double nu = 0.25,
           double rho = -0.75);

  swaption(double expiry, double tenor, double forward, double annuity,
           double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
           CalibrationTarget quoteType = CalibrationTarget::Price);

  swaption(double expiry, double tenor, double forward, double annuity,
           double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
           CalibrationTarget quoteType, CellMatrix &initialParams);

  ~swaption() { delete swaptionSABR_; };

  double swaptionFairValue(double strike,
                           OptionType callORput = OptionType::Call) const;
  double swapFairValue(double strike) const;

  double simulation(double corrRN);

  double getImpliedVol(double strike) const;

  double getParameterAlpha() const;
  double getParameterNu() const;
  double getParameterRho() const;

private:
  double annuity_;
  Sabr *swaptionSABR_;
};

} // namespace velesquant
#endif