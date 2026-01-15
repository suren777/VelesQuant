//		cms.h

#ifndef CMS_H
#define CMS_H

#include <string>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class cms {
public:
  cms(double expirySR, double tenorSR, double forwardSR, double annuitySR,
      double payCMS, double discountCMS, double beta = 0.85, double alpha = 0.5,
      double nu = 0.25, double rho = -0.75);

  // Added function by James Carey on 19-11-2015 to allow for alternative
  // convexity adjustment calculation
  cms(double expirySR, double tenorSR, double forwardSR, double freqSR,
      double freqCMS, double payCMS, double discountCMS, double beta,
      CellMatrix &strikes, CellMatrix &marketQuotes,
      CalibrationTarget quoteType = CalibrationTarget::Price);

  cms(double expirySR, double tenorSR, double forwardSR, double annuitySR,
      double payCMS, double discountCMS, double beta, CellMatrix &strikes,
      CellMatrix &marketQuotes,
      CalibrationTarget quoteType = CalibrationTarget::Price);

  cms(double expirySR, double tenorSR, double forwardSR, double annuitySR,
      double payCMS, double discountCMS, double beta, CellMatrix &strikes,
      CellMatrix &marketQuotes, CalibrationTarget quoteType,
      CellMatrix &initialParams);

  ~cms() { delete cmsSABR_; };

  double fairValue(double strike,
                   OptionType callORput = OptionType::Call) const;

  double simulation(double corrRN);

  double getForward() const;
  double getDiscountCMS() const;

  double getMaturity() const;
  double getATMvol() const;
  double getImpliedVol(double strike) const;

  double getParameterAlpha() const;
  double getParameterNu() const;
  double getParameterRho() const;

private:
  double discountCMS_;
  Sabr *cmsSABR_;

  double convAdj(double expirySR, double tenorSR, double forwardSR,
                 double annuitySR, double discountCMS, double atmVol) const;

  // Added function by James Carey on 19-11-2015 to allow for alternative
  // convexity adjustment calculation
  double convAdjAlt(double expirySR, double tenorSR, double forwardSR,
                    double freqSR, double freqCMS, double atmVol) const;
};

} // namespace velesquant
#endif