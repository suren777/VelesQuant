//		quantoedCMS.h

#ifndef QUANTOEDCMS_H
#define QUANTOEDCMS_H

#include <string>
#include <velesquant/volatility/sabr.h>
#include <velesquant/types.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class quantoedCMS {
public:
  quantoedCMS(double expirySR, double tenorSR, double forwardSR,
              double annuitySR, double payCMS, double discountCMS, double corFX,
              double atmVolFX, double beta, CellMatrix &strikes,
              CellMatrix &marketQuotes,
              CalibrationTarget quoteType = CalibrationTarget::Price);

  ~quantoedCMS() { delete quantoedCMSSABR_; };

  double fairValue(double strike,
                   OptionType callORput = OptionType::Call) const;
  double simulation(double corrRN);
  double getForward() const;

private:
  Sabr *quantoedCMSSABR_;

  double quantoAdj(double expirySR, double forwardSR, double atmVol,
                   double corFX, double atmVolFX) const;
  double convAdj(double expirySR, double tenorSR, double forwardSR,
                 double quantoedForwardSR, double annuitySR, double discountCMS,
                 double atmVol) const;
};

} // namespace velesquant
#endif