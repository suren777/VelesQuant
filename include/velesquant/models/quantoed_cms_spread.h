//      quantoedCMSspread.h

#ifndef QUANTOEDCMSSPREAD_H
#define QUANTOEDCMSSPREAD_H

#include <string>
#include <vector>
#include <velesquant/models/quantoed_cms.h>
#include <velesquant/types.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class quantoedCMSspread {
public:
  quantoedCMSspread(double expirySR1, double tenorSR1, double forwardSR1,
                    double annuitySR1, double payCMS1, double discountCMS1,
                    double corFX1, double atmVolFX1, double beta1,
                    CellMatrix &strikes1, CellMatrix &marketQuotes1,
                    CalibrationTarget quoteType1, double expirySR2,
                    double tenorSR2, double forwardSR2, double annuitySR2,
                    double payCMS2, double discountCMS2, double corFX2,
                    double atmVolFX2, double beta2, CellMatrix &strikes2,
                    CellMatrix &marketQuotes2, CalibrationTarget quoteType2,
                    double corr);

  ~quantoedCMSspread();

  std::vector<double> simulationQuantoedCMSs(double cr1, double cr2);

private:
  quantoedCMS *sabrQuantoedCMS1_;
  quantoedCMS *sabrQuantoedCMS2_;
};

} // namespace velesquant
#endif