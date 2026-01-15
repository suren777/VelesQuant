//      cms_spread.h

#ifndef CMSSPREAD_H
#define CMSSPREAD_H

#include <string>
#include <vector>
#include <velesquant/models/cms.h>
#include <velesquant/types.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class cms_spread {
public:
  cms_spread(double expirySR1, double tenorSR1, double forwardSR1,
             double annuitySR1, double payCMS1, double discountCMS1,
             double beta1, CellMatrix &strikes1, CellMatrix &marketQuotes1,
             CalibrationTarget quoteType1, double expirySR2, double tenorSR2,
             double forwardSR2, double annuitySR2, double payCMS2,
             double discountCMS2, double beta2, CellMatrix &strikes2,
             CellMatrix &marketQuotes2, CalibrationTarget quoteType2,
             double corr);

  // Added function by James Carey on 19-11-2015 to allow for alternative
  // convexity adjustment calculation
  cms_spread(double expirySR1, double tenorSR1, double forwardSR1,
             double freqSR1, double freqCMS1, double payCMS1,
             double discountCMS1, double beta1, CellMatrix &strikes1,
             CellMatrix &marketQuotes1, CalibrationTarget quoteType1,
             double expirySR2, double tenorSR2, double forwardSR2,
             double freqSR2, double freqCMS2, double payCMS2,
             double discountCMS2, double beta2, CellMatrix &strikes2,
             CellMatrix &marketQuotes2, CalibrationTarget quoteType2,
             double corr);

  ~cms_spread();

  double margrabe() const;
  double spreadOption(double K, double a, double b) const;

  std::vector<double> simulationCMSs();
  std::vector<double> simulationCMSs(double corr);
  std::vector<double> simulationCMSs(double cr1, double cr2);

private:
  double corr_;
  cms *sabrCMS1_;
  cms *sabrCMS2_;

  double lstar(double epsilon, double K, double alpha, double sigma1,
               double beta, double sigma2, double rho) const;
  double intFun(double epsilon, int sign, double K, double alpha, double sigma1,
                double beta, double sigma2, double rho) const;
};

} // namespace velesquant
#endif