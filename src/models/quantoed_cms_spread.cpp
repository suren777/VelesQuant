#include <velesquant/models/quantoed_cms_spread.h>

#include <algorithm>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <ql/quantlib.hpp>
#include <string>
#include <vector>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/quantoed_cms.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>

using namespace std;
using namespace velesquant::xlw;
using namespace QuantLib;

namespace velesquant {

quantoedCMSspread::quantoedCMSspread(
    double expirySR1, double tenorSR1, double forwardSR1, double annuitySR1,
    double payCMS1, double discountCMS1, double corFX1, double atmVolFX1,
    double beta1, CellMatrix &strikes1, CellMatrix &marketQuotes1,
    CalibrationTarget quoteType1, double expirySR2, double tenorSR2,
    double forwardSR2, double annuitySR2, double payCMS2, double discountCMS2,
    double corFX2, double atmVolFX2, double beta2, CellMatrix &strikes2,
    CellMatrix &marketQuotes2, CalibrationTarget quoteType2, double corr) {
  sabrQuantoedCMS1_ = new quantoedCMS(
      expirySR1, tenorSR1, forwardSR1, annuitySR1, payCMS1, discountCMS1,
      corFX1, atmVolFX1, beta1, strikes1, marketQuotes1, quoteType1);
  sabrQuantoedCMS2_ = new quantoedCMS(
      expirySR2, tenorSR2, forwardSR2, annuitySR2, payCMS2, discountCMS2,
      corFX2, atmVolFX2, beta2, strikes2, marketQuotes2, quoteType2);
};

quantoedCMSspread::~quantoedCMSspread() {
  delete sabrQuantoedCMS1_;
  delete sabrQuantoedCMS2_;
};

vector<double> quantoedCMSspread::simulationQuantoedCMSs(double cr1,
                                                         double cr2) {
  vector<double> quantoedCMSs(2);
  quantoedCMSs[0] = sabrQuantoedCMS1_->simulation(cr1);
  quantoedCMSs[1] = sabrQuantoedCMS2_->simulation(cr2);
  return quantoedCMSs;
};
} // namespace velesquant