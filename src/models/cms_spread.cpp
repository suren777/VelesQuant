#include <velesquant/models/cms_spread.h>

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <ql/quantlib.hpp>
#include <string>
#include <vector>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/cms.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>

using namespace std;
using namespace velesquant::xlw;
using namespace QuantLib;
using namespace boost::placeholders;

namespace velesquant {

cms_spread::cms_spread(double expirySR1, double tenorSR1, double forwardSR1,
                       double annuitySR1, double payCMS1, double discountCMS1,
                       double beta1, CellMatrix &strikes1,
                       CellMatrix &marketQuotes1, CalibrationTarget quoteType1,
                       double expirySR2, double tenorSR2, double forwardSR2,
                       double annuitySR2, double payCMS2, double discountCMS2,
                       double beta2, CellMatrix &strikes2,
                       CellMatrix &marketQuotes2, CalibrationTarget quoteType2,
                       double corr)
    : corr_(corr) {
  sabrCMS1_ = new cms(expirySR1, tenorSR1, forwardSR1, annuitySR1, payCMS1,
                      discountCMS1, beta1, strikes1, marketQuotes1, quoteType1);
  sabrCMS2_ = new cms(expirySR2, tenorSR2, forwardSR2, annuitySR2, payCMS2,
                      discountCMS2, beta2, strikes2, marketQuotes2, quoteType2);
};

// Added function by James Carey on 19-11-2015 to allow for alternative
// convexity adjustment calculation
cms_spread::cms_spread(double expirySR1, double tenorSR1, double forwardSR1,
                       double freqSR1, double freqCMS1, double payCMS1,
                       double discountCMS1, double beta1, CellMatrix &strikes1,
                       CellMatrix &marketQuotes1, CalibrationTarget quoteType1,
                       double expirySR2, double tenorSR2, double forwardSR2,
                       double freqSR2, double freqCMS2, double payCMS2,
                       double discountCMS2, double beta2, CellMatrix &strikes2,
                       CellMatrix &marketQuotes2, CalibrationTarget quoteType2,
                       double corr)
    : corr_(corr) {
  sabrCMS1_ =
      new cms(expirySR1, tenorSR1, forwardSR1, freqSR1, freqCMS1, payCMS1,
              discountCMS1, beta1, strikes1, marketQuotes1, quoteType1);
  sabrCMS2_ =
      new cms(expirySR2, tenorSR2, forwardSR2, freqSR2, freqCMS2, payCMS2,
              discountCMS2, beta2, strikes2, marketQuotes2, quoteType2);
};

cms_spread::~cms_spread() {
  delete sabrCMS1_;
  delete sabrCMS2_;
};

double cms_spread::margrabe() const {
  double discountCMS = sabrCMS1_->getDiscountCMS();
  double Maturity = sabrCMS1_->getMaturity();
  double atmVol1 = sabrCMS1_->getATMvol();
  double atmVol2 = sabrCMS2_->getATMvol();
  double cmsFwd1 = sabrCMS1_->getForward();
  double cmsFwd2 = sabrCMS2_->getForward();
  double sigmaT = sqrt(Maturity) * sqrt(atmVol1 * atmVol1 + atmVol2 * atmVol2 -
                                        2 * corr_ * atmVol1 * atmVol2);
  double d1 = log(cmsFwd1 / cmsFwd2) / sigmaT + 0.5 * sigmaT;
  double d2 = d1 - sigmaT;
  return discountCMS * (cmsFwd1 * cdf_normal(d1) - cmsFwd2 * cdf_normal(d2));
};

double cms_spread::spreadOption(double K, double a, double b) const // p151
{
  double discountCMS = sabrCMS1_->getDiscountCMS();
  double sU = sqrt(sabrCMS1_->getMaturity());
  double sT = sqrt(sabrCMS2_->getMaturity());
  double sigma1 = sabrCMS1_->getATMvol() * sU;
  double sigma2 = sabrCMS2_->getATMvol() * sT;
  double rho = corr_ * sU / sT;
  double alpha = a * sabrCMS1_->getForward() * exp(-0.5 * sigma1 * sigma1);
  double beta = b * sabrCMS2_->getForward() * exp(-0.5 * sigma2 * sigma2);
  int pn1 = 1;
  if (b > 0)
    pn1 = -1;
  boost::function<double(double)> fcn;
  fcn = boost::bind(&cms_spread::intFun, this, _1, pn1, K, alpha, sigma1, beta,
                    sigma2, rho);
  GaussHermiteIntegration gHerInt(16);
  double integr = gHerInt(fcn);
  /*
          double absAcc = 1.0E-5; int maxEval = 1000;
          double low = -1.0E+10, upp = 1.0E+10;
          SimpsonIntegral numInt(absAcc,maxEval);
          double integr = numInt(fcn,low,upp);
  */
  return discountCMS * integr;
};

double cms_spread::intFun(double epsilon, int pn1, double K, double alpha,
                          double sigma1, double beta, double sigma2,
                          double rho) const {
  double d1 = lstar(epsilon, K, alpha, sigma1, beta, sigma2, rho);
  double d2 = d1 - sigma2 * sqrt(1 - rho * rho);
  double int1 = (alpha * exp(sigma1 * epsilon) - K) * cdf_normal(pn1 * d1);
  double int2 = exp(rho * sigma2 * epsilon) * cdf_normal(pn1 * d2);
  return pdf_normal(epsilon) *
         (int1 + beta * exp(0.5 * sigma2 * sigma2 * (1 - rho * rho)) * int2);
};

double cms_spread::lstar(double epsilon, double K, double alpha, double sigma1,
                         double beta, double sigma2, double rho) const {
  double pn = (K - alpha * exp(sigma1 * epsilon)) /
              (beta * exp(rho * sigma2 * epsilon));
  double ls = -1.0E+20;
  if (pn > 0.0)
    ls = log(pn) / (sigma2 * sqrt(1 - rho * rho));
  return ls;
};

vector<double> cms_spread::simulationCMSs() {
  double r1 = random_normal();
  double r2 = corr_ * r1 + sqrt(1 - corr_ * corr_) * random_normal();
  vector<double> cmss(2);
  cmss[0] = sabrCMS1_->simulation(r1);
  cmss[1] = sabrCMS2_->simulation(r2);
  return cmss;
};
vector<double> cms_spread::simulationCMSs(double corr) {
  double r1 = random_normal();
  double r2 = corr * r1 + sqrt(1 - corr * corr) * random_normal();
  vector<double> cmss(2);
  cmss[0] = sabrCMS1_->simulation(r1);
  cmss[1] = sabrCMS2_->simulation(r2);
  return cmss;
};
vector<double> cms_spread::simulationCMSs(double cr1, double cr2) {
  vector<double> cmss(2);
  cmss[0] = sabrCMS1_->simulation(cr1);
  cmss[1] = sabrCMS2_->simulation(cr2);
  return cmss;
};

} // namespace velesquant