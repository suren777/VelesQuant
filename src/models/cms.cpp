#include <velesquant/models/cms.h>

#include <algorithm>
#include <string>
#include <vector>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>

using namespace std;
using namespace velesquant::xlw;

namespace velesquant {

cms::cms(double expirySR, double tenorSR, double forwardSR, double annuitySR,
         double /*payCMS*/, double discountCMS, double beta, double alpha,
         double nu, double rho)
    : discountCMS_(discountCMS) {
  Sabr *swaptionSABR = new Sabr(expirySR, forwardSR, beta, alpha, nu, rho);
  double atmVol = swaptionSABR->impliedVol(forwardSR);
  delete swaptionSABR;
  double forwardCMS =
      convAdj(expirySR, tenorSR, forwardSR, annuitySR, discountCMS, atmVol);
  cmsSABR_ = new Sabr(expirySR, forwardCMS, beta, alpha, nu, rho);
};

cms::cms(double expirySR, double tenorSR, double forwardSR, double annuitySR,
         double /*payCMS*/, double discountCMS, double beta,
         CellMatrix &strikes, CellMatrix &marketQuotes,
         CalibrationTarget quoteType)
    : discountCMS_(discountCMS) {
  Sabr *swaptionSABR = new Sabr(expirySR, forwardSR, beta);
  int m = strikes.RowsInStructure();
  vector<double> theStrikes(m), theQuotes(m);
  for (int i = 0; i < m; ++i) {
    theStrikes[i] = strikes(i, 0).NumericValue();
    theQuotes[i] = marketQuotes(i, 0).NumericValue();
  }
  swaptionSABR->calibrator(theStrikes, theQuotes, quoteType);
  double alpha = swaptionSABR->getParameterAlpha();
  double nu = swaptionSABR->getParameterNu();
  double rho = swaptionSABR->getParameterRho();
  double atmVol = swaptionSABR->impliedVol(forwardSR);
  delete swaptionSABR;
  double forwardCMS =
      convAdj(expirySR, tenorSR, forwardSR, annuitySR, discountCMS, atmVol);
  cmsSABR_ = new Sabr(expirySR, forwardCMS, beta, alpha, nu, rho);
};

// Added function by James Carey on 19-11-2015 to allow for alternative
// convexity adjustment calculation
cms::cms(double expirySR, double tenorSR, double forwardSR, double freqSR,
         double freqCMS, double /*payCMS*/, double discountCMS, double beta,
         CellMatrix &strikes, CellMatrix &marketQuotes,
         CalibrationTarget quoteType)
    : discountCMS_(discountCMS) {
  Sabr *swaptionSABR = new Sabr(expirySR, forwardSR, beta);
  int m = strikes.RowsInStructure();
  vector<double> theStrikes(m), theQuotes(m);
  for (int i = 0; i < m; ++i) {
    theStrikes[i] = strikes(i, 0).NumericValue();
    theQuotes[i] = marketQuotes(i, 0).NumericValue();
  }
  swaptionSABR->calibrator(theStrikes, theQuotes, quoteType);
  double alpha = swaptionSABR->getParameterAlpha();
  double nu = swaptionSABR->getParameterNu();
  double rho = swaptionSABR->getParameterRho();
  double atmVol = swaptionSABR->impliedVol(forwardSR);
  delete swaptionSABR;

  double forwardCMS =
      convAdjAlt(expirySR, tenorSR, forwardSR, freqSR, freqCMS, atmVol);
  cmsSABR_ = new Sabr(expirySR, forwardCMS, beta, alpha, nu, rho);
};

cms::cms(double expirySR, double tenorSR, double forwardSR, double annuitySR,
         double /*payCMS*/, double discountCMS, double beta,
         CellMatrix &strikes, CellMatrix &marketQuotes,
         CalibrationTarget quoteType, CellMatrix &initialParams)
    : discountCMS_(discountCMS) {
  Sabr *swaptionSABR = new Sabr(expirySR, forwardSR, beta);
  int m = strikes.RowsInStructure();
  vector<double> theStrikes(m), theQuotes(m);
  for (int i = 0; i < m; ++i) {
    theStrikes[i] = strikes(i, 0).NumericValue();
    theQuotes[i] = marketQuotes(i, 0).NumericValue();
  }
  vector<double> theParams(3);
  theParams[0] = initialParams(0, 0).NumericValue();
  theParams[1] = initialParams(1, 0).NumericValue();
  theParams[2] = initialParams(2, 0).NumericValue();
  // swaptionSABR->calibratorWithInitial(theStrikes,theQuotes,quoteType,theParams);
  swaptionSABR->calibratorWithInitial(theStrikes, theQuotes, quoteType);
  double alpha = swaptionSABR->getParameterAlpha();
  double nu = swaptionSABR->getParameterNu();
  double rho = swaptionSABR->getParameterRho();
  double atmVol = swaptionSABR->impliedVol(forwardSR);
  delete swaptionSABR;
  double forwardCMS =
      convAdj(expirySR, tenorSR, forwardSR, annuitySR, discountCMS, atmVol);
  cmsSABR_ = new Sabr(expirySR, forwardCMS, beta, alpha, nu, rho);
};

double cms::fairValue(double strike, OptionType callORput) const {
  return discountCMS_ * cmsSABR_->premiumBlackScholes(strike, callORput);
};

double cms::simulation(double corrRN) { return cmsSABR_->simulation(corrRN); };

double cms::getForward() const { return cmsSABR_->getForward(); };

double cms::getDiscountCMS() const { return discountCMS_; };

double cms::getMaturity() const { return cmsSABR_->getMaturity(); };

double cms::getATMvol() const {
  double forwardCMS = cmsSABR_->getForward();
  return cmsSABR_->impliedVol(forwardCMS);
};

double cms::getImpliedVol(double strike) const {
  return cmsSABR_->impliedVol(strike);
};

double cms::getParameterAlpha() const { return cmsSABR_->getParameterAlpha(); };

double cms::getParameterNu() const { return cmsSABR_->getParameterNu(); };

double cms::getParameterRho() const { return cmsSABR_->getParameterRho(); };

double cms::convAdj(double expirySR, double tenorSR, double forwardSR,
                    double annuitySR, double discountCMS, double atmVol) const {
  double durationReciprocal = 1.0 / tenorSR;
  double convexityCoeff =
      (discountCMS / annuitySR - durationReciprocal) / forwardSR;
  if (atmVol * atmVol * expirySR < 1.0 * atan(1.0))
    return forwardSR *
           (durationReciprocal +
            convexityCoeff * forwardSR * exp(atmVol * atmVol * expirySR)) /
           (durationReciprocal + convexityCoeff * forwardSR);
  else
    return forwardSR *
           (durationReciprocal +
            convexityCoeff * forwardSR * (1.0 + atmVol * atmVol * expirySR)) /
           (durationReciprocal + convexityCoeff * forwardSR);
};

// Added function by James Carey on 19-11-2015 to allow for alternative
// convexity adjustment calculation
double cms::convAdjAlt(double expirySR, double tenorSR, double forwardSR,
                       double freqSR, double freqCMS, double atmVol) const {
  double numPeriods = tenorSR * freqSR;
  double forwardOverFreq = forwardSR / freqSR;
  double adjustmentFactor =
      1 - forwardOverFreq *
              (freqSR / freqCMS +
               numPeriods / (pow(1.0 + forwardOverFreq, numPeriods) - 1.0)) /
              (1.0 + forwardOverFreq);
  double convexityAdj =
      forwardSR * adjustmentFactor * atmVol * atmVol * expirySR;
  return forwardSR + convexityAdj;
};
} // namespace velesquant