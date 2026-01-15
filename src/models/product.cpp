
#include <boost/function.hpp>
#include <ql/quantlib.hpp>
#include <string>
#include <vector>
#include <velesquant/models/cms.h>
#include <velesquant/models/cms_spread.h>
#include <velesquant/models/log_basket.h>
#include <velesquant/models/object_cache.h>
#include <velesquant/models/quantoed_cms.h>
#include <velesquant/models/quantoed_cms_spread.h>
#include <velesquant/models/swaption.h>
#include <velesquant/models/utility.h>
// Removed: ArgList no longer needed
#include <algorithm>
#include <velesquant/types.h>

using namespace std;
using namespace velesquant::xlw;
using namespace QuantLib;

#pragma warning(disable : 4996)

namespace velesquant {

CellMatrix simulationQuantoedCMSs(
    double expirySR1, double tenorSR1, double forwardSR1, double annuitySR1,
    double payCMS1, double discountCMS1, double corFX1, double atmVolFX1,
    double beta1, CellMatrix &strikes1, CellMatrix &marketQuotes1,
    string quoteType1, double expirySR2, double tenorSR2, double forwardSR2,
    double annuitySR2, double payCMS2, double discountCMS2, double corFX2,
    double atmVolFX2, double beta2, CellMatrix &strikes2,
    CellMatrix &marketQuotes2, string quoteType2, double corr, int NP) {
  quantoedCMSspread *model = new quantoedCMSspread(
      expirySR1, tenorSR1, forwardSR1, annuitySR1, payCMS1, discountCMS1,
      corFX1, atmVolFX1, beta1, strikes1, marketQuotes1,
      stringToCalibrationTarget(quoteType1), expirySR2, tenorSR2, forwardSR2,
      annuitySR2, payCMS2, discountCMS2, corFX2, atmVolFX2, beta2, strikes2,
      marketQuotes2, stringToCalibrationTarget(quoteType2), corr);
  CellMatrix theQuantoedCMSs(2, NP);
  double cor2 = sqrt(1 - corr * corr);
  for (int j = 0; j < NP; j++) {
    double cr1 = random_normal();
    double cr2 = corr * cr1 + cor2 * random_normal();
    vector<double> oneCMS = model->simulationQuantoedCMSs(cr1, cr2);
    theQuantoedCMSs(0, j) = oneCMS[0];
    theQuantoedCMSs(1, j) = oneCMS[1];
  }
  delete model;
  return theQuantoedCMSs;
}

std::string
CreateQuantoedCMSObj(const std::string &theObjName, double expirySR,
                     double tenorSR, double forwardSR, double annuitySR,
                     double payCMS, double discountCMS, double corFX,
                     double atmVolFX, double beta, CellMatrix &strikesSR,
                     CellMatrix &marketQuotesSR, std::string quoteTypeSR) {
  return placeObject(theObjName,
                     new quantoedCMS(expirySR, tenorSR, forwardSR, annuitySR,
                                     payCMS, discountCMS, corFX, atmVolFX, beta,
                                     strikesSR, marketQuotesSR,
                                     stringToCalibrationTarget(quoteTypeSR)));
}
double quantoedCMSForward(const std::string &theName) {
  return getObject<quantoedCMS>(theName)->getForward();
}
double optionQuantoedCMS(const std::string &theName, double theStrike,
                         std::string callORput) {
  return getObject<quantoedCMS>(theName)->fairValue(
      theStrike, stringToOptionType(callORput));
}
CellMatrix simulationQuantoedCMS(const std::string &theName, int NP) {
  quantoedCMS *cms = getObject<quantoedCMS>(theName);
  CellMatrix quantoedCMSpath(1, NP);
  for (int j = 0; j < NP; j++)
    quantoedCMSpath(0, j) = cms->simulation(random_normal());
  return quantoedCMSpath;
}

double cmsSpreadOption(double expirySR1, double tenorSR1, double forwardSR1,
                       double annuitySR1, double payCMS1, double discountCMS1,
                       double beta1, CellMatrix &strikes1,
                       CellMatrix &marketQuotes1, string quoteType1,
                       double expirySR2, double tenorSR2, double forwardSR2,
                       double annuitySR2, double payCMS2, double discountCMS2,
                       double beta2, CellMatrix &strikes2,
                       CellMatrix &marketQuotes2, string quoteType2,
                       double corr, double strike, double coef1, double coef2) {
  cms_spread *model = new cms_spread(
      expirySR1, tenorSR1, forwardSR1, annuitySR1, payCMS1, discountCMS1, beta1,
      strikes1, marketQuotes1, stringToCalibrationTarget(quoteType1), expirySR2,
      tenorSR2, forwardSR2, annuitySR2, payCMS2, discountCMS2, beta2, strikes2,
      marketQuotes2, stringToCalibrationTarget(quoteType2), corr);
  double result = model->spreadOption(strike, coef1, coef2);
  delete model;
  return result;
}

double cmsMargrabe(double expirySR1, double tenorSR1, double forwardSR1,
                   double annuitySR1, double payCMS1, double discountCMS1,
                   double beta1, CellMatrix &strikes1,
                   CellMatrix &marketQuotes1, string quoteType1,
                   double expirySR2, double tenorSR2, double forwardSR2,
                   double annuitySR2, double payCMS2, double discountCMS2,
                   double beta2, CellMatrix &strikes2,
                   CellMatrix &marketQuotes2, string quoteType2, double corr) {
  cms_spread *model = new cms_spread(
      expirySR1, tenorSR1, forwardSR1, annuitySR1, payCMS1, discountCMS1, beta1,
      strikes1, marketQuotes1, stringToCalibrationTarget(quoteType1), expirySR2,
      tenorSR2, forwardSR2, annuitySR2, payCMS2, discountCMS2, beta2, strikes2,
      marketQuotes2, stringToCalibrationTarget(quoteType2), corr);
  double result = model->margrabe();
  delete model;
  return result;
}

CellMatrix simulationCMS(double expirySR1, double tenorSR1, double forwardSR1,
                         double annuitySR1, double payCMS1, double discountCMS1,
                         double beta1, CellMatrix &strikes1,
                         CellMatrix &marketQuotes1, string quoteType1,
                         double expirySR2, double tenorSR2, double forwardSR2,
                         double annuitySR2, double payCMS2, double discountCMS2,
                         double beta2, CellMatrix &strikes2,
                         CellMatrix &marketQuotes2, string quoteType2,
                         double corr, int NP) {
  cms_spread *model = new cms_spread(
      expirySR1, tenorSR1, forwardSR1, annuitySR1, payCMS1, discountCMS1, beta1,
      strikes1, marketQuotes1, stringToCalibrationTarget(quoteType1), expirySR2,
      tenorSR2, forwardSR2, annuitySR2, payCMS2, discountCMS2, beta2, strikes2,
      marketQuotes2, stringToCalibrationTarget(quoteType2), corr);

  CellMatrix theCMSs(2, NP);
  double cor2 = sqrt(1 - corr * corr);
  for (int j = 0; j < NP; j++) {
    double cr1 = random_normal();
    double cr2 = corr * cr1 + cor2 * random_normal();
    vector<double> oneCMS = model->simulationCMSs(cr1, cr2);
    theCMSs(0, j) = oneCMS[0];
    theCMSs(1, j) = oneCMS[1];
  }
  delete model;
  return theCMSs;
}

CellMatrix simulationCMSAltConvAdj(
    double expirySR1, double tenorSR1, double forwardSR1, double freqSR1,
    double freqCMS1, double payCMS1, double discountCMS1, double beta1,
    CellMatrix &strikes1, CellMatrix &marketQuotes1, string quoteType1,
    double expirySR2, double tenorSR2, double forwardSR2, double freqSR2,
    double freqCMS2, double payCMS2, double discountCMS2, double beta2,
    CellMatrix &strikes2, CellMatrix &marketQuotes2, string quoteType2,
    double corr, int NP) {
  cms_spread *model = new cms_spread(
      expirySR1, tenorSR1, forwardSR1, freqSR1, freqCMS1, payCMS1, discountCMS1,
      beta1, strikes1, marketQuotes1, stringToCalibrationTarget(quoteType1),
      expirySR2, tenorSR2, forwardSR2, freqSR2, freqCMS2, payCMS2, discountCMS2,
      beta2, strikes2, marketQuotes2, stringToCalibrationTarget(quoteType2),
      corr);

  CellMatrix theCMSs(2, NP);
  double cor2 = sqrt(1 - corr * corr);
  for (int j = 0; j < NP; j++) {
    double cr1 = random_normal();
    double cr2 = corr * cr1 + cor2 * random_normal();
    vector<double> oneCMS = model->simulationCMSs(cr1, cr2);
    theCMSs(0, j) = oneCMS[0];
    theCMSs(1, j) = oneCMS[1];
  }
  delete model;
  return theCMSs;
}

std::string CreateCMSObj(const std::string &theObjName, double expirySR,
                         double tenorSR, double forwardSR, double annuitySR,
                         double payCMS, double discountCMS, double beta,
                         CellMatrix &strikesSR, CellMatrix &marketQuotesSR,
                         std::string quoteTypeSR) {
  return placeObject(theObjName,
                     new cms(expirySR, tenorSR, forwardSR, annuitySR, payCMS,
                             discountCMS, beta, strikesSR, marketQuotesSR,
                             stringToCalibrationTarget(quoteTypeSR)));
}

// Added function by James Carey on 19-11-2015 to allow for alternative
// convexity adjustment calculation
std::string
CreateCMSObjAltConvAdj(const std::string &theObjName, double expirySR,
                       double tenorSR, double forwardSR, double freqSR,
                       double freqCMS, double payCMS, double discountCMS,
                       double beta, CellMatrix &strikesSR,
                       CellMatrix &marketQuotesSR, string quoteTypeSR) {
  return placeObject(
      theObjName, new cms(expirySR, tenorSR, forwardSR, freqSR, freqCMS, payCMS,
                          discountCMS, beta, strikesSR, marketQuotesSR,
                          stringToCalibrationTarget(quoteTypeSR)));
}

std::string
CreateCMSObjectWithInitial(const std::string &theObjName, double expirySR,
                           double tenorSR, double forwardSR, double annuitySR,
                           double payCMS, double discountCMS, double beta,
                           CellMatrix &strikesSR, CellMatrix &marketQuotesSR,
                           std::string quoteTypeSR, CellMatrix &initialParams) {
  return placeObject(theObjName,
                     new cms(expirySR, tenorSR, forwardSR, annuitySR, payCMS,
                             discountCMS, beta, strikesSR, marketQuotesSR,
                             stringToCalibrationTarget(quoteTypeSR),
                             initialParams));
}
double valueCMS(const std::string &theName, double theStrike,
                std::string callORput) {
  return getObject<cms>(theName)->fairValue(theStrike,
                                            stringToOptionType(callORput));
}
double forwardCMS(const std::string &theName) {
  return getObject<cms>(theName)->getForward();
}
CellMatrix cmsSABRparameters(const std::string &theName) {
  cms *c = getObject<cms>(theName);
  CellMatrix sabrPara(3, 1);
  sabrPara(0, 0) = c->getParameterAlpha();
  sabrPara(1, 0) = c->getParameterNu();
  sabrPara(2, 0) = c->getParameterRho();
  return sabrPara;
}
double cmsOptionImpliedVol(const std::string &theName, double theStrike) {
  return getObject<cms>(theName)->getImpliedVol(theStrike);
}

std::string CreateSwaptionObj(const std::string &theObjName, double expiry,
                              double tenor, double forward, double annuity,
                              double beta, CellMatrix &strikes,
                              CellMatrix &marketQuotes, std::string quoteType) {
  return placeObject(theObjName,
                     new swaption(expiry, tenor, forward, annuity, beta,
                                  strikes, marketQuotes,
                                  stringToCalibrationTarget(quoteType)));
}
std::string CreateSwaptionObjectWithInitial(
    const std::string &theObjName, double expiry, double tenor, double forward,
    double annuity, double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
    std::string quoteType, CellMatrix &initialParams) {
  return placeObject(
      theObjName,
      new swaption(expiry, tenor, forward, annuity, beta, strikes, marketQuotes,
                   stringToCalibrationTarget(quoteType), initialParams));
}
double valueSwaption(const std::string &theName, double theStrike,
                     std::string callORput) {
  return getObject<swaption>(theName)->swaptionFairValue(
      theStrike, stringToOptionType(callORput));
}
double valueSwap(const std::string &theName, double theStrike) {
  return getObject<swaption>(theName)->swapFairValue(theStrike);
}
CellMatrix swaptionSABRparameters(const std::string &theName) {
  swaption *s = getObject<swaption>(theName);
  CellMatrix sabrPara(3, 1);
  sabrPara(0, 0) = s->getParameterAlpha();
  sabrPara(1, 0) = s->getParameterNu();
  sabrPara(2, 0) = s->getParameterRho();
  return sabrPara;
}
double swaptionImpliedVol(const std::string &theName, double theStrike) {
  return getObject<swaption>(theName)->getImpliedVol(theStrike);
}

int qldate() {
  QuantLib::Calendar myCal = QuantLib::UnitedKingdom();
  QuantLib::Date NYeve(31, QuantLib::Dec, 2008);

  return myCal.isBusinessDay(NYeve);
}

std::string testQL() {
  Calendar myCAL = QuantLib::UnitedKingdom();
  Date date1(28, QuantLib::Dec, 2008);
  QuantLib::Date date2(04, QuantLib::Jan, 2009);

  return myCAL.name();
}

std::string CreateBasketObj(const std::string &theObjName, CellMatrix &Spot,
                            CellMatrix &Strike, CellMatrix &Maturities,
                            CellMatrix &Forwards, CellMatrix &IV,
                            CellMatrix &correlation) {

  int n = Spot.RowsInStructure();
  int m = Maturities.ColumnsInStructure();
  vector<double> spot(n), strike(n), maturities(m);
  vector<vector<double>> forwards(n, vector<double>(m)),
      iv(n, vector<double>(m)), corr(n, vector<double>(n));
  for (int i = 0; i < n; i++) {
    spot[i] = Spot(i, 0).NumericValue();
    strike[i] = Strike(i, 0).NumericValue();
    for (int j = 0; j < m; j++) {
      if (i == 0)
        maturities[j] = Maturities(0, j).NumericValue();
      forwards[i][j] = Forwards(i, j).NumericValue();
      iv[i][j] = IV(i, j).NumericValue();
    }
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      corr[i][j] = correlation(i, j).NumericValue();

  return placeObject(theObjName,
                     new lBasket(spot, strike, maturities, forwards, iv, corr));
}

CellMatrix basketSimulation(const std::string &theName, CellMatrix Schedule,
                            int seed) {
  lBasket *b = getObject<lBasket>(theName);
  int m = Schedule.RowsInStructure();
  std::vector<double> sched(m);
  for (int i = 0; i < m; i++) {
    if (!Schedule(i, 0).IsANumber())
      break;
    sched[i] = Schedule(i, 0).NumericValue();
  }
  int Nass = b->get_nassets();
  CellMatrix thePaths(Nass, m);
  b->update_seed(seed);
  vector<vector<double>> pathF = b->simulate_basket(sched);
  for (int j = 0; j < Nass; j++)
    for (int i = 0; i < m; i++)
      thePaths(j, i) = pathF[j][i];

  return thePaths;
}

CellMatrix WorstAutocalSimulation(const std::string &theName,
                                  CellMatrix Schedule, double Coupon,
                                  double Notional, double H, double B, int NP,
                                  int seed) {
  lBasket *b = getObject<lBasket>(theName);
  int m = Schedule.RowsInStructure();
  std::vector<double> sched(m);
  for (int i = 0; i < m; i++) {
    if (!Schedule(i, 0).IsANumber())
      break;
    sched[i] = Schedule(i, 0).NumericValue();
  }
  int Nass = b->get_nassets();
  CellMatrix payoffs(NP, m);
  b->update_seed(seed);
  for (int i = 0; i < NP; i++) {
    int call_flag = 1;
    vector<double> pathF = b->simulate_basketWR(sched);
    for (int j = 0; j < m; j++) {
      payoffs(i, j) = Notional * Coupon * (pathF[j] >= B) * double(call_flag);
      if (pathF[j] >= H) {
        payoffs(i, j) =
            payoffs(i, j).NumericValue() + Notional * double(call_flag);
        call_flag = 0;
      }
      if ((j == m - 1) && (call_flag == 1))
        payoffs(i, j) = payoffs(i, j).NumericValue() + Notional;
    }
  }
  return payoffs;
}

CellMatrix WorstAutocalWithRemovalSimulation(
    const std::string &theName, CellMatrix Schedule, CellMatrix Barriers,
    CellMatrix Coupons, CellMatrix InitialPrice, CellMatrix Ishares,
    double Notional, int Nsteps, int NP, int seed, int called) {
  lBasket *b = getObject<lBasket>(theName);
  int m = Schedule.ColumnsInStructure();
  std::vector<double> sched(m + 1);
  std::vector<int> ishares(Ishares.RowsInStructure());
  std::vector<double> barriers(Barriers.RowsInStructure());
  std::vector<double> coupons(Coupons.RowsInStructure());
  std::vector<double> initial(InitialPrice.RowsInStructure());
  std::vector<double> spots(Ishares.RowsInStructure());
  for (int i = 0; i < int(barriers.size()); i++)
    barriers[i] = Barriers(i, 0).NumericValue();
  for (int i = 0; i < int(coupons.size()); i++)
    coupons[i] = Coupons(i, 0).NumericValue();
  for (int i = 0; i < int(initial.size()); i++)
    initial[i] = InitialPrice(i, 0).NumericValue();
  for (int i = 0; i < int(ishares.size()); i++)
    ishares[i] = int(Ishares(i, 0).NumericValue());
  sched[0] = 0.0;
  for (int i = 0; i < m; i++) {
    if (!Schedule(0, i).IsANumber())
      break;
    sched[i + 1] = Schedule(0, i).NumericValue();
  }
  int Nass = b->get_nassets();
  CellMatrix payoffs(NP, m);
  b->update_seed(seed);
  for (int i = 0; i < NP; i++) {
    std::vector<int> aishares(ishares.size());
    aishares = ishares;
    b->get_spots(spots);
    for (int j = 0; j < m; j++) {

      payoffs(i, j) =
          Notional * b->sim_basket_with_removal(
                         sched[j], sched[j + 1] - sched[j], aishares, barriers,
                         coupons, initial, spots, Nsteps, called);
    }
  }
  return payoffs;
}
} // namespace velesquant