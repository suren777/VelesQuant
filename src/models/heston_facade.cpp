#include <velesquant/volatility/s_vol.h>
#include <velesquant/models/heston_facade.h>
#include <velesquant/models/object_cache.h>
#include <velesquant/models/utility.h>

namespace velesquant {

using namespace std;

std::string CreateHestonObj(const std::string &theObjName, double spot,
                            double var0, double kappa, double theta, double xi,
                            double rho, int seed) {
  return placeObject(theObjName,
                     new sVol(spot, var0, kappa, theta, xi, rho, seed));
}

double hestonPrice(std::string theName, double maturity, double forward,
                   double strike) {
  return getObject<sVol>(theName)->hestonPriceCF(maturity, forward, strike,
                                                 "call");
}

Matrix hestonCalibrator(const std::string &theName, Matrix theMaturitys,
                        Matrix theForwards, Matrix theStrikes, Matrix theQuotes,
                        string quoteType) {
  sVol *h = getObject<sVol>(theName);
  int m = theStrikes.rows();
  int k = theQuotes.rows();
  int l = theQuotes.cols();
  std::vector<double> maturitys, forwards, strikes, quotes;
  if ((m == k) && (l == 1)) {
    for (int i = 0; i < m; i++) {
      maturitys.push_back(theMaturitys(i, 0));
      forwards.push_back(theForwards(i, 0));
      strikes.push_back(theStrikes(i, 0));
      quotes.push_back(theQuotes(i, 0));
    }
  } else {
    m = theStrikes.cols();
    if (m != l)
      throw std::runtime_error("Error! Strikes dont match vols");

    for (int i = 0; i < k; i++)
      for (int j = 0; j < l; j++)
        if (theQuotes(i, j) > 1e-10) {
          maturitys.push_back(theMaturitys(i, 0));
          forwards.push_back(theForwards(i, 0));
          strikes.push_back(theStrikes(0, j));
          quotes.push_back(theQuotes(i, j));
        }
  }
  h->calibrator(maturitys, forwards, strikes, quotes, quoteType);
  Matrix hestonPara(5, 1);
  hestonPara(0, 0) = h->getParameterVar0();
  hestonPara(1, 0) = h->getParameterKappa();
  hestonPara(2, 0) = h->getParameterTheta();
  hestonPara(3, 0) = h->getParameterXi();
  hestonPara(4, 0) = h->getParameterRho();
  return hestonPara;
}

Matrix hestonSimulation(const std::string &theName, Matrix theTimes,
                        Matrix theForwards, int NP) {
  sVol *h = getObject<sVol>(theName);

  int m = theTimes.rows();
  std::vector<double> times(m), forwards(m);
  for (int i = 0; i < m; i++) {
    times[i] = theTimes(i, 0);
    forwards[i] = theForwards(i, 0);
  }
  Matrix thePaths(m, NP);
  for (int j = 0; j < NP; j++) {
    vector<double> pathF = h->simulationHeston(times, forwards);
    for (int i = 0; i < m; i++)
      thePaths(i, j) = pathF[i];
  }
  return thePaths;
}

// Additional Heston exports if needed (DNT, Cliquet etc)
// For now adhering to declarations in existing model.h as a guide for what's
// public. model.h listed only: hestonCalibrator, hestonSimulation, hestonPrice,
// CreateHestonObj So sticking to those for now to minimize surface area unless
// user requested specific extras.

} // namespace velesquant
