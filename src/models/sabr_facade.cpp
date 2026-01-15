#include <cmath>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/object_cache.h>
#include <velesquant/models/sabr_facade.h>
#include <velesquant/models/utility.h>

namespace velesquant {

using namespace std;

std::string CreateSABRObj(const std::string &theObjName, double maturity,
                          double forward, double beta, double alpha, double nu,
                          double rho) {
  return placeObject(theObjName,
                     new Sabr(maturity, forward, beta, alpha, nu, rho));
}

std::string dvlCreateSABRObjShifted(const std::string &theObjName,
                                    double maturity, double forward,
                                    double beta, double alpha, double nu,
                                    double rho, double shift) {
  return placeObject(theObjName,
                     new Sabr(maturity, forward, beta, alpha, nu, rho, shift));
}

double sabrVolatility(double strike, double forward, double maturity,
                      double alpha, double beta, double nu, double rho) {
  // basic validation
  if (!((alpha > 0.0) && (beta >= 0.0 && beta <= 1.0) && (nu >= 0.0) &&
        (rho * rho <= 1.0)))
    throw std::runtime_error("SABR parameters are not valid");

  const double mean = std::pow(forward * strike, 0.5 * (1.0 - beta));
  const double logR = std::log(forward / strike);
  const double tem = (1.0 - beta) * (1.0 - beta) * logR * logR;
  const double dem = (1.0 + tem / 24.0 + tem * tem / 1920.0) * mean;
  const double coe =
      1.0 + maturity * ((1.0 - beta) * (1.0 - beta) * alpha * alpha /
                            (mean * mean) / 24.0 +
                        rho * beta * nu * alpha / mean / 4.0 +
                        (2.0 - 3.0 * rho * rho) * nu * nu / 24.0);
  const double z = (nu / alpha) * mean * logR;
  double multiplier;
  if (std::fabs(z * z) > 1.0E-20) {
    const double xz = std::log(
        (std::sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho));
    multiplier = z / xz;
  } else {
    multiplier = 1.0 - 0.5 * rho * z - (3.0 * rho * rho - 2.0) * z * z / 12.0;
  }
  return (alpha / dem) * multiplier * coe;
}

Matrix sabrCalibrator(const std::string &theName, Matrix theStrikes,
                      Matrix theQuotes, std::string quoteType) {
  Sabr *s = getObject<Sabr>(theName);
  int m = theStrikes.rows();
  int k = 0;
  std::vector<double> strikes(m);
  std::vector<double> quotes(m);
  for (int i = 0; i < m; i++) {
    if (!std::isnan(theStrikes(i, 0)) && !std::isnan(theQuotes(i, 0))) {
      k++;
      strikes[k - 1] = theStrikes(i, 0);
      quotes[k - 1] = theQuotes(i, 0);
    };
  }
  if (k == 0)
    throw std::runtime_error("Not enough data"); // VEL_RAISE
  strikes.resize(k);
  quotes.resize(k);
  s->calibrator(strikes, quotes, stringToCalibrationTarget(quoteType));
  Matrix sabrPara(3, 1);
  sabrPara(0, 0) = s->getParameterAlpha();
  sabrPara(1, 0) = s->getParameterNu();
  sabrPara(2, 0) = s->getParameterRho();
  return sabrPara;
}

Matrix sabrCalibratorWithInitial(const std::string &theName, Matrix theStrikes,
                                 Matrix theQuotes, std::string quoteType,
                                 std::string atmFlag) {
  Sabr *s = getObject<Sabr>(theName);
  int m = theStrikes.rows();

  std::vector<double> strikes(m);
  std::vector<double> quotes(m);
  double forward = s->getForward();
  double T = s->getMaturity();

  for (int i = 0; i < m; i++) {
    strikes[i] = theStrikes(i, 0);
    quotes[i] = theQuotes(i, 0);

    // get ATM vol
    if (std::abs(strikes[i] / forward - 1.0) < 1e-8) { // float compare
      if (quoteType == "premium" && atmFlag == "Yes") {
        double atmVol = implied_vol(T, forward, strikes[i], quotes[i]);
        s->setATMvol(atmVol);
      } else {
        s->setATMvol(quotes[i]);
      }
    }
  }

  if (atmFlag == "Yes") {
    s->calibratorWithInitialATM(strikes, quotes,
                                stringToCalibrationTarget(quoteType));
  } else {
    s->calibratorWithInitial(strikes, quotes,
                             stringToCalibrationTarget(quoteType));
  }

  Matrix sabrPara(3, 1);
  sabrPara(0, 0) = s->getParameterAlpha();
  sabrPara(1, 0) = s->getParameterNu();
  sabrPara(2, 0) = s->getParameterRho();
  return sabrPara;
}

} // namespace velesquant
