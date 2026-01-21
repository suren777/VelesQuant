
#include <velesquant/models/utility.h>

#include <algorithm>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>
#include <ql/quantlib.hpp>
#include <string>
#include <vector>
#include <velesquant/models/black_formula.h>
#include <velesquant/models/solver.h>
#include <velesquant/numerics/interpolation.h>

// using namespace std;
// using namespace QuantLib;

namespace velesquant {

boost::mt19937 rng;
boost::normal_distribution<> nd(0.0, 1.0);
double random_normal() { return nd(rng); }

void set_seed(int seed) { rng.seed(seed); }

boost::math::normal_distribution<> normal(0, 1);
double cdf_normal(double p) { return cdf(normal, p); }

double pdf_normal(double p) { return pdf(normal, p); }

double implied_vol(double Maturity, double Forward, double Strike,
                   double Price) {
  BlackCall theCall(Maturity, Forward, Strike);
  return myBisection(Price, 1.0E-12, 1.0E+2, 1.0E-12, theCall);
};

double implied_vol2(double Maturity, double Forward, double Strike,
                    double Price) {
  BlackCall theCall(Maturity, Forward, Strike);
  return myNewtonRaphson<BlackCall, &BlackCall::Value, &BlackCall::Vega>(
      Price, 1.0E-1, 1.0E-5, theCall);
};

double annuity(double Expiry, double Tenor, double Freq,
               const std::vector<double> &DF, const std::vector<double> &T) {
  double ti, DFi;
  double annuity = 0.0;

  for (double i = 0; i < Tenor / Freq; i++) {
    ti = Expiry + (i + 1.0) * Freq;
    DFi = interpolation("LinearInterpolation", T, DF, ti);

    annuity = annuity + Freq * DFi;
  }

  return annuity;
};

double fwdSR(double Expiry, double Tenor, double Freq,
             const std::vector<double> &DF, const std::vector<double> &T) {
  double D0, DN;

  D0 = interpolation("LinearInterpolation", T, DF, Expiry);
  DN = interpolation("LinearInterpolation", T, DF, Expiry + Tenor);

  return (D0 - DN) / annuity(Expiry, Tenor, Freq, DF, T);
};

std::vector<std::vector<double>>
cholesky(const std::vector<std::vector<double>> &corrMatrix) {
  int N = int(corrMatrix.size());
  QuantLib::Matrix AAA(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      AAA[i][j] = corrMatrix[i][j];
  QuantLib::Matrix BBB(N, N);
  BBB = QuantLib::CholeskyDecomposition(AAA);
  std::vector<std::vector<double>> matrixCholesky(N, std::vector<double>(N));
  // matrixCholesky.resize(N);
  for (int i = 0; i < N; i++) {
    //		matrixCholesky[i].resize(N);
    for (int j = 0; j < N; j++)
      matrixCholesky[i][j] = BBB[i][j];
  }
  return matrixCholesky;
};

double option_vega(double F, double K, double sigma, double T) {
  double sig = sigma * std::sqrt(T);
  double d1 = (std::log(F / K) + sig * sig) / sig;
  return F * pdf_normal(d1);
};

void lAxy(const std::vector<std::vector<double>> &lowerTriangular,
          const std::vector<double> &input, std::vector<double> &result) {
  int n = input.size();
  for (int i = 0; i < n; i++) {
    result[i] = 0;
    for (int j = 0; j <= i; j++)
      result[i] += lowerTriangular[i][j] * input[j];
  }
};

// Converts DF curve to spot rates
std::vector<double> DFtoR(const std::vector<double> &DF,
                          const std::vector<double> &T, int n) {
  int m = DF.size();
  std::vector<double> R0s(m);

  for (int i = 0; i < m; i++) {
    if (n == 0) {
      R0s[i] = -std::log(DF[i]) / T[i];
    } else
      R0s[i] = n * (std::pow(1.0 / DF[i], 1 / (n * T[i])) - 1.0);
  }

  return R0s;
};

// Calculates Div yield based on forward prices, spot and DF curve
std::vector<double> DFtoDiv(const std::vector<double> &DF,
                            const std::vector<double> &dfT,
                            const std::vector<double> &Fwd,
                            const std::vector<double> &fwdT, double Spot,
                            int n) {
  int m = Fwd.size();
  std::vector<double> Divs(m), R0s(m);

  double interpR0;

  R0s = DFtoR(DF, dfT, n);
  QuantLib::LinearInterpolation interp(dfT.begin(), dfT.end(), R0s.begin());

  for (int i = 0; i < m; i++) {
    interpR0 = interp(fwdT[i]);
    Divs[i] = interpR0 - std::log(Fwd[i] / Spot) / fwdT[i];
  }

  return Divs;
};

// Calculates forward rates based on the DF curve
std::vector<double> DFtoFwd(const std::vector<double> &DF,
                            const std::vector<double> &T, double delta, int n) {
  int m = DF.size();
  std::vector<double> Fwds(m), R1(m);
  double R2, F1, F2;

  R1 = DFtoR(DF, T, n);
  QuantLib::LinearInterpolation interp(T.begin(), T.end(), R1.begin());
  bool allowExtrapolation = true;

  for (int i = 0; i < m; i++) {
    R2 = interp(T[i] + delta, allowExtrapolation);

    if (n == 0) {
      Fwds[i] = (R2 * (T[i] + delta) - R1[i] * T[i]) / delta;
    } else {
      F1 = std::pow(1.0 + R1[i] / n, n * T[i]);
      F2 = std::pow(1.0 + R2 / n, n * (T[i] + delta));
      Fwds[i] = n * (std::pow(F2 / F1, 1.0 / (n * delta)) - 1.0);
    }
  }

  return Fwds;
};

// Calculates the forward price
double FwdPrice(double Spot, double R0, double Div, double T) {
  return Spot * std::exp((R0 - Div) * T);
};

double asinh(double a) { return std::log(a + std::sqrt(a * a + 1.0)); };

// Calculates the log return volatilty for a MC simulation.
// This should converge towards the implied volatility.
std::vector<double>
MCimpVol(double freq, const std::vector<double> &theTimes,
         const std::vector<std::vector<double>> &thelogReturns) {
  int m = thelogReturns.size();
  int n = thelogReturns[0].size();

  std::vector<double> var(m), impVol(m);
  double dt;

  // Initial time step
  dt = theTimes[0];

  // calculate log return standard deviation
  double sum =
      std::accumulate(thelogReturns[0].begin(), thelogReturns[0].end(), 0.0);
  double mean = sum / n;
  double sq_sum =
      std::inner_product(thelogReturns[0].begin(), thelogReturns[0].end(),
                         thelogReturns[0].begin(), 0.0);
  double stDev = std::sqrt(sq_sum / n - mean * mean);

  // calculates implied volatility
  double fwdVol = stDev * std::sqrt(freq);

  var[0] = fwdVol * fwdVol * dt;
  impVol[0] = std::sqrt(var[0] / theTimes[0]);

  // Remaining time steps
  for (int j = 1; j < m; j++) {
    dt = theTimes[j] - theTimes[j - 1];

    // calculate log return standard deviation
    double sum2 =
        std::accumulate(thelogReturns[j].begin(), thelogReturns[j].end(), 0.0);
    double mean2 = sum2 / n;
    double sq_sum2 =
        std::inner_product(thelogReturns[j].begin(), thelogReturns[j].end(),
                           thelogReturns[j].begin(), 0.0);
    double stDev2 = std::sqrt(sq_sum2 / n - mean2 * mean2);

    // calculates implied volatility
    double fwdVol2 = stDev2 * std::sqrt(freq);

    var[j] = var[j - 1] + fwdVol2 * fwdVol2 * dt;
    impVol[j] = std::sqrt(var[j] / theTimes[j]);
  }

  return impVol;
};

OptionType stringToOptionType(const std::string &type) {
  if (type == "put" || type == "Put" || type == "p" || type == "P") {
    return OptionType::Put;
  }
  return OptionType::Call;
}

CalibrationTarget stringToCalibrationTarget(const std::string &type) {
  if (type == "premium" || type == "Premium" || type == "price" ||
      type == "Price") {
    return CalibrationTarget::Price;
  } else if (type == "impliedVol" || type == "ImpliedVol" || type == "vol" ||
             type == "Vol" || type == "volatility" || type == "Volatility") {
    return CalibrationTarget::Volatility;
  } else {
    return CalibrationTarget::Price;
  }
}

} // namespace velesquant