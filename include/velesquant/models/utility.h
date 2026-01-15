//      utility.h

#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>

namespace velesquant {

// using namespace std; // REMOVED by refactor

void set_seed(int seed);

double random_normal();

double cdf_normal(double p);
double pdf_normal(double p);
double option_vega(double F, double K, double sigma, double T);

double implied_vol(double Maturity, double Forward, double Strike,
                   double Price);

std::vector<std::vector<double>>
cholesky(const std::vector<std::vector<double>> &corrMatrix);

void lAxy(const std::vector<std::vector<double>> &A,
          const std::vector<double> &x, std::vector<double> &y);

std::vector<double> DFtoR(const std::vector<double> &DF,
                          const std::vector<double> &T, int n);
std::vector<double> DFtoDiv(const std::vector<double> &DF,
                            const std::vector<double> &dfT,
                            const std::vector<double> &Fwd,
                            const std::vector<double> &fwdT, double Spot,
                            int n);
std::vector<double> DFtoFwd(const std::vector<double> &DF,
                            const std::vector<double> &T, double delta, int n);
double FwdPrice(double Spot, double R0, double Div, double T);
double asinh(double);

std::vector<double>
MCimpVol(double freq, const std::vector<double> &theTimes,
         const std::vector<std::vector<double>> &thelogReturns);

double annuity(double Expiry, double Tenor, double Freq,
               const std::vector<double> &DF, const std::vector<double> &T);
double fwdSR(double Expiry, double Tenor, double Freq,
             const std::vector<double> &DF, const std::vector<double> &T);

struct defSwap {
  double Expiry;
  double Tenor;
  double Frequency;
  double SwapRate;
  double VolATM;
  double Value;
};

struct volQuote {
  double Expiry;
  double Strike;
  double IV;
};

enum class OptionType { Call, Put };

OptionType stringToOptionType(const std::string &type);

enum class CalibrationTarget { Price, Volatility };

CalibrationTarget stringToCalibrationTarget(const std::string &type);

} // namespace velesquant
#endif