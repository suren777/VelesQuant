
#include <cmath>
#include <cstdlib>
#include <velesquant/models/black_scholes.h>
#include <velesquant/models/utility.h>
using namespace std;
#include <algorithm>

namespace velesquant {

double BlackScholesCall(double Spot, double Strike, double r, double d,
                        double Vol, double Expiry) {
  double standardDeviation = Vol * sqrt(Expiry);
  double d1 = (log(Spot / Strike) + (r - d) * Expiry) / standardDeviation +
              0.5 * standardDeviation;
  double d2 = d1 - standardDeviation;
  return Spot * exp(-d * Expiry) * cdf_normal(d1) -
         Strike * exp(-r * Expiry) * cdf_normal(d2);
}

double BlackScholesCallVega(double Spot, double Strike, double r, double d,
                            double Vol, double Expiry) {
  double standardDeviation = Vol * sqrt(Expiry);
  double d1 = (log(Spot / Strike) + (r - d) * Expiry) / standardDeviation +
              0.5 * standardDeviation;
  return Spot * exp(-d * Expiry) * sqrt(Expiry) * pdf_normal(d1);
}
} // namespace velesquant