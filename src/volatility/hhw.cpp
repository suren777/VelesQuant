#include <ql/quantlib.hpp>
#include <velesquant/models/utility.h>
#include <velesquant/volatility/hhw.h>
#include <velesquant/volatility/lm.h>
using namespace std;
#include <algorithm>
#include <cmath>

const double PI = 3.14159265358979323846264338327950288;
using namespace QuantLib;

namespace velesquant {

HHW::HHW(double s0, double v0, double r0, double kappa, double eta, double rho,
         double sigma1, double sigma2, double a)
    : s0_(s0), r0_(r0), v0_(v0), kappa_(kappa), eta_(eta), rho_(rho),
      sigma1_(sigma1), sigma2_(sigma2), a_(a) {};

double HHW::ib(const double &t, const double &T) const {
  return T / 2 * b(t * T / 2 + T / 2) *
         (1 - exp(-a_ * T - (t * T / 2 + T / 2)));
};

double HHW::intb(const double &T) const {
  std::function<double(double)> myInt = [this, T](double t) {
    return ib(t, T);
  };
  QuantLib::GaussChebyshevIntegration gChebInt(64);
  return gChebInt(myInt);
};

double HHW::HHWIntegrand(double y, double T, double K, int type) const {
  double delta = (type == 1) ? 0 : 1;
  Cdoub li(0, 1);
  Cdoub H = (li * y - delta) / a_ * (1 - exp(-T * a_));
  double alpha = kappa_ * eta_;
  double beta = (type == 1) ? (kappa_ - rho_ * sigma1_) : kappa_;
  double gamma = (type == 1) ? 0.5 : -0.5;
  double sigma12 = (sigma1_ * sigma1_);
  Cdoub brs = beta - li * rho_ * sigma1_ * y;
  Cdoub d = std::sqrt(brs * brs - sigma12 * (2.0 * li * gamma * y - y * y));
  Cdoub g = (brs + d) / (brs - d);
  Cdoub G = (brs + d) / sigma12 * ((1.0 - exp(d * T)) / (1.0 - g * exp(d * T)));
  Cdoub F =
      alpha / sigma12 *
          ((brs + d) * T - 2.0 * log((1.0 - g * exp(d * T)) / (1.0 - g))) +
      (li * y - delta) * intb(T) +
      sigma2_ * sigma2_ / 2 * (li * y - delta) * (li * y - delta) / (a_ * a_) *
          (T + 2 / a_ * exp(-a_ * T) - .5 / a_ * exp(-2 * a_ * T) -
           3 / 2.0 / a_);
  Cdoub pre;
  if (type == 1)
    pre = exp(F + G * v0_ + H * r0_ + li * log(s0_) * y);
  else
    pre = exp(F + G * v0_ + H * r0_ + li * log(s0_) * y - c(T));

  return real(exp(-li * log(K)) * pre / (li * y));
};
double HHW::c(double T) const {
  return -r0_ / a_ * (1 - exp(-a_ * T)) - intb(T) +
         sigma2_ * sigma2_ / 2 *
             (T + 2 / a_ * exp(-a_ * T) - .5 / a_ * exp(-2 * a_ * T) -
              3 / 2.0 / a_);
};

double HHW::HHWPrice(double maturity, double strike) const {
  std::function<double(double)> fcn1 = [this, maturity, strike](double y) {
    return HHWIntegrand(y, maturity, strike, 1);
  };
  std::function<double(double)> fcn2 = [this, maturity, strike](double y) {
    return HHWIntegrand(y, maturity, strike, 2);
  };
  GaussLaguerreIntegration gLegInt(16);
  double integr1 = gLegInt(fcn1);
  double integr2 = gLegInt(fcn2);
  return s0_ * (.5 + integr1 / PI) -
         strike * exp(c(maturity)) * (.5 + integr2 / PI);
}
} // namespace velesquant