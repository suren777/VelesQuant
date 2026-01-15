#include <algorithm>

#include <cmath>
#include <complex>
#include <ql/quantlib.hpp>
#include <velesquant/models/utility.h>
#include <velesquant/volatility/lm.h>
#include <velesquant/volatility/schobzhu.h>

using namespace std;
using namespace QuantLib;

namespace velesquant {

const double PI = 3.14159265358979323846264338327950288;
const double DT = 0.001;
const double SDT = sqrt(DT);

SchobelZhu::SchobelZhu(double spot, double sigma0, double kappa, double theta,
                       double xi, double rho)
    : s0_(spot), sigma0_(sigma0), kappa_(kappa), theta_(theta), xi_(xi),
      rho_(rho) {};

double SchobelZhu::SchobelIntegrand(double u, double maturity, double forward,
                                    double strike, OptionType type) const {
  Cdoub li(0, 1);
  Cdoub X = (type == OptionType::Call) ? li * u : (1.0 + li * u);

  double xi2 = xi_ * xi_;

  Cdoub a = X * (1 - 2 * kappa_ * rho_ / xi_ - X * (1 - rho_ * rho_));
  Cdoub b = X * kappa_ * theta_ * rho_ / xi_;
  Cdoub c = X * rho_ / xi_;
  Cdoub alpha = std::sqrt(xi2 * a + kappa_ * kappa_);
  Cdoub beta = (kappa_ - xi2 * c) / alpha;
  Cdoub gamma = kappa_ * kappa_ * theta_ - xi2 * b;

  Cdoub saT = std::sinh(alpha * maturity);
  Cdoub caT = std::cosh(alpha * maturity);
  Cdoub D =
      1.0 / xi2 * (kappa_ - alpha * (saT + beta * caT) / (caT + beta * saT));
  Cdoub B =
      1.0 / xi2 / alpha *
      ((kappa_ * theta_ * alpha - beta * gamma + gamma * (saT + beta * caT)) /
           (caT + beta * saT) -
       kappa_ * theta_ * alpha);
  Cdoub C =
      .5 * (kappa_ * maturity - log(caT + beta * saT)) +
      (.5 *
           (kappa_ * theta_ * alpha * kappa_ * theta_ * alpha - gamma * gamma) *
           (saT / (caT + beta * saT) - alpha * maturity) +
       (kappa_ * theta_ * alpha - beta * gamma) * gamma * (caT - 1.0) /
           (caT + beta * saT)) /
          (xi2 * alpha * alpha * alpha);
  Cdoub A = X * (sigma0_ * sigma0_ + xi2 * maturity) * rho_ / xi_;

  Cdoub Psi = exp(-.5 * A + .5 * D * sigma0_ * sigma0_ + B * sigma0_ + C);

  return std::real(exp(li * u * log(forward / strike)) * Psi / (li * u));
}

double SchobelZhu::SchobelPrice(double maturity, double forward,
                                double strike) const {
  // Using lambdas instead of boost::bind for modern C++
  auto Phi0 = [this, maturity, forward, strike](double u) {
    return this->SchobelIntegrand(u, maturity, forward, strike,
                                  OptionType::Call);
  };
  auto Phi1 = [this, maturity, forward, strike](double u) {
    return this->SchobelIntegrand(u, maturity, forward, strike,
                                  OptionType::Put);
  };
  GaussLaguerreIntegration gLegInt(16);
  double callPrice = s0_ * (0.5 + (1 / PI) * gLegInt(Phi1)) -
                     strike * (0.5 + (1 / PI) * gLegInt(Phi0));
  return callPrice;
}

Vdoub SchobelZhu::simulation(const Vdoub &times, const Vdoub &forwards) const {
  int N = times.size();
  vector<double> pathF(N);
  int T = 0;
  double V = sigma0_, S = s0_, rho2 = sqrt(1 - rho_ * rho_);
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * V * SDT * z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * DT + xi_ * SDT * z2; // Euler scheme
    if (i * DT <= times[T] && (i + 1) * DT > times[T]) {
      pathF[T] = forwards[T] * S / s0_;
      T++;
    }
  }
  return pathF;
}

void SchobelZhu::objFcn(int m, int /*n*/, double *x, double *fvec,
                        int * /*iflag*/) {
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterTheta(x[2]);
  setParameterXi(x[3]);
  setParameterRho(x[4]);
  double diff = 0.0;
  for (int i = 0; i < m; i++) {
    double modelPrice = SchobelPrice(maturitys_[i], forwards_[i], strikes_[i]);
    if (target_ == CalibrationTarget::Volatility) {
      double vol =
          implied_vol(maturitys_[i], forwards_[i], strikes_[i], modelPrice);
      diff = vol - marketQuotes_[i];
    } else {
      diff = modelPrice - marketQuotes_[i];
    }
    fvec[i] = diff;
  }
};

void SchobelZhu::calibrator(const Vdoub &maturitys, const Vdoub &forwards,
                            const Vdoub &strikes, const Vdoub &marketQuotes,
                            CalibrationTarget target) {
  target_ = target;
  int m = static_cast<int>(strikes.size()); // no. of observations
  maturitys_.resize(m);
  maturitys_ = maturitys;
  forwards_.resize(m);
  forwards_ = forwards;
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;

  int n = 5;                // no. of SZ model paremeters
  std::vector<double> x(n); // initial estimate of parameters vector
  x[0] = getParameterVar0();
  x[1] = getParameterKappa();
  x[2] = getParameterTheta();
  x[3] = getParameterXi();
  x[4] = getParameterRho();

  std::vector<double> fvec(m); // no need to populate
  double ftol = 1e-08;         // tolerance
  double xtol = 1e-08;         // tolerance
  double gtol = 1e-08;         // tolerance
  int maxfev = 400;            // maximum function evaluations
  double epsfcn = 1e-08;       // tolerance
  std::vector<double> diag(n); // some internal thing
  int mode = 1;                // some internal thing
  double factor = 1;           // a default recommended value
  int nprint = 0;              // don't know what it does
  int info = 0;                // output variable
  int nfev = 0; // output variable will store no. of function evals
  std::vector<double> fjac(m * n); // output array of jacobian
  int ldfjac = m;                  // recommended setting
  std::vector<int> ipvt(n);        // for internal use
  std::vector<double> qtf(n);      // for internal use
  std::vector<double> wa1(n);      // for internal use
  std::vector<double> wa2(n);      // for internal use
  std::vector<double> wa3(n);      // for internal use
  std::vector<double> wa4(m);      // for internal use

  auto fcnSZ = [this](int m, int n, double *x, double *fvec, int *iflag) {
    this->objFcn(m, n, x, fvec, iflag);
  };
  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcnSZ);

  QL_ENSURE(info >= 1 && info <= 4,
            "Schob Zhu Model Calibration Fails: " << getLmdifMessage(info)
                                                  << " (info=" << info << ")");
  // the below is output result
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterTheta(x[2]);
  setParameterXi(x[3]);
  setParameterRho(x[4]);
  // RAII: vectors automatically cleaned up
};
} // namespace velesquant