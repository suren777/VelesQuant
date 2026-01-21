
#include <algorithm>
#include <cmath>
#include <complex>
#include <ql/quantlib.hpp>
#include <velesquant/models/utility.h>
#include <velesquant/volatility/lm.h>
#include <velesquant/volatility/svolt.h>

// using namespace std;
// using namespace QuantLib;

namespace velesquant {

typedef std::vector<double> Vdoub;
typedef std::complex<double> Cdoub;

const double PI = 3.14159265358979323846264338327950288;
const double DT = 0.0001;
const double SDT = std::sqrt(DT);

double sVolT::interpolateF(std::vector<double> T, std::vector<double> F,
                           double t) const {
  if (t < T.front())
    return spot_ + (F.front() - spot_) / T.front() * t;

  auto it = std::upper_bound(T.begin(), T.end(), t);
  if (it == T.end()) {
    size_t n = T.size();
    return F[n - 2] +
           (t - T[n - 2]) / (T[n - 1] - T[n - 2]) * (F[n - 1] - F[n - 2]);
  }

  size_t i = std::distance(T.begin(), it);
  return F[i - 1] + (F[i] - F[i - 1]) / (T[i] - T[i - 1]) * (t - T[i - 1]);
}

sVolT::sVolT(double spot, double var0, double kappa, double rho, Vdoub &time,
             Vdoub &theta, Vdoub &xit, int seed)
    : spot_(spot), var0_(var0), kappa_(kappa), rho_(rho) {
  if (time.size() == theta.size() && theta.size() == xit.size()) {
    createT(time);
    createXitT(xit);
    createthetat(theta);
    model_ = full;
  } else
    throw("Wrong Input Parameters!!!\n");
  if (seed == 0)
    update_seed(std::time(NULL));
  else if (seed > 0)
    update_seed(seed);
};

sVolT::sVolT(double spot, double var0, double kappa, double rho, Vdoub &time,
             Vdoub &theta, double xit, int seed)
    : spot_(spot), var0_(var0), kappa_(kappa), rho_(rho) {
  if (time.size() == theta.size()) {
    createT(time);
    createXitT(xit);
    createthetat(theta);
    model_ = meanrev;
  } else
    throw("Wrong Input Parameters!!!\n");
  if (seed == 0)
    update_seed(std::time(NULL));
  else if (seed > 0)
    update_seed(seed);
};

sVolT::sVolT(double spot, double var0, double kappa, double rho, Vdoub &time,
             double theta, Vdoub &xit, int seed)
    : spot_(spot), var0_(var0), kappa_(kappa), rho_(rho) {
  if (time.size() == xit.size()) {
    createT(time);
    createXitT(xit);
    createthetat(theta);
    model_ = volvol;
  } else
    throw("Wrong Input Parameters!!!\n");
  if (seed == 0)
    update_seed(std::time(NULL));
  else if (seed > 0)
    update_seed(seed);
};

double sVolT::trapz(Vdoub x, Vdoub y) const {
  int n = x.size();
  double sum = 0;
  for (int i = 1; i <= n - 1; i++)
    sum += 0.5 * (x[i] - x[i - 1]) * (y[i - 1] + y[i]);
  return sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double sVolT::hestonIntegrandTd(Cdoub k, double maturity, double forward,
                                double strike, int choice) const {
  Cdoub i(0.0, 1.0);
  Cdoub a, b;
  Cdoub d, r1, r2, g, Y, F1(0, 0), F2(0, 0), gn;
  double I, xi, xi2, dT, theta;
  int fl = 0;

  if (maturity > t_[nT_ - 1])
    throw("Maturity is out of range");

  if (choice == 0) {
    a = 0.5 * k * (k + i);
    I = 0;
  } else {
    a = 0.5 * k * (k - i); // OK
    I = 1.0;
  }

  for (int j = 0; j < nT_; j++) {
    xi = xit_[j];
    theta = thetat_[j];
    xi2 = xi * xi;
    if (j == 0)
      dT = ((t_[j] < maturity) ? t_[j] : maturity);
    else
      dT = ((t_[j] < maturity) ? (t_[j] - t_[j - 1]) : (maturity - t_[j - 1]));

    fl = (t_[j] >= maturity) ? 1 : 0;
    b = kappa_ - rho_ * xi * (k * i + I);
    d = std::sqrt(b * b + 2.0 * a * xi2); // change
    r1 = (b + d) / (xi * xi);             //
    r2 = (b - d) / (xi * xi);             //
    g = r2 / r1;
    gn = (b - d - F2 * xi2) / (b + d - F2 * xi2);
    Y = std::exp(-d * dT);
    F1 += kappa_ * theta / xi2 *
          ((b - d) * dT - 2.0 * std::log((1.0 - gn * Y) / (1.0 - gn)));
    F2 = r1 * (g - gn * Y) / (1.0 - gn * Y);
    if (fl == 1)
      break;
  }

  Cdoub CF = std::exp(F1 + F2 * var0_);
  Cdoub integrand = std::exp(i * k * std::log(forward / strike)) * CF / (k * i);
  if (std::real(integrand) != std::real(integrand))
    return 0.0;
  else
    return std::real(integrand);
}

std::complex<double> sVolT::hestonCFTd(std::complex<double> k,
                                       double maturity) const {
  Cdoub i(0.0, 1.0);
  Cdoub a, b;
  Cdoub d, r1, r2, g, Y, F1(0, 0), F2(0, 0), gn;
  double I, xi, xi2, dT, theta;
  int fl = 0;

  if (maturity > t_[nT_ - 1])
    throw("Maturity is out of range");
  a = 0.5 * k * (k + i);
  I = 0;
  for (int j = 0; j < nT_; j++) {
    xi = xit_[j];
    theta = thetat_[j];
    xi2 = xi * xi;
    if (j == 0)
      dT = ((t_[j] < maturity) ? t_[j] : maturity);
    else
      dT = ((t_[j] < maturity) ? (t_[j] - t_[j - 1]) : (maturity - t_[j - 1]));

    fl = (t_[j] >= maturity) ? 1 : 0;
    b = kappa_ - rho_ * xi * (k * i + I);
    d = std::sqrt(b * b + 2.0 * a * xi2); // change
    r1 = (b + d) / (xi * xi);             //
    r2 = (b - d) / (xi * xi);             //
    g = r2 / r1;
    gn = (b - d - F2 * xi2) / (b + d - F2 * xi2);
    Y = std::exp(-d * dT);
    F1 += kappa_ * theta / xi2 *
          ((b - d) * dT - 2.0 * std::log((1.0 - gn * Y) / (1.0 - gn)));
    F2 = r1 * (g - gn * Y) / (1.0 - gn * Y);
    if (fl == 1)
      break;
  }

  Cdoub CF = std::exp(F1 + F2 * var0_);
  if (CF != CF)
    return 0.0;
  else
    return CF;
}

#ifdef _MSC_VER
#endif
double sVolT::hestonPriceTd(double maturity, double forward, double strike,
                            OptionType optType) const {
  auto fcn1 = [this, maturity, forward, strike](double k) {
    return this->hestonIntegrandTd(k, maturity, forward, strike, 0);
  };
  auto fcn2 = [this, maturity, forward, strike](double k) {
    return this->hestonIntegrandTd(k, maturity, forward, strike, 1);
  };
  QuantLib::GaussLaguerreIntegration gLegInt(128);
  double integr1 = gLegInt(fcn1);
  double integr2 = gLegInt(fcn2);
  double callPrice =
      forward * (0.5 + integr2 / PI) - strike * (0.5 + integr1 / PI);
  if (optType == OptionType::Call)
    return callPrice;
  else
    return callPrice + strike - forward;
}

double sVolT::intCFfun(double u, double ki, double X, double maturity) const {
  Cdoub ret =
      std::exp(Cdoub(0.0, 1.0) * u * X) * hestonCFTd(Cdoub(u, ki), maturity);
  return ret.real() / (u * u + ki * ki);
}

#ifdef _MSC_VER
#endif

double sVolT::hestonPriceTdCF(double maturity, double forward, double strike,
                              OptionType optType) const {
  if (strike == 0.0)
    return forward;
  double X = log(forward / strike);
  double ki = -0.5;
  auto fcn1 = [this, ki, X, maturity](double u) {
    return this->intCFfun(u, ki, X, maturity);
  };
  QuantLib::GaussLaguerreIntegration gLegInt(64);
  double integr1 = gLegInt(fcn1);
  double callPrice = forward - sqrt(strike * forward) * integr1 / PI;
  if (callPrice <= 0.0)
    return 0;
  if (optType == OptionType::Call)
    return callPrice;
  else
    return callPrice + strike - forward;
}

void sVolT::objFcnTd(int m, int /*n*/, double *x, double *fvec,
                     int * /*iflag*/) {
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterRho(x[2]);
  switch (model_) {
  case 0:
    for (int i = 0; i < nT_; i++)
      setParameterTheta(x[3], i);
    for (int i = 0; i < nT_; i++)
      setParameterXi(x[4 + i], i);
    break;
  case 1:
    for (int i = 0; i < nT_; i++)
      setParameterTheta(x[3 + i], i);
    for (int i = 0; i < nT_; i++)
      setParameterXi(x[3 + nT_], i);
    break;
  case 2:
    for (int i = 0; i < nT_; i++)
      setParameterTheta(x[3 + i], i);
    for (int i = 0; i < nT_; i++)
      setParameterXi(x[3 + i + nT_], i);
  }
  for (int i = 0; i < m; i++)
    fvec[i] = (1.0 - hestonPriceTdCF(maturities_[i], forwards_[i], strikes_[i],
                                     OptionType::Call) /
                         marketQuotes_[i]) *
              100;
};

void sVolT::calibratorLM(Vdoub maturities, Vdoub forwards, Vdoub strikes,
                         Vdoub marketQuotes, CalibrationTarget quoteType) {
  int m = strikes.size(); // no. of observations
  maturities_.resize(m);
  maturities_ = maturities;
  forwards_.resize(m);
  forwards_ = forwards;
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  if (quoteType == CalibrationTarget::Volatility) {
    for (int i = 0; i < m; i++) {
      double vol = marketQuotes[i] * std::sqrt(maturities[i]);
      double d1 = std::log(forwards[i] / strikes[i]) / vol + 0.5 * vol;
      double d2 = d1 - vol;
      marketQuotes_[i] =
          forwards[i] * cdf_normal(d1) - strikes[i] * cdf_normal(d2);
    }
  }
  int n;
  if (model_ == 2)
    n = 3 + nT_ * 2;
  else
    n = 4 + nT_;
  // no. of HESTON model paremeters
  if (n > int(maturities.size()))
    throw("Need more data points!!!");
  std::vector<double> x(n); // initial estimate of parameters vector
  x[0] = getParameterVar0();
  x[1] = getParameterKappa();
  x[2] = getParameterRho();

  switch (model_) {
  case 0:
    x[3] = getParameterTheta(0);
    for (int i = 0; i < nT_; i++)
      x[4 + i] = getParameterXi(i);
    break;
  case 1:
    for (int i = 0; i < nT_; i++)
      x[3 + i] = getParameterTheta(i);
    x[3 + nT_] = getParameterXi(0);
    break;
  case 2:
    for (int i = 0; i < nT_; i++)
      x[3 + i] = getParameterTheta(i);
    for (int i = 0; i < nT_; i++)
      x[3 + i + nT_] = getParameterXi(i);
  }
  std::vector<double> fvec(m);     // no need to populate
  double ftol = 1e-08;             // tolerance
  double xtol = 1e-08;             // tolerance
  double gtol = 1e-08;             // tolerance
  int maxfev = 400;                // maximum function evaluations
  double epsfcn = 1e-08;           // tolerance
  std::vector<double> diag(n);     // some internal thing
  int mode = 1;                    // some internal thing
  double factor = 1;               // a default recommended value
  int nprint = 0;                  // don't know what it does
  int info = 0;                    // output variable
  int nfev = 0;                    // output variable
  std::vector<double> fjac(m * n); // output array of jacobian
  int ldfjac = m;                  // recommended setting
  std::vector<int> ipvt(n);        // for internal use
  std::vector<double> qtf(n);      // for internal use
  std::vector<double> wa1(n);      // for internal use
  std::vector<double> wa2(n);      // for internal use
  std::vector<double> wa3(n);      // for internal use
  std::vector<double> wa4(m);      // for internal use

  auto fcnheston = [this](int m, int n, double *x, double *fvec, int *iflag) {
    this->objFcnTd(m, n, x, fvec, iflag);
  };
  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcnheston);
  QL_ENSURE(info >= 1 && info <= 4,
            "Heston Model Calibration Fails: " << getLmdifMessage(info)
                                               << " (info=" << info << ")");

  // the below is output result
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterRho(x[2]);
  switch (model_) {
  case 0:
    for (int i = 0; i < nT_; i++)
      setParameterTheta(x[3], i);
    for (int i = 0; i < nT_; i++)
      setParameterXi(x[4 + i], i);
    break;
  case 1:
    for (int i = 0; i < nT_; i++)
      setParameterTheta(x[3 + i], i);
    for (int i = 0; i < nT_; i++)
      setParameterXi(x[3 + nT_], i);
    break;
  case 2:
    for (int i = 0; i < nT_; i++)
      setParameterTheta(x[3 + i], i);
    for (int i = 0; i < nT_; i++)
      setParameterXi(x[3 + i + nT_], i);
  }
};

Vdoub sVolT::simulationHestonTd(Vdoub times, Vdoub forwards) const {
  int T = 0, N = times.size(), j = 0;
  Vdoub pathF(N);

  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_), xin = xit_[0],
         theta = thetat_[0], accT = 0;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    if (accT > t_[j]) {
      j++;
      if (j >= nT_) {
        theta = thetat_[nT_ - 1];
        xin = xit_[nT_ - 1];
      } else {
        theta = thetat_[j];
        xin = xit_[j];
      }
    }
    accT += DT;
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * std::sqrt(std::abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta - V) * DT +
         xin * std::sqrt(std::abs(V)) * SDT * z2; // Euler scheme
    if ((i * DT <= times[T]) && ((i + 1) * DT > times[T])) {
      pathF[T] = forwards[T] * S / spot_;
      T++;
    }
  }
  return pathF;
};

Vdoub sVolT::simulationHestonTdMax(Vdoub times, Vdoub forwards) const {
  int T = 0, N = times.size(), j = 0;
  Vdoub pathF(N);

  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_), xin = xit_[0],
         theta = thetat_[0], accT = 0;
  double Smax = S;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    if (accT > t_[j]) {
      j++;
      if (j >= nT_) {
        theta = thetat_[nT_ - 1];
        xin = xit_[nT_ - 1];
      } else {
        theta = thetat_[j];
        xin = xit_[j];
      }
    }
    accT += DT;
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * std::sqrt(std::abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta - V) * DT +
         xin * std::sqrt(std::abs(V)) * SDT * z2; // Euler scheme
    Smax = (Smax < S) ? S : Smax;
    if ((i * DT <= times[T]) && ((i + 1) * DT > times[T])) {
      pathF[T] = forwards[T] * Smax / spot_;
      T++;
    }
  }
  return pathF;
};

Vdoub sVolT::simulationHestonTdCliq(Vdoub times, Vdoub forwards, double gcap,
                                    double gfloor, double lcap, double lfloor,
                                    double alpha) const {
  int T = 0, N = times.size(), j = 0;
  Vdoub pathF(N);

  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_), xin = xit_[0],
         theta = thetat_[0], accT = 0;
  double R, S1, S0 = S;
  int Day = 1;
  double dt = 1.0 / 365;
  double cRet = 0.0;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    if (accT > t_[j]) {
      j++;
      if (j >= nT_) {
        theta = thetat_[nT_ - 1];
        xin = xit_[nT_ - 1];
      } else {
        theta = thetat_[j];
        xin = xit_[j];
      }
    }
    accT += DT;
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * std::sqrt(std::abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta - V) * DT +
         xin * std::sqrt(std::abs(V)) * SDT * z2; // Euler scheme
    if (accT > Day * dt) {
      Day++;
      S1 = interpolateF(times, forwards, Day * dt) * S / spot_;
      R = S1 / (alpha * S0) - 1.0;
      cRet += std::max(lfloor, std::min(R, lcap));
      S0 = S1;
    }
    if ((i * DT <= times[T]) && ((i + 1) * DT > times[T])) {
      pathF[T] = std::min(std::max(cRet, gfloor), gcap);
      T++;
    }
  }
  return pathF;
};

std::vector<double> sVolT::simulationHestonTdCliq(
    std::vector<double> times, std::vector<double> /*forwards*/, double gcap,
    double gfloor, double lcap, double lfloor, double alpha, int tDay,
    Termstructure *mycurve) const {
  int T = 0, N = times.size(), j = 0;
  Vdoub pathF(N);

  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_), xin = xit_[0],
         theta = thetat_[0], accT = 0;
  double S0 = S;

  int Day = 1;
  double dt = 1.0 / 365;
  double cRet = 0.0;
  double r, q;
  r = mycurve->rate(tDay, times[0] / dt);
  q = mycurve->divident(tDay, times[0] / dt);
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    if (accT > t_[j]) {
      j++;
      if (j >= nT_) {
        theta = thetat_[nT_ - 1];
        xin = xit_[nT_ - 1];
      } else {
        theta = thetat_[j];
        xin = xit_[j];
      }
    }
    accT += DT;
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * std::sqrt(std::abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta - V) * DT +
         xin * std::sqrt(std::abs(V)) * SDT * z2; // Euler scheme
    if (accT > Day * dt) {
      Day++;
      cRet += std::max(lfloor, std::min(S / (alpha * S0) - 1.0, lcap));
      S0 = S;
    }
    if ((i * DT <= times[T]) && ((i + 1) * DT > times[T])) {
      pathF[T] = std::min(std::max(cRet, gfloor), gcap) *
                 std::exp((r - q) * (times[T + 1] - times[T]));
      /* S01 = S; */
      T++;
      r = mycurve->rate(tDay, (times[T] - times[T - 1]) / dt);
      q = mycurve->divident(tDay, (times[T] - times[T - 1]) / dt);
    }
  }
  return pathF;
};

void sVolT::update_seed(int seed) { set_seed(seed); };

/*
void sVolT::calibratorNM(Vdoub maturities, Vdoub forwards, Vdoub strikes,
Vdoub marketQuotes, string quoteType)
{
int m = strikes.size();		//no. of observations
maturities_.resize(m);
maturities_ = maturities;
forwards_.resize(m);
forwards_ = forwards;
strikes_.resize(m);
strikes_ = strikes;
marketQuotes_.resize(m);
marketQuotes_ = marketQuotes;
if(quoteType == "impliedVol")
{
for(int i=0; i<m; i++)
{
double vol = marketQuotes[i] * sqrt(maturities[i]);
double d1 = std::log(forwards[i]/strikes[i])/vol +0.5*vol;
double d2 = d1-vol;
marketQuotes_[i] = forwards[i]*cdf_normal(d1) - strikes[i]*cdf_normal(d2);
}
}

int n = 4+nT_;					//no. of HESTON model paremeters
Vdoub x(n);	//initial estimate of parameters vector
x[0] = getParameterVar0();
x[1] = getParameterKappa();
x[2] = getParameterTheta();
x[3] = getParameterRho();
for( int i = 0; i < nT_;i++) x[4+i]=getParameterXin(i);


double ftol = 1e-08; //tolerance
NMSimplex solve_me(ftol);

boost::function<double (sVol*, Vdoub)> objheston = &sVol::objFcnNM;
auto fcnheston = boost::bind(objheston, this, _1);
x = solve_me.minimize(x, .1, fcnheston);

//the below is output result
setParameterVar0(x[0]);
setParameterKappa(x[1]);
setParameterTheta(x[2]);
setParameterRho(x[3]);
for (int i = 0 ; i < nT_ ;i++)
xit_[i]=x[i+4];
};

double sVolT::objFcnNM(Vdoub x)
{
setParameterVar0(abs(x[0]));
setParameterKappa(x[1]);
setParameterTheta(x[2]);
setParameterRho(x[3]);
for (int i = 0 ; i < nT_ ;i++)
xit_[i]=x[i+4];
double sum = 0.0;
for(int i=0; i<(int)maturities_.size() ; i++)
{
double aux = 0.0;
aux = hestonPrice1(maturities_[i], forwards_[i], strikes_[i], "call") -
marketQuotes_[i]; sum += aux*aux;
}
return sum;
};


*/

} // namespace velesquant
