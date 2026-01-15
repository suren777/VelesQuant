
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <cmath>
#include <complex>
#include <ql/quantlib.hpp>
#include <velesquant/volatility/lm.h>
#include <velesquant/volatility/s_vol.h>
#include <velesquant/models/utility.h>

using namespace std;
using namespace QuantLib;
using namespace boost::placeholders;

namespace velesquant {

const double PI = 3.14159265358979323846264338327950288;
const double DT = .001;
const double SDT = sqrt(DT);

double sVol::interpolateF(vector<double> T, vector<double> F, double t) const {
  int i, N = T.size();
  if (t == 0.0)
    return spot_;
  for (i = 0; i < N; i++)
    if (t < T[i])
      break;
  if (i == 0)
    return spot_ + (F[0] - spot_) / T[0] * t;
  else if ((i > 0) && (i < N))
    return F[i - 1] + (F[i] - F[i - 1]) / (T[i] - T[i - 1]) * (t - T[i - 1]);
  else
    return F[N - 2] +
           (t - T[N - 2]) / (T[N - 1] - T[N - 2]) * (F[N - 1] - F[N - 2]);
}

vector<double> sVol::simulationHeston(vector<double> times,
                                      vector<double> forwards) const {
  int N = times.size();
  vector<double> pathF(N);
  int T = 0;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * DT + xi_ * sqrt(abs(V)) * SDT * z2 +
         0.25 * xi_ * xi_ * DT * (z2 * z2 - 1); // Milstein scheme
    // V += kappa_*(theta_-V)*DT +xi_*sqrt(abs(V))*SDT*z2;  //Euler scheme
    if (i * DT <= times[T] && (i + 1) * DT > times[T]) {
      pathF[T] = forwards[T] * S / spot_;
      T++;
    }
  }
  return pathF;
}

vector<double> sVol::simulationHestonDO(vector<double> times,
                                        vector<double> forwards,
                                        double barrier) const {
  int N = times.size();
  vector<double> pathF(2);
  int T = 0;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  double flag = 1.0;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    if (flag > 0.0) {
      S += S * sqrt(abs(V)) * SDT *
           z1; // This line (S) must put BEFORE line (V)!!!
      V += kappa_ * (theta_ - V) * DT + xi_ * sqrt(abs(V)) * SDT * z2 +
           0.25 * xi_ * xi_ * DT * (z2 * z2 - 1); // Milstein scheme
      // V += kappa_*(theta_-V)*DT +xi_*sqrt(abs(V))*SDT*z2;  //Euler scheme
      if ((i * DT <= times[T]) && ((i + 1) * DT > times[T]) && (S <= barrier))
        flag = 0.0;
    }
  }
  pathF[0] = forwards[N - 1] * S / spot_;
  pathF[1] = flag;
  return pathF;
}

void sVol::update_seed(int seed) { set_seed(seed); };

sVol::sVol(double spot, double var0, double kappa, double theta, double xi,
           double rho, int seed)
    : spot_(spot), var0_(var0), kappa_(kappa), theta_(theta), xi_(xi),
      rho_(rho) {
  if (seed == 0)
    update_seed(int(std::time(NULL)));
  else if (seed > 0)
    update_seed(seed);
};

double sVol::trapz(vector<double> x, vector<double> y) const {
  int n = x.size();
  double sum = 0;
  for (int i = 1; i <= n - 1; i++)
    sum += 0.5 * (x[i] - x[i - 1]) * (y[i - 1] + y[i]);
  return sum;
}

#ifdef _MSC_VER
#pragma warning(disable : 4244)
#endif

double sVol::hestonPrice(double maturity, double forward, double strike,
                         string optType) const {
  // if ((forward/strike>1e5)||(strike==0.0)) return forward-strike;
  if (strike == 0.0)
    return forward;
  else {
    boost::function<double(double)> fcn1, fcn2;

    fcn1 = boost::bind(&sVol::hestonIntegrand, this, _1, maturity, forward,
                       strike, 0);
    fcn2 = boost::bind(&sVol::hestonIntegrand, this, _1, maturity, forward,
                       strike, 1);
    GaussLaguerreIntegration gLegInt(32);
    double integr1 = gLegInt(fcn1);
    double integr2 = gLegInt(fcn2);
    double callPrice =
        forward * (0.5 + integr2 / PI) - strike * (0.5 + integr1 / PI);

    if (callPrice <= 0.0)
      return 0;

    if (optType == "Call" || optType == "call")
      return callPrice;
    else
      return callPrice + strike - forward;
  }
}

double sVol::hestonIntegrand(complex<double> k, double maturity, double forward,
                             double strike, int choice) const {
  complex<double> i(0.0, 1.0);
  complex<double> a, b;
  if (choice == 0) {
    a = 0.5 * k * (k + i);
    b = kappa_ - rho_ * xi_ * k * i;
  } else {
    a = 0.5 * k * (k - i);
    b = kappa_ - rho_ * xi_ * (k * i + 1.0);
  }
  complex<double> d = sqrt(b * b + 2.0 * a * xi_ * xi_);
  complex<double> r1 = (b + d) / (xi_ * xi_);
  complex<double> r2 = (b - d) / (xi_ * xi_);
  complex<double> g = r2 / r1;
  complex<double> Y = exp(-d * maturity);
  complex<double> F1 =
      kappa_ *
      (r2 * maturity - 2.0 / (xi_ * xi_) * log((1.0 - g * Y) / (1.0 - g)));
  complex<double> F2 = r2 * (1.0 - Y) / (1.0 - g * Y);
  complex<double> CF = exp(F1 * theta_ + F2 * var0_);
  complex<double> integrand = exp(i * k * log(forward / strike)) * CF / (k * i);
  if (real(integrand) != real(integrand))
    return 0.0;
  else
    return real(integrand);
}

complex<double> sVol::hestonCF(complex<double> k, double maturity) const {
  complex<double> i(0.0, 1.0);
  complex<double> a = 0.5 * k * (k + i);
  complex<double> b = kappa_ - rho_ * xi_ * k * i;
  complex<double> d = sqrt(b * b + 2.0 * a * xi_ * xi_);
  complex<double> r1 = (b + d) / (xi_ * xi_);
  complex<double> r2 = (b - d) / (xi_ * xi_);
  complex<double> g = r2 / r1;
  complex<double> Y = exp(-d * maturity);
  complex<double> F1 =
      kappa_ *
      (r2 * maturity - 2.0 / (xi_ * xi_) * log((1.0 - g * Y) / (1.0 - g)));
  complex<double> F2 = r2 * (1.0 - Y) / (1.0 - g * Y);
  complex<double> CF = exp(F1 * theta_ + F2 * var0_);
  if (CF != CF)
    return 0.0;
  else
    return CF;
}

double sVol::hestonPriceCF(double maturity, double forward, double strike,
                           string optType) const {
  if (strike == 0.0)
    return forward;
  double X = log(forward / strike);
  double ki = -0.5;
  boost::function<double(double)> fcn1;
  fcn1 = boost::bind(&sVol::intCFfun, this, _1, ki, X, maturity);
  GaussLaguerreIntegration gLegInt(32);
  double integr1 = gLegInt(fcn1);
  double callPrice = forward - sqrt(strike * forward) * integr1 / PI;
  if (callPrice <= 0.0)
    return 0;
  if (optType == "Call" || optType == "call")
    return callPrice;
  else
    return callPrice + strike - forward;
}

double sVol::intCFfun(double u, double ki, double X, double maturity) const {
  complex<double> ret = exp(complex<double>(0.0, 1.0) * u * X) *
                        hestonCF(complex<double>(u, ki), maturity);
  return ret.real() / (u * u + ki * ki);
}
void sVol::calibrator(vector<double> maturitys, vector<double> forwards,
                      vector<double> strikes, std::vector<double> marketQuotes,
                      string quoteType) {
  int m = strikes.size(); // no. of observations
  weights_.resize(m);
  for (int i = 0; i < m; i++) {
    weights_[i] = exp(-abs(spot_ - strikes[i]) / spot_ * 100 * maturitys[i]);
  }
  double sum = 0;
  for (int i = 0; i < m; i++)
    sum += weights_[i];
  for (int i = 0; i < m; i++)
    weights_[i] /= sum;
  maturitys_.resize(m);
  maturitys_ = maturitys;
  forwards_.resize(m);
  forwards_ = forwards;
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  if (quoteType == "impliedVol") {
    for (int i = 0; i < m; i++) {
      double vol = marketQuotes[i] * sqrt(maturitys[i]);
      double d1 = std::log(forwards[i] / strikes[i]) / vol + 0.5 * vol;
      double d2 = d1 - vol;
      marketQuotes_[i] =
          forwards[i] * cdf_normal(d1) - strikes[i] * cdf_normal(d2);
    }
  }

  int n = 5;                 // no. of HESTON model paremeters
  double *x = new double[n]; // initial estimate of parameters vector
  x[0] = getParameterVar0();
  x[1] = getParameterKappa();
  x[2] = getParameterTheta();
  x[3] = getParameterXi();
  x[4] = getParameterRho();

  double *fvec = new double[m]; // no need to populate
  double ftol = 1e-08;          // tolerance
  double xtol = 1e-08;          // tolerance
  double gtol = 1e-08;          // tolerance
  int maxfev = 400;             // maximum function evaluations
  double epsfcn = 1e-08;        // tolerance
  double *diag = new double[n]; // some internal thing
  int mode = 1;                 // some internal thing
  double factor = 1;            // a default recommended value
  int nprint = 0;               // don't know what it does
  int info = 0;                 // output variable
  int nfev = 0; // output variable will store no. of function evals
  double *fjac = new double[m * n]; // output array of jacobian
  int ldfjac = m;                   // recommended setting
  int *ipvt = new int[n];           // for internal use
  double *qtf = new double[n];      // for internal use
  double *wa1 = new double[n];      // for internal use
  double *wa2 = new double[n];      // for internal use
  double *wa3 = new double[n];      // for internal use
  double *wa4 = new double[m];      // for internal use

  boost::function<void(sVol *, int, int, double *, double *, int *)> objheston =
      &sVol::objFcn;
  lmfcn fcnheston = boost::bind(objheston, this, _1, _2, _3, _4, _5);
  lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
        nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4,
        fcnheston);
  //*   info is an integer output variable. if the user has terminated
  // execution, info is set to the (negative)
  //*     value of iflag. see description of fcn. otherwise, info is set as
  // follows.
  //	*     info = 0  improper input parameters.
  //	*     info = 1  both actual and predicted relative reductions in the sum
  // of squares are at most ftol.
  //	*     info = 2  relative error between two consecutive iterates is at
  // most xtol.
  //	*     info = 3  conditions for info = 1 and info = 2 both hold.
  //	*     info = 4  the cosine of the angle between fvec and any column of
  // the jacobian is at most gtol in absolute value.
  //	*     info = 5  number of calls to fcn has reached or exceeded maxfev.
  //	*     info = 6  ftol is too small. no further reduction in the sum of
  // squares is possible.
  //	*     info = 7  xtol is too small. no further improvement in the
  // approximate solution x is possible.
  //	*     info = 8  gtol is too small. fvec is orthogonal to the columns of
  // the jacobian to machine precision.
  QL_ENSURE(info != 4, "Heston Model Calibration Fails " << info);
  // the below is output result
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterTheta(x[2]);
  setParameterXi(x[3]);
  setParameterRho(x[4]);

  delete[] x;
  delete[] fvec;
  delete[] diag;
  delete[] fjac;
  delete[] ipvt;
  delete[] qtf;
  delete[] wa1;
  delete[] wa2;
  delete[] wa3;
  delete[] wa4;
};
void sVol::IVcalibrator(vector<double> maturitys, vector<double> forwards,
                        vector<double> strikes,
                        std::vector<double> marketQuotes, string quoteType) {
  int m = strikes.size(); // no. of observations
  weights_.resize(m);
  for (int i = 0; i < m; i++) {
    weights_[i] = exp(-abs(spot_ - strikes[i]) / spot_ * 100 * maturitys[i]);
  }
  double sum = 0;
  for (int i = 0; i < m; i++)
    sum += weights_[i];
  for (int i = 0; i < m; i++)
    weights_[i] /= sum;
  maturitys_.resize(m);
  maturitys_ = maturitys;
  forwards_.resize(m);
  forwards_ = forwards;
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  if (quoteType != "impliedVol")
    calibrator(maturitys, forwards, strikes, marketQuotes, quoteType);
  else {
    int n = 5;                 // no. of HESTON model paremeters
    double *x = new double[n]; // initial estimate of parameters vector
    x[0] = getParameterVar0();
    x[1] = getParameterKappa();
    x[2] = getParameterTheta();
    x[3] = getParameterXi();
    x[4] = getParameterRho();

    double *fvec = new double[m]; // no need to populate
    double ftol = 1e-08;          // tolerance
    double xtol = 1e-08;          // tolerance
    double gtol = 1e-08;          // tolerance
    int maxfev = 400;             // maximum function evaluations
    double epsfcn = 1e-08;        // tolerance
    double *diag = new double[n]; // some internal thing
    int mode = 1;                 // some internal thing
    double factor = 1;            // a default recommended value
    int nprint = 0;               // don't know what it does
    int info = 0;                 // output variable
    int nfev = 0; // output variable will store no. of function evals
    double *fjac = new double[m * n]; // output array of jacobian
    int ldfjac = m;                   // recommended setting
    int *ipvt = new int[n];           // for internal use
    double *qtf = new double[n];      // for internal use
    double *wa1 = new double[n];      // for internal use
    double *wa2 = new double[n];      // for internal use
    double *wa3 = new double[n];      // for internal use
    double *wa4 = new double[m];      // for internal use

    boost::function<void(sVol *, int, int, double *, double *, int *)>
        objheston = &sVol::objFcnIV;
    lmfcn fcnheston = boost::bind(objheston, this, _1, _2, _3, _4, _5);
    lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
          nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4,
          fcnheston);
    //*   info is an integer output variable. if the user has terminated
    // execution, info is set to the (negative)
    //*     value of iflag. see description of fcn. otherwise, info is set as
    // follows.
    //	*     info = 0  improper input parameters.
    //	*     info = 1  both actual and predicted relative reductions in the sum
    // of squares are at most ftol.
    //	*     info = 2  relative error between two consecutive iterates is at
    // most xtol.
    //	*     info = 3  conditions for info = 1 and info = 2 both hold.
    //	*     info = 4  the cosine of the angle between fvec and any column of
    // the jacobian is at most gtol in absolute value.
    //	*     info = 5  number of calls to fcn has reached or exceeded maxfev.
    //	*     info = 6  ftol is too small. no further reduction in the sum of
    // squares is possible.
    //	*     info = 7  xtol is too small. no further improvement in the
    // approximate solution x is possible.
    //	*     info = 8  gtol is too small. fvec is orthogonal to the columns of
    // the jacobian to machine precision.
    QL_ENSURE(info != 4, "Heston Model Calibration Fails " << info);
    // the below is output result
    setParameterVar0(x[0]);
    setParameterKappa(x[1]);
    setParameterTheta(x[2]);
    setParameterXi(x[3]);
    setParameterRho(x[4]);

    delete[] x;
    delete[] fvec;
    delete[] diag;
    delete[] fjac;
    delete[] ipvt;
    delete[] qtf;
    delete[] wa1;
    delete[] wa2;
    delete[] wa3;
    delete[] wa4;
  }
};

void sVol::FXcalibrator(vector<double> maturitys, vector<double> forwards,
                        vector<double> strikes,
                        std::vector<double> marketQuotes,
                        string /*quoteType*/) {
  int m = strikes.size(); // no. of observations
  // weights_.resize(m);
  // for (int i=0;i<m;i++)
  //	{
  //		weights_[i]=exp(-abs(spot_-strikes[i])/spot_*100*maturitys[i]);
  //	}
  // double sum = 0;
  // for (int i=0;i<m;i++) sum+=weights_[i];
  // for (int i=0;i<m;i++) weights_[i]/=sum;
  maturitys_.resize(m);
  maturitys_ = maturitys;
  forwards_.resize(m);
  forwards_ = forwards;
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  for (int i = 0; i < m; i++) {
    double vol = marketQuotes[i] * sqrt(maturitys[i]);
    double d1 = std::log(forwards[i] / strikes[i]) / vol + 0.5 * vol;
    double d2 = d1 - vol;
    marketQuotes_[i] =
        forwards[i] * cdf_normal(d1) - strikes[i] * cdf_normal(d2);
  }
  int n = 4;                 // no. of HESTON FX model paremeters
  double *x = new double[n]; // initial estimate of parameters vector
  x[0] = getParameterVar0();
  x[1] = getParameterTheta();
  x[2] = getParameterXi();
  x[3] = getParameterRho();

  double *fvec = new double[m]; // no need to populate
  double ftol = 1e-08;          // tolerance
  double xtol = 1e-08;          // tolerance
  double gtol = 1e-08;          // tolerance
  int maxfev = 400;             // maximum function evaluations
  double epsfcn = 1e-08;        // tolerance
  double *diag = new double[n]; // some internal thing
  int mode = 1;                 // some internal thing
  double factor = 1;            // a default recommended value
  int nprint = 0;               // don't know what it does
  int info = 0;                 // output variable
  int nfev = 0; // output variable will store no. of function evals
  double *fjac = new double[m * n]; // output array of jacobian
  int ldfjac = m;                   // recommended setting
  int *ipvt = new int[n];           // for internal use
  double *qtf = new double[n];      // for internal use
  double *wa1 = new double[n];      // for internal use
  double *wa2 = new double[n];      // for internal use
  double *wa3 = new double[n];      // for internal use
  double *wa4 = new double[m];      // for internal use

  boost::function<void(sVol *, int, int, double *, double *, int *)> objheston =
      &sVol::objFcnFX;
  lmfcn fcnheston = boost::bind(objheston, this, _1, _2, _3, _4, _5);
  lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
        nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4,
        fcnheston);
  //*   info is an integer output variable. if the user has terminated
  // execution, info is set to the (negative)
  //*     value of iflag. see description of fcn. otherwise, info is set as
  // follows.
  //	*     info = 0  improper input parameters.
  //	*     info = 1  both actual and predicted relative reductions in the sum
  // of squares are at most ftol.
  //	*     info = 2  relative error between two consecutive iterates is at
  // most xtol.
  //	*     info = 3  conditions for info = 1 and info = 2 both hold.
  //	*     info = 4  the cosine of the angle between fvec and any column of
  // the jacobian is at most gtol in absolute value.
  //	*     info = 5  number of calls to fcn has reached or exceeded maxfev.
  //	*     info = 6  ftol is too small. no further reduction in the sum of
  // squares is possible.
  //	*     info = 7  xtol is too small. no further improvement in the
  // approximate solution x is possible.
  //	*     info = 8  gtol is too small. fvec is orthogonal to the columns of
  // the jacobian to machine precision.
  QL_ENSURE(info != 4, "Heston Model FX Calibration Fails " << info);
  // the below is output result
  setParameterVar0(x[0]);
  setParameterTheta(x[1]);
  setParameterXi(x[2]);
  setParameterRho(x[3]);
  delete[] x;
  delete[] fvec;
  delete[] diag;
  delete[] fjac;
  delete[] ipvt;
  delete[] qtf;
  delete[] wa1;
  delete[] wa2;
  delete[] wa3;
  delete[] wa4;
};

// void sVol::calibrator(vector<double> maturitys, vector<double> forwards,
// vector<double> strikes, 	std::vector<double> marketQuotes, string
// quoteType)
//	{
//	int m = strikes.size();		//no. of observations
//
//	maturitys_.resize(m);
//	maturitys_ = maturitys;
//	forwards_.resize(m);
//	forwards_ = forwards;
//	strikes_.resize(m);
//	strikes_ = strikes;
//	marketQuotes_.resize(m);
//	marketQuotes_ = marketQuotes;
//	if(quoteType == "impliedVol")
//		{
//		for(int i=0; i<m; i++)
//			{
//			double vol = marketQuotes[i] * sqrt(maturitys[i]);
//			double d1 = std::log(forwards[i]/strikes[i])/vol
//+0.5*vol; 			double d2 = d1-vol;
//marketQuotes_[i] = forwards[i]*cdf_normal(d1) - strikes[i]*cdf_normal(d2);
//			}
//		}
//
//	int n = 5;					//no. of HESTON model
// paremeters 	double* x = new double[n];	//initial estimate of parameters
// vector 	x[0] = getParameterVar0(); 	x[1] = getParameterKappa(); 	x[2] =
//getParameterTheta(); 	x[3] = getParameterXi(); 	x[4] = getParameterRho();
//
//	double* fvec = new double[m]; //no need to populate
//	double ftol = 1e-08; //tolerance
//	double xtol = 1e-08; //tolerance
//	double gtol = 1e-08; //tolerance
//	int maxfev = 400; //maximum function evaluations
//	double epsfcn = 1e-08; //tolerance
//	double* diag = new double[n]; //some internal thing
//	int mode = 1; //some internal thing
//	double factor = 1; // a default recommended value
//	int nprint = 0; //don't know what it does
//	int info = 0; //output variable
//	int nfev = 0; //output variable will store no. of function evals
//	double* fjac = new double[m*n]; //output array of jacobian
//	int ldfjac = m; //recommended setting
//	int* ipvt = new int[n]; //for internal use
//	double* qtf = new double[n]; //for internal use
//	double* wa1 = new double[n]; //for internal use
//	double* wa2 = new double[n]; //for internal use
//	double* wa3 = new double[n]; //for internal use
//	double* wa4 = new double[m]; //for internal use
//
//	boost::function<void (sVol*, int,int,double*,double*,int*)> objheston =
//&sVol::objFcn; 	lmfcn fcnheston = boost::bind(objheston, this,
//_1,_2,_3,_4,_5); 	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag,
//mode, factor, 		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf,
// wa1, wa2, wa3, wa4, fcnheston);
//	//*   info is an integer output variable. if the user has terminated
// execution, info is set to the (negative)
//	//*     value of iflag. see description of fcn. otherwise, info is set
// as follows.
//	//	*     info = 0  improper input parameters.
//	//	*     info = 1  both actual and predicted relative reductions in
// the sum of squares are at most ftol.
//	//	*     info = 2  relative error between two consecutive iterates
// is at most xtol.
//	//	*     info = 3  conditions for info = 1 and info = 2 both hold.
//	//	*     info = 4  the cosine of the angle between fvec and any
// column of the jacobian is at most gtol in absolute value.
//	//	*     info = 5  number of calls to fcn has reached or exceeded
// maxfev.
//	//	*     info = 6  ftol is too small. no further reduction in the
// sum of squares is possible.
//	//	*     info = 7  xtol is too small. no further improvement in the
// approximate solution x is possible.
//	//	*     info = 8  gtol is too small. fvec is orthogonal to the
// columns of the jacobian to machine precision. 	QL_ENSURE(info != 4,
// "Heston Model Calibration Fails " << info);
//	//the below is output result
//	setParameterVar0(x[0]);
//	setParameterKappa(x[1]);
//	setParameterTheta(x[2]);
//	setParameterXi(x[3]);
//	setParameterRho(x[4]);
//
//	delete[] x;
//	delete[] fvec;
//	delete[] diag;
//	delete[] fjac;
//	delete[] ipvt;
//	delete[] qtf;
//	delete[] wa1;
//	delete[] wa2;
//	delete[] wa3;
//	delete[] wa4;
//	};

// void sVol::objFcn(int m, int n, double* x, double* fvec, int* iflag)
//	{
//	setParameterVar0(x[0]);
//	setParameterKappa(x[1]);
//	setParameterTheta(x[2]);
//	setParameterXi(x[3]);
//	setParameterRho(x[4]);
//	for(int i=0; i<m; i++)
//		fvec[i] = hestonPrice(maturitys_[i], forwards_[i], strikes_[i],
//"call") - marketQuotes_[i];
//	};

void sVol::objFcn(int m, int /*n*/, double *x, double *fvec, int * /*iflag*/) {
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterTheta(x[2]);
  setParameterXi(x[3]);
  setParameterRho(x[4]);
  for (int i = 0; i < m; i++)
    //	fvec[i] = (hestonPrice(maturitys_[i], forwards_[i], strikes_[i], "call")
    //- marketQuotes_[i])*weights_[i];
    fvec[i] = (hestonPrice(maturitys_[i], forwards_[i], strikes_[i], "call") -
               marketQuotes_[i]);
  // fvec[i] = (hestonPriceCF(maturitys_[i], forwards_[i], strikes_[i], "call")
  // - marketQuotes_[i]);
};
void sVol::objFcnIV(int m, int /*n*/, double *x, double *fvec,
                    int * /*iflag*/) {
  setParameterVar0(x[0]);
  setParameterKappa(x[1]);
  setParameterTheta(x[2]);
  setParameterXi(x[3]);
  setParameterRho(x[4]);
  double A = 0.0;
  if (x[1] * x[2] * 2.0 < x[3] * x[3])
    A = 1000000;

  for (int i = 0; i < m; i++)
    fvec[i] = (implied_vol(maturitys_[i], forwards_[i], strikes_[i],
                           hestonPrice(maturitys_[i], forwards_[i], strikes_[i],
                                       "call")) -
               marketQuotes_[i]) +
              A;
};
void sVol::objFcnFX(int m, int /*n*/, double *x, double *fvec,
                    int * /*iflag*/) {
  setParameterVar0(x[0]);
  setParameterTheta(x[1]);
  setParameterXi(x[2]);
  setParameterRho(x[3]);
  double A = 0.0;
  if (kappa_ * x[1] * 2.0 < x[2] * x[2])
    A = 10000;
  for (int i = 0; i < m; i++)
    fvec[i] =
        1.0 -
        (hestonPrice(maturitys_[i], forwards_[i], strikes_[i], "call") + A) /
            marketQuotes_[i];
};

// NOT USED!!!
double sVol::hestonIntegrandNET(complex<double> k, double maturity,
                                double forward, double strike) const {
  complex<double> i(0.0, 1.0);
  double b = 2.0 / (xi_ * xi_) * (kappa_ + rho_ * xi_);
  complex<double> d = sqrt(b * b + 4.0 * (k * k - i * k) / (xi_ * xi_));
  complex<double> r1 = (b + d) / 2.0;
  complex<double> g = (b + d) / (b - d);
  complex<double> Y = exp(0.5 * xi_ * xi_ * maturity * d);
  complex<double> F1 =
      kappa_ * theta_ *
      (r1 * maturity + 2.0 / (xi_ * xi_) * log((1.0 - g) / (1.0 - g * Y)));
  complex<double> F2 = r1 * (1.0 - Y) / (1.0 - g * Y);
  complex<double> integrand =
      exp(-i * k * log(forward / strike) + F1 + F2 * var0_) / (k * k - i * k);
  if (real(integrand) != real(integrand))
    return 0.0;
  else
    return real(integrand);
}
double sVol::hestonPriceNET(double maturity, double forward, double strike,
                            string optType) const {
  const double ki = 0.5;
  int kmax = std::max(1000.0, std::ceil(10.0 / std::sqrt(var0_ * maturity)));
  vector<double> int_x(kmax * 5), int_y(kmax * 5);
  int count = 0;
  for (double phi = 0.000001; phi < kmax; phi += 0.2) {
    int_x[count] = phi;
    complex<double> pass_phi(phi, ki);
    int_y[count] = hestonIntegrandNET(pass_phi, maturity, forward, strike);
    count += 1;
  }
  // computing the price
  double callPrice = forward - strike * trapz(int_x, int_y) / PI;
  if (optType == "Call" || optType == "call")
    return callPrice;
  else
    return callPrice + strike - forward;
}
vector<double> sVol::simulationHestonMax(vector<double> times,
                                         vector<double> forwards) const {
  int N = times.size();
  vector<double> pathF(N);
  int T = 0;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  double Smax = S;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * DT + xi_ * sqrt(abs(V)) * SDT * z2 +
         0.25 * xi_ * xi_ * DT * (z2 * z2 - 1); // Milstein scheme
    Smax = (Smax < S) ? S : Smax;
    // V += kappa_*(theta_-V)*DT +xi_*sqrt(abs(V))*SDT*z2;  //Euler scheme
    if (i * DT <= times[T] && (i + 1) * DT > times[T]) {
      pathF[T] = forwards[T] * Smax / spot_;
      T++;
    }
  }
  return pathF;
}

vector<double> sVol::simulationHestonCliq(vector<double> times,
                                          vector<double> forwards, double gcap,
                                          double gfloor, double lcap,
                                          double lfloor, double alpha) const {
  int N = times.size();
  vector<double> pathF(N);
  int T = 0;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  double S1, S0 = S;
  double cRet = 0.0;
  double R = 0.0;
  int Day = 1;
  double dt = 1.0 / 365;
  double accT = 0.0;

  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * SDT *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * DT + xi_ * sqrt(abs(V)) * SDT * z2 +
         0.25 * xi_ * xi_ * DT * (z2 * z2 - 1); // Milstein scheme
    accT += DT;
    if (accT >= Day * dt) {
      Day++;
      S1 = interpolateF(times, forwards, Day * dt) * S / spot_;
      R = S1 / (alpha * S0) - 1.0;
      cRet += std::max(lfloor, std::min(R, lcap));
      S0 = S1;
    }
    if (i * DT <= times[T] && (i + 1) * DT > times[T]) {
      pathF[T] = std::min(std::max(cRet, gfloor), gcap);
      T++;
    }
  }
  return pathF;
}

vector<double> sVol::simulationHestonDNT(vector<double> times,
                                         vector<double> forwards, double UP,
                                         double DOWN) const {
  int N = times.size();
  vector<double> pathF(N, 0.0);
  int T = 0;
  int pr = 1e6;
  int iUP = UP * pr;
  int iDOWN = DOWN * pr;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  int S1;
  double dt = 0.0009;
  double sdt = sqrt(dt);

  // for(int i=0; i<=int(times[N-1]/DT); i++)
  //{
  //	double z1 = random_normal();
  //	double z2 = rho_*z1 + rho2*random_normal();
  //	S += S*sqrt(abs(V))*SDT*z1;  // This line (S) must put BEFORE line
  //(V)!!! 	V += kappa_*(theta_-V)*DT +xi_*sqrt(abs(V))*SDT*z2
  //+0.25*xi_*xi_*DT*(z2*z2-1);  //Milstein scheme 	S1 =
  // interpolateF(times,forwards,i*DT)*S/spot_*pr; 	if
  // (!((S1<=iUP)&&(S1>=iDOWN))) break; 	if(i*DT <= times[T] && (i+1)*DT
  // > times[T])
  //	{
  //		pathF[T] = 1.0;
  //		T++;
  //	}
  // }
  for (int i = 0; i <= int(times[N - 1] / dt); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * sdt *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * dt + xi_ * sqrt(abs(V)) * sdt * z2 +
         0.25 * xi_ * xi_ * dt * (z2 * z2 - 1); // Milstein scheme
    S1 = interpolateF(times, forwards, i * dt) * S / spot_ * pr;
    if (!((S1 <= iUP) && (S1 >= iDOWN)))
      break;
    if (i * dt <= times[T] && (i + 1) * dt > times[T]) {
      pathF[T] = 1.0;
      T++;
    }
  }
  return pathF;
}

double sVol::simulationHestonDNTdt(vector<double> times,
                                   vector<double> forwards, double UP,
                                   double DOWN, double dt) const {
  int N = times.size();

  int Nstep = (times[N - 1] / dt);
  int pr = 1e6;
  int iUP = UP * pr;
  int iDOWN = DOWN * pr;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  int S1;
  double sdt = sqrt(dt);
  for (int i = 0; i <= Nstep; i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * sdt *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * dt + xi_ * sqrt(abs(V)) * sdt * z2 +
         0.25 * xi_ * xi_ * dt * (z2 * z2 - 1); // Milstein scheme
    S1 = interpolateF(times, forwards, i * dt) * S / spot_ * pr;
    if (!((S1 <= iUP) && (S1 >= iDOWN)))
      return 0.0;
  }
  return 1.0;
}
int sVol::simulationHestonDNTdtS(vector<double> times, vector<double> forwards,
                                 double maturity, double UP, double DOWN,
                                 double dt) const {
  int Nstep = maturity / dt;
  int pr = 1e6;
  int iUP = UP * pr;
  int iDOWN = DOWN * pr;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  int S1;
  double sdt = sqrt(dt);
  for (int i = 0; i <= Nstep; i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * sdt *
         z1; // This line (S) must put BEFORE line (V)!!!
    V += kappa_ * (theta_ - V) * dt + xi_ * sqrt(abs(V)) * sdt * z2 +
         0.25 * xi_ * xi_ * dt * (z2 * z2 - 1); // Milstein scheme
    S1 = interpolateF(times, forwards, i * dt) * S / spot_ * pr;
    if (!((S1 <= iUP) && (S1 >= iDOWN)))
      return 0;
  }
  return 1;
}
int sVol::simulationHestonDNTdtE(vector<double> times, vector<double> forwards,
                                 double maturity, double UP, double DOWN,
                                 double dt) const {
  int Nstep = maturity / dt;
  int pr = 1e6;
  int iUP = UP * pr;
  int iDOWN = DOWN * pr;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  int S1;
  double sdt = sqrt(dt);
  for (int i = 0; i <= Nstep; i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * sdt *
         z1; // This line (S) must put BEFORE line (V)!!!
    // V += kappa_*(theta_-V)*dt +xi_*sqrt(abs(V))*sdt*z2
    // +0.25*xi_*xi_*dt*(z2*z2-1);  //Milstein scheme
    V += kappa_ * (theta_ - V) * dt +
         xi_ * sqrt(abs(V)) * sdt * z2; // Euler scheme
    S1 = interpolateF(times, forwards, i * dt) * S / spot_ * pr;
    if (!((S1 <= iUP) && (S1 >= iDOWN)))
      return 0;
  }
  return 1;
}

vector<double> sVol::simulationHestonUpnOut(vector<double> times,
                                            vector<double> forwards,
                                            double BARRIER) const {

  int N = times.size();
  vector<double> pathF(N, 0.0);
  if (BARRIER <= spot_)
    return pathF;
  int T = 0;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  double S1;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * SDT * z1;
    V += kappa_ * (theta_ - V) * DT + xi_ * sqrt(abs(V)) * SDT * z2 +
         0.25 * xi_ * xi_ * DT * (z2 * z2 - 1); // Milstein scheme
    S1 = interpolateF(times, forwards, i * DT) * S / spot_;
    if (S1 > BARRIER)
      break;
    if (i * DT <= times[T] && (i + 1) * DT > times[T]) {
      pathF[T] = 1.0;
      T++;
    }
  }
  return pathF;
}

vector<double> sVol::simulationHestonDownnOut(vector<double> times,
                                              vector<double> forwards,
                                              double BARRIER) const {

  int N = times.size();
  vector<double> pathF(N, 0.0);
  if (BARRIER >= spot_)
    return pathF;
  int T = 0;
  double V = var0_, S = spot_, rho2 = sqrt(1 - rho_ * rho_);
  double S1;
  for (int i = 0; i <= int(times[N - 1] / DT); i++) {
    double z1 = random_normal();
    double z2 = rho_ * z1 + rho2 * random_normal();
    S += S * sqrt(abs(V)) * SDT * z1;
    V += kappa_ * (theta_ - V) * DT + xi_ * sqrt(abs(V)) * SDT * z2 +
         0.25 * xi_ * xi_ * DT * (z2 * z2 - 1); // Milstein scheme
    S1 = interpolateF(times, forwards, i * DT) * S / spot_;
    if (S1 < BARRIER)
      break;
    if (i * DT <= times[T] && (i + 1) * DT > times[T]) {
      pathF[T] = 1.0;
      T++;
    }
  }
  return pathF;
}
} // namespace velesquant
