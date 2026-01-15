
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <omp.h>
#include <ql/errors.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <velesquant/volatility/lm.h>
#include <velesquant/pde_solvers/cyclic_reduction.h>
#include <velesquant/pde_solvers/short_rate_1f_pde.h>
using namespace std;
#pragma warning(disable : 4715)
#pragma warning(disable : 4018)
#define MR 256

namespace velesquant {

static const int DAYS_A_YEAR = 4001;
static const double delT = 10.0 / (DAYS_A_YEAR - 1.0);

void ShortRate1FPDE::buildGrid(double Rmax, double factor) {
  gridR_ = grMesh(MR - 1, Rmax, factor);
  a_.resize(MR);
  b_.resize(MR);
  c_.resize(MR);
  e_.resize(MR);
  f_.resize(MR);
  g_.resize(MR);

#pragma omp parallel for
  for (int i = 0; i < MR; i++) {

    double R = gridR_[i];
    double Ru = (i == MR - 1) ? vright_ : gridR_[i + 1];
    double Rl = (i == 0) ? vleft_ : gridR_[i - 1];
    a_[i] = 2.0 / (R - Rl) / (Ru - Rl);
    b_[i] = -2.0 / (R - Rl) / (Ru - R);
    c_[i] = 2.0 / (Ru - R) / (Ru - Rl);
    e_[i] = -(Ru - R) / (R - Rl) / (Ru - Rl);
    f_[i] = (Ru - 2 * R + Rl) / (R - Rl) / (Ru - R);
    g_[i] = (R - Rl) / (Ru - Rl) / (Ru - R);
  }
  for (int r = 0; r < MR; r++) {
    if (gridR_[r] >= R0_) {
      iR0_ = r;
      break;
    }
  }
};

void ShortRate1FPDE::discountBack(double t0, double Tn, vector<double> &f) {
  int Nt = int(Tn / delT + 1);
  int N0 = int(t0 / delT + 1);
  for (int t = Nt - 2; t >= N0; t--)
    oneStepBackward(t, f);
};

double ShortRate1FPDE::pricingZB(double Maturity) {
  vector<double> payoff(MR, 1.0);
  discountBack(0, Maturity, payoff);
  return payoff[iR0_];
};

void ShortRate1FPDE::pricingCouponBondt(double t0, double Tn, double Coupon,
                                        double PayFrequency,
                                        vector<double> &f) {
  int Nt = int(Tn / delT + 1); // working daily time grid
  int iExpiry = int(t0 / delT + 1);
  int Ncoupon = int((Tn - t0) / PayFrequency + 0.5);
  for (int t = Nt - 2; t >= iExpiry; t--) // swap
  {
    double couponTime = t0 + Ncoupon * PayFrequency;
    if (couponTime > t * delT && couponTime <= (t + 1) * delT) {
      for (int r = 0; r < MR; r++)
        f[r] += PayFrequency * Coupon;
      Ncoupon--;
    }
    oneStepBackward(t, f); // from t+1 to t
  }
};

double ShortRate1FPDE::pricingCouponBond(double Expiry, double Tenor,
                                         double Coupon, double PayFrequency) {
  vector<double> payoff(MR, 1);
  pricingCouponBondt(Expiry, Expiry + Tenor, Coupon, PayFrequency, payoff);
  if (Expiry > 0)
    discountBack(0, Expiry, payoff);
  return payoff[iR0_];
};

double ShortRate1FPDE::pricingCBO(double Expiry, double Tenor, double Coupon,
                                  double Strike, double PayFrequency,
                                  const std::string &type) {
  vector<double> f(MR, 1);
  int Nt = int(Expiry / delT + 1);
  pricingCouponBondt(Expiry, Expiry + Tenor, Coupon, PayFrequency,
                     f); // Price CB
  for (int i = 0; i < MR; i++)
    f[i] = (type == "Call") ? max(0.0, f[i] - Strike)
                            : max(0.0, Strike - f[i]); // Apply Payoff
  if (Expiry > 0)
    discountBack(0, Expiry, f); // Discount to time t=0
  return f[iR0_];               // Find the correct value
};

double ShortRate1FPDE::pricingZBO(double Expiry, double Maturity, double Strike,
                                  const std::string &type) {
  int Nt = int(Expiry / delT + 1); // working daily time grid?
  vector<double> f(MR, 1);
  // Price ZB from T to t
  discountBack(Expiry, Expiry + Maturity, f);
  // Apply Payoff to f: max(0,f-K) || max(0,K-f)
  for (int i = 0; i < MR; i++)
    f[i] =
        (type == "Call") ? max(0.0, f[i] - Strike) : max(0.0, -f[i] + Strike);
  if (Expiry > 0)
    discountBack(0, Expiry, f); // Discount to time t=0
  return f[iR0_];
};

double ShortRate1FPDE::pricingSwap(double Expiry, double Tenor, double Strike,
                                   double PayFrequency) {
  vector<double> payoff(MR, -1);
  pricingCouponBondt(Expiry, Expiry + Tenor, -Strike, PayFrequency, payoff);
#pragma omp parallel for
  for (int r = 0; r < MR; r++)
    payoff[r] += 1.0; // swap payoff
  if (Expiry > 0)
    discountBack(0, Expiry, payoff);
  return payoff[iR0_];
};

double ShortRate1FPDE::pricingCallableSwap(double Expiry, double Tenor,
                                           vector<double> Exercises,
                                           double Coupon, double Strike,
                                           double PayFrequency,
                                           const std::string &type) {
  vector<double> payoff(MR, -1), call_value(MR);
  Exercises.insert(Exercises.begin(), Expiry);
  Exercises.push_back(Expiry + Tenor);
  int Ne = Exercises.size();
  for (int e = Ne - 2; e >= 0; e--) {
    pricingCouponBondt(Exercises[e], Exercises[e + 1], -Coupon, PayFrequency,
                       payoff);
    if (e == Ne - 2) {
#pragma omp parallel for
      for (int r = 0; r < MR; r++)
        call_value[r] = (type == "Call")
                            ? max(0.0, (payoff[r] + 1.0) - Strike)
                            : max(0.0, -(payoff[r] + 1.0) + Strike); // payoff
    } else {
      discountBack(Exercises[e], Exercises[e + 1], call_value);
#pragma omp parallel for
      for (int r = 0; r < MR; r++)
        call_value[r] = (type == "Call")
                            ? max(0.0, (payoff[r] + 1.0) - Strike)
                            : max(0.0, -(payoff[r] + 1.0) + Strike); // payoff
    }
  }
  if (Expiry > 0)
    discountBack(0, Expiry, call_value);
  return call_value[iR0_];
};

double ShortRate1FPDE::pricingSwaption(double Expiry, double Tenor,
                                       double Strike, double PayFrequency) {
  vector<double> payoff(MR, -1.0);
  pricingCouponBondt(Expiry, Expiry + Tenor, -Strike, PayFrequency, payoff);
  for (int r = 0; r < MR; r++)
    payoff[r] = max(0.0, payoff[r] + 1.0); // swaption payoff
  if (Expiry > 0)
    discountBack(0, Expiry, payoff);
  return payoff[iR0_];
};

double ShortRate1FPDE::pricingBermudan(double Expiry, double Tenor,
                                       vector<double> Exercises, double Strike,
                                       double PayFrequency) {
  std::vector<double> payoff(MR, -1.0), swaption(MR);
  Exercises.insert(Exercises.begin(), Expiry);
  Exercises.push_back(Expiry + Tenor);
  int Ne = Exercises.size();
  for (int e = Ne - 2; e >= 0; e--) {
    pricingCouponBondt(Exercises[e], Exercises[e + 1], -Strike, PayFrequency,
                       payoff);
    if (e == Ne - 2) {
#pragma omp parallel for
      for (int r = 0; r < MR; r++)
        swaption[r] = max(0.0, payoff[r] + 1.0); // swaption payoff
    } else {
      discountBack(Exercises[e], Exercises[e + 1], swaption);
#pragma omp parallel for
      for (int r = 0; r < MR; r++)
        swaption[r] =
            max(swaption[r], payoff[r] + 1.0); // bermudan swaption payoff
    }
  }
  if (Expiry > 0)
    discountBack(0, Expiry, swaption);
  return swaption[iR0_];
};

vector<double> ShortRate1FPDE::calculateDFs(vector<double> &timePoints) {
  int n = timePoints.size();
  vector<int> iTs(n);
  for (int i = 0; i < n; i++)
    iTs[i] = int(timePoints[i] / delT + 1);
  // buildGrid();
  vector<double> inV(MR, 0.0), outV(MR), DFs(n, 0.0);
  for (int r = 1; r < MR - 1; r++)
    if (gridR_[r] >= R0_) {
      inV[r] = 2.0 / (gridR_[r + 1] -
                      gridR_[r - 1]); // Dirac delta as initial distribution
      break;
    }
  for (int i = 0; i < n; i++) {
    int sT = 0;
    if (i > 0)
      sT = iTs[i - 1] - 1;
    for (int t = sT; t < iTs[i] - 1; t++)
      oneStepForward(t, inV);
    DFs[i] = trapezoidal(inV);
    QL_ENSURE(DFs[i] <= 1.001, "pricingDF Functor Fails " << DFs[i]);
  }
  return DFs;
};

void ShortRate1FPDE::calibrator(vector<double> timeDFs, vector<double> DFs,
                                vector<defSwap> swapQuotes) {
  QL_REQUIRE(!swapQuotes.empty(), "Swap quotes cannot be empty");
  QL_REQUIRE(!timeDFs.empty(), "Time DFs cannot be empty");
  DFs_ = DFs;
  timeDFs_ = timeDFs;
  quoteSwap_ = swapQuotes;

  timeThetas_.resize(1);
  timeThetas_[0] = 0.5;
  double dt = timeDFs_[timeDFs_.size() - 1] / 20;
  timeThetas_[0] = dt;
  for (double i = timeThetas_[0] + dt; i <= timeDFs_[timeDFs_.size() - 1];
       i += dt)
    timeThetas_.push_back(i);
  thetas_.resize(timeThetas_.size());
  for (int i = 0; i < timeThetas_.size(); i++)
    thetas_[i] = .001;
  DFinterp_.resize(timeThetas_.size());
  for (int i = 0; i < timeThetas_.size(); i++)
    DFinterp_[i] = getDFinterp(timeThetas_[i]);

  timeSigmas_.resize(0);
  timeSigmas_.push_back(quoteSwap_[0].Expiry);
  double pos = timeSigmas_[0];
  for (int i = 1; i < quoteSwap_.size(); i++) {
    if (quoteSwap_[i].Expiry > pos) {
      timeSigmas_.push_back(quoteSwap_[i].Expiry);
      pos = quoteSwap_[i].Expiry;
    }
  }
  sigmas_.resize(timeSigmas_.size(), 0.01);

  int ns = sigmas_.size();
  int n = ns + 1;            // no. of model volatility term paremeters
  double *x = new double[n]; // initial estimate of parameters vector
  for (int i = 0; i < ns; i++)
    x[i] = sigmas_[i]; // sigmas initial value
  calibrateKappa();

  x[n - 1] = alpha_ / 1000;  // beta initial value
  int m = quoteSwap_.size(); // no. of observations
  QL_ENSURE(m >= n, "too much freedom in Calibration  " << m - n);

  for (int i = 0; i < m; i++)
    quoteSwap_[i].Value = swaptionATM(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                                      quoteSwap_[i].VolATM);

  double *fvec = new double[m]; // no need to populate
  double ftol = 1e-6;           // tolerance
  double xtol = 1e-6;           // tolerance
  double gtol = 1e-6;           // tolerance
  int maxfev = 800;             // maximum function evaluations
  double epsfcn = 1e-7;         // tolerance
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
  boost::function<void(ShortRate1FPDE *, int, int, double *, double *, int *)>
      obj = &ShortRate1FPDE::objFcnCalibrator;

  lmfcn fcn = boost::bind(obj, this, _1, _2, _3, _4, _5);
  in_calibration_ = true;
  lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
        nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);
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
  in_calibration_ = false;
  QL_ENSURE(info != 4, "Model Calibration Fails " << info);
  // the below is output result
  for (int i = 0; i < ns; i++)
    sigmas_[i] = fabs(x[i]);                      // sigmas final value
  beta_ = (alpha_ - 1e-32) * sin(x[n - 1]) / 500; // beta final value

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

double ShortRate1FPDE::pen_fun(double *x, int n) {
  double penalty = 0.0;
  double lambda = 1e6;
  // for(int i=0; i<n-2;
  // i++)penalty+=((x[i]<=0)||(x[i]>=1))?lambda*pow(x[i],2):0;
  penalty += (beta_ >= alpha_) ? lambda * beta_ * beta_ : 0;
  return penalty;
}

void ShortRate1FPDE::objFcnCalibrator(int m, int n, double *x, double *fvec,
                                      int *iflag) {
  int ns = sigmas_.size();
  for (int i = 0; i < ns; i++)
    sigmas_[i] = fabs(x[i]);
  beta_ = (alpha_ - 1e-32) * sin(x[n - 1]) / 500;
  termStructureCalibrator();
  for (int i = 0; i < m; i++) {
    double model =
        pricingSwaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                        quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
    double diff = 1 - quoteSwap_[i].Value / model;
    fvec[i] = diff; //+ penalty;
  }
};

void ShortRate1FPDE::termStructureCalibrator() {
  int n = timeThetas_.size();
  vector<int> iTs(n);
  for (int i = 0; i < n; i++)
    iTs[i] = int(timeThetas_[i] / delT + 1);

  vector<double> inV(MR, 0.0), outV(MR), lastV(MR);
  for (int r = 1; r < MR - 1; r++)
    if (gridR_[r] >= R0_) {
      inV[r] = 2.0 / (gridR_[r + 1] -
                      gridR_[r - 1]); // Dirac delta as initial distribution
      break;
    }
#pragma omp parallel for
  for (int i = 0; i < MR; i++)
    lastV[i] = inV[i]; // remember the starting density
  for (int i = 0; i < n; i++) {
    int sT = 0;
    if (i > 0)
      sT = iTs[i - 1] - 1;
    double df = 0.0;
    int niter = 0;
    do {
      niter++;
#pragma omp parallel for
      for (int i = 0; i < MR; i++)
        inV[i] = lastV[i]; // reset to the starting density
      if (thetas_[i] == 0.0)
        thetas_[i] = 0.0005; // case of 0 value handler
      thetas_[i] *= 1.001;   // 0.1% up for theta parameter
      for (int t = sT; t < iTs[i] - 1; t++)
        oneStepForward(t, inV);
      double dfUP = trapezoidal(inV);
#pragma omp parallel for
      for (int i = 0; i < MR; i++)
        inV[i] = lastV[i]; // reset to the starting density
      thetas_[i] /= 1.001; // back to theta
      for (int t = sT; t < iTs[i] - 1; t++)
        oneStepForward(t, inV);
      df = trapezoidal(inV);
      thetas_[i] *=
          (1 + 0.001 * (DFinterp_[i] - df) / (dfUP - df)); // Newton iteration
      thetas_[i] =
          max(-0.0050, thetas_[i]); // Cutoff the nagative value (BETTER)
    } while (fabs(1.0 - df / DFinterp_[i]) >= 1.0E-6 && niter < 7);
#pragma omp parallel for
    for (int i = 0; i < MR; i++)
      inV[i] = lastV[i]; // reset to the starting density

    for (int t = sT; t < iTs[i] - 1; t++)
      oneStepForward(t, inV);
#pragma omp parallel for
    for (int i = 0; i < MR; i++)
      lastV[i] = inV[i]; // remember the starting density
  };
  return;
};

void ShortRate1FPDE::oneStepBackward(const int t, vector<double> &inV) {
  array<double, MR> l, c, u, d;
  double tm = (t + .5) * delT;
  double sigma = whichSigma(tm);
  double sigma2 = sigma * sigma;
  double theta = whichTheta(tm);
  double ndt = 1.0 / delT;
#pragma omp parallel for
  for (int r = 0; r < MR; r++) {
    double f1 = ((r == 0) || (r == MR - 1)) ? 0 : (theta - kappa_ * gridR_[r]);
    double f2 = sigma2 * pow(alpha_ + beta_ * gridR_[r], 2 * gamma_);
    double f3 = -gridR_[r];
    l[r] = 0.5 * (f1 * e_[r] + f2 / 2 * a_[r]);
    c[r] = -ndt + 0.5 * (f1 * f_[r] + f2 / 2 * b_[r] + f3);
    u[r] = 0.5 * (f1 * g_[r] + f2 / 2 * c_[r]);
    if (r == 0) {
      u[0] += l[0];
      d[0] = -((2 * ndt + c[r]) * inV[0] + u[r] * inV[1]);
    } else if (r == MR - 1) {
      l[MR - 1] += u[MR - 1];
      d[MR - 1] = -(l[MR - 1] * inV[MR - 2] + (2.0 * ndt + c[r]) * inV[MR - 1]);
    } else
      d[r] = -(l[r] * inV[r - 1] + (2.0 * ndt + c[r]) * inV[r] +
               u[r] * inV[r + 1]);
  }
  l[0] = 0.0;
  u[MR - 1] = 0.0;
  TriDiagonalSolve(MR, l, c, u, d, inV);
};

void ShortRate1FPDE::oneStepForward(const int T, vector<double> &inV) {
  array<double, MR> l, c, u, d;
  double Tm = (T + .5) * delT;
  double sigma = whichSigma(Tm);
  double theta = whichTheta(Tm);
  double sigma2 = sigma * sigma;
  double p1 = 2 * sigma2 * beta_ * gamma_;
  double p2 = 0.5 * p1 * beta_ * (2 * gamma_ - 1);
  double ndt = 1.0 / delT;

#pragma omp parallel for
  for (int r = 0; r < MR; r++) {
    double f1 = ((r == 0) || (r == MR - 1)) ? 0 : F1(r, p1, theta);
    double f2 = -F2(r, sigma2);
    double f3 = -F3(r, p2);
    l[r] = 0.5 * (f1 * e_[r] + f2 / 2 * a_[r]);
    c[r] = ndt + 0.5 * (f1 * f_[r] + f2 / 2 * b_[r] + f3);
    u[r] = 0.5 * (f1 * g_[r] + f2 / 2 * c_[r]);
    if (r == 0) {
      u[r] += l[r];
      d[r] = -((-2 * ndt + c[r]) * inV[r] + u[r] * inV[r + 1]);
    } else if (r == MR - 1) {
      l[r] += u[r];
      d[r] = -(l[r] * inV[r - 1] + (-2.0 * ndt + c[r]) * inV[r]);
    } else
      d[r] = -(l[r] * inV[r - 1] + (-2.0 * ndt + c[r]) * inV[r] +
               u[r] * inV[r + 1]);
  }
  l[0] = 0.0;
  u[MR - 1] = 0.0;

  // Solve PDE
  TriDiagonalSolve(MR, l, c, u, d, inV);
};

double ShortRate1FPDE::F1(int j, double p, double theta) {
  return theta - kappa_ * gridR_[j] -
         p * pow(alpha_ + beta_ * gridR_[j], 2 * gamma_ - 1);
}
double ShortRate1FPDE::F2(int j, double sigma2) {
  return sigma2 * pow(alpha_ + beta_ * gridR_[j], 2 * gamma_);
}

double ShortRate1FPDE::F3(int j, double p) {
  return p * pow(alpha_ + beta_ * gridR_[j], 2 * gamma_ - 2) + kappa_ -
         gridR_[j];
}

double ShortRate1FPDE::trapezoidal(vector<double> &inV) {
  double value = 0.0;
#pragma omp parallel for reduction(+ : value)
  for (int r = 1; r < MR; r++)
    value += 0.5 * (gridR_[r] - gridR_[r - 1]) * (inV[r] + inV[r - 1]);
  return value;
};

double ShortRate1FPDE::whichSigma(double t) const {
  int n = sigmas_.size();
  if (t < timeSigmas_[0])
    return sigmas_[0];
  if (t >= timeSigmas_[n - 1])
    return sigmas_[n - 1];
  for (int i = 1; i < n; i++)
    if ((t >= timeSigmas_[i - 1]) && (t < timeSigmas_[i]))
      return sigmas_[i];
};
double ShortRate1FPDE::whichTheta(double t) const {
  int n = thetas_.size();
  if (t < timeThetas_[0])
    return thetas_[0];
  if (t >= timeThetas_[n - 1])
    return thetas_[n - 1];
  for (int i = 1; i < n; i++)
    if (t < timeThetas_[i])
      return thetas_[i];
};

double ShortRate1FPDE::getDFinterp(double t) {
  QuantLib::LogLinearInterpolation interp(timeDFs_.begin(), timeDFs_.end(),
                                          DFs_.begin());
  double a = interp(t, false);
  return a;
};

vector<double> ShortRate1FPDE::grMesh(int Mv, double Rmax, double factor) {
  int N = Mv;
  vector<double> Mesh(Mv + 1);
  double d3 = Rmax / (N * factor);
  double start = asinh((-Rmax - R0_) / d3);
  double dx = 1.0 / N * (asinh((Rmax - R0_) / d3) - start);

  for (int i = 0; i <= N; i++)
    Mesh[i] = R0_ + d3 * sinh(start + i * dx);
  vleft_ = Mesh[0] - (Mesh[1] - Mesh[0]);
  vright_ = Mesh[Mv] + (Mesh[Mv] - Mesh[Mv - 1]);
  for (int i = 0; i < N; i++)
    if ((Mesh[i] < R0_) && (R0_ < Mesh[i + 1])) {
      Mesh[i + 1] = R0_;
      break;
    }
  return Mesh;
};

void ShortRate1FPDE::calibrateKappa() {
  int MAX_ITER = 1000;
  double X_TOL = 1e-5;
  double F_TOL = 1e-5;
  double a0 = .1;
  double a1 = a0;
  double dF, F;
  double min = 1e-6;
  double max = 1;
  int N = quoteSwap_.size();
  vector<double> iv_ratio;
  vector<double> bond_ratio;
  vector<int> index;
  for (int i = 0; i < N - 1; i++) {
    if (quoteSwap_[i].Expiry == quoteSwap_[i + 1].Expiry) {
      double ex = quoteSwap_[i].Expiry;
      double t1 = quoteSwap_[i].Tenor + ex;
      double t2 = quoteSwap_[i + 1].Tenor + ex;
      index.push_back(i);
      iv_ratio.push_back(quoteSwap_[i + 1].VolATM / quoteSwap_[i].VolATM);
      bond_ratio.push_back((getDFinterp(ex) - getDFinterp(t1)) /
                           (getDFinterp(ex) - getDFinterp(t2)));
    }
  }

  a0 = kappa_;
  double delta;
  for (int i = 0; i < MAX_ITER; i++) {
    F = objKappa(a0, iv_ratio, bond_ratio, index);
    double Fu = objKappa(a0 + X_TOL, iv_ratio, bond_ratio, index);
    double Fd = objKappa(a0 - X_TOL, iv_ratio, bond_ratio, index);
    dF = 0.5 * (Fu - Fd) / X_TOL;
    if ((fabs(F) < F_TOL) || (fabs(dF) < F_TOL))
      break;
    a1 = a0;
    delta = F / dF;
    a1 -= delta;
    if (a1 <= min) {
      delta = 0.5F * (a0 - min);
      a1 = a0 - delta;
      if ((a1 == min) || (a1 == max))
        break;
    } else if (a1 >= max) {
      delta = 0.5F * (a0 - max);
      a1 = a0 - delta;
      if ((a1 == min) || (a1 == max))
        break;
    }
    if (delta > 0)
      max = a0;
    else
      min = a0;
    if (abs(a1 - a0) < X_TOL * fabs(a1) || (a1 != a1))
      break;
    a0 = a1;
    if (fabs(max - min) < X_TOL)
      break;
  }
  QL_ENSURE(a0 == a0, "Failed to calibrate Kappa \n");
  kappa_ = a0;
};

double ShortRate1FPDE::objKappa(double a, const vector<double> &ivr,
                                const vector<double> &br,
                                const vector<int> &ind) {
  int N = ivr.size();
  double sum = 0.0;
  for (int i = 1; i < N; i++) {
    double ex = quoteSwap_[ind[i]].Expiry;
    double t1 = quoteSwap_[ind[i]].Tenor + ex;
    double t2 = quoteSwap_[ind[i] + 1].Tenor + ex;
    double vs = fabs(br[i] * Bratio(a, ex, t1, t2)) - ivr[i];
    sum += vs * vs;
  }
  return sum;
};

double ShortRate1FPDE::Bratio(double a, double Mi, double Tk, double Tj) {
  if (a == 0.0)
    return 100000;
  else {
    double t1 = (Tj - Mi);
    double t2 = (Tk - Mi);
    return (1.0 - exp(-a * t1)) / (1.0 - exp(-a * t2));
  }
};

vector<double> ShortRate1FPDE::SwaptionDiagnostic(double Expiry, double Tenor,
                                                  double Strike,
                                                  double PayFrequency) {
  int Nt = int((Expiry + Tenor) / delT + 1); // working daily time grid
  int iExpiry = int(Expiry * 2520 + 0.5) / 10;
  int Ncoupon = int(Tenor / PayFrequency + 0.5);
  vector<double> payoff(MR, -1.0), pv(MR);
  for (int t = Nt - 2; t >= iExpiry; t--) // swap
  {
    double couponTime = Expiry + Ncoupon * PayFrequency;
    if (couponTime > t * delT && couponTime <= (t + 1) * delT) {
#pragma omp parallel for
      for (int r = 0; r < MR; r++)
        payoff[r] -= PayFrequency * Strike;
      Ncoupon--;
    }
    oneStepBackward(t, payoff);
  }
#pragma omp parallel for

  for (int r = 0; r < MR; r++) // swaption payoff
    payoff[r] = max(0.0, payoff[r] + 1.0);
  for (int t = iExpiry - 1; t >= 0; t--)
    oneStepBackward(t, payoff);
  return payoff;
};

vector<double> ShortRate1FPDE::RiskNeutralDensity(double T1, double T2) {
  int iTs[] = {int(T1 / delT + 1), int(T2 / delT + 1)};
  vector<double> inV(MR, 0.0);
  for (int r = 1; r < MR - 1; r++)
    if (gridR_[r] >= R0_) {
      inV[r] = 2.0 / (gridR_[r + 1] -
                      gridR_[r - 1]); // Dirac delta as initial distribution
      break;
    }
  for (int t = iTs[0]; t < iTs[1] - 1; t++)
    oneStepForward(t, inV);
  return inV;
};

double ShortRate1FPDE::getSwapRate(double Expiry, double Tenor,
                                   double PayFrequency) {
  double SwapLow = pricingSwap(Expiry, Tenor, 0.01, PayFrequency);
  double SwapHigh = pricingSwap(Expiry, Tenor, 0.02, PayFrequency);
  double SwapRate = 0.01 + SwapLow / (SwapLow - SwapHigh) * 0.01;
  QL_ENSURE(SwapRate >= -0.0005,
            "getSwapRate Functor Fails " << SwapRate); // Negative 50BP Rate
  return SwapRate;
};

double ShortRate1FPDE::getImpVolATM(double Expiry, double Tenor,
                                    double PayFrequency) {
  double SwapLow = pricingSwap(Expiry, Tenor, 0.00, PayFrequency);
  double SwapHigh = pricingSwap(Expiry, Tenor, 0.01, PayFrequency);
  double Level = (SwapLow - SwapHigh) / 0.01;
  double SwapRate = SwapLow / Level;
  double swaptionATM = pricingSwaption(Expiry, Tenor, SwapRate, PayFrequency);
  double Lo = 0.001; // 0.1 PP
  double Hi = 4.999; // 500 PP
  double Vol = 0.5 * (Lo + Hi);
  double swaptionVal =
      Level * SwapRate * (2.0 * cdf_normal(0.5 * Vol * sqrt(Expiry)) - 1.0);
  int Niter = 0;
  do {
    Niter++;
    if (swaptionVal < swaptionATM)
      Lo = Vol;
    if (swaptionVal > swaptionATM)
      Hi = Vol;
    Vol = 0.5 * (Lo + Hi);
    swaptionVal =
        Level * SwapRate * (2.0 * cdf_normal(0.5 * Vol * sqrt(Expiry)) - 1.0);
  } while (fabs(1.0 - swaptionVal / swaptionATM) >= 1.0E-6 && Niter < 13);
  return Vol;
};
} // namespace velesquant