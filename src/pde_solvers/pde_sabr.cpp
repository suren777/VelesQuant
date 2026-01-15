#include <velesquant/pde_solvers/cyclic_reduction.h>
#include <velesquant/pde_solvers/pde_sabr.h>
#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#endif
#include <velesquant/volatility/lm.h>

#include <boost/bind.hpp>
#include <omp.h>
#include <ql/errors.hpp>
using namespace std;
#include <algorithm>

namespace velesquant {

void pdeSABR::setZbounds() {
  zmax_ = zmin_ = -3 * sqrt(T_);
  zmax_ *= -1;
  if (beta_ < 1) {
    double aux = -pow(f_ + shift_, 1 - beta_) * nu_ / alpha_;
    double zb = -1 / nu_ *
                log((sqrt(1 - rho_ * rho_ + (rho_ + aux) * (rho_ + aux)) -
                     rho_ - nu_ * aux / alpha_) /
                    (1 - rho_));
    zmin_ = (zb > zmin_) ? zb : zmin_;
  }
}

double pdeSABR::Y(double z) {
  return alpha_ / nu_ * (sinh(nu_ * z) + rho_ * (cosh(nu_ * z) - 1));
};
vector<double> pdeSABR::Y(vector<double> z) {
  vector<double> y;
  for (size_t i = 0; i + 1 < z.size(); i++)
    y.push_back(Y(z[i]));
  return y;
}

double pdeSABR::F(double z) {
  return pow(Lm(f_) + (1 - beta_) * z, 1.0 / (1 - beta_));
};
vector<double> pdeSABR::F(vector<double> z) {
  vector<double> y;
  for (size_t i = 0; i + 1 < z.size(); i++)
    y.push_back(F(z[i]));
  return y;
}

double pdeSABR::C(double y, double f) {
  return sqrt(alpha_ * alpha_ + 2 * rho_ * alpha_ * nu_ * y +
              nu_ * nu_ * y * y) *
         L(f);
};
vector<double> pdeSABR::C(vector<double> y, vector<double> f) {
  vector<double> aux;
  for (size_t i = 0; i + 1 < y.size(); i++)
    y.push_back(C(y[i], f[i]));
  return aux;
}

double pdeSABR::G(double F) {
  return (F == f_) ? beta_ / Lm(f_) : (L(F) - L(f_)) / (F - f_);
};
vector<double> pdeSABR::G(vector<double> z) {
  vector<double> y;
  for (size_t i = 0; i + 1 < z.size(); i++)
    y.push_back(F(z[i]));
  return y;
}

double pdeSABR::L(double F) { return pow(F + shift_, beta_); }
double pdeSABR::Lm(double F) { return pow(F + shift_, 1 - beta_); }
double pdeSABR::Mshifted(double F, double t) {

  t = 0;
  double l = L(F);
  double z = (Lm(F) - Lm(f_)) / (eps_ * alpha_ * (1 - beta_));
  double G = (l - L(f_)) / (F - f_);
  //	double ex = exp(eps_*eps_*rho_*nu_*alpha_*G*(T_-t));
  double ex = exp(eps_ * eps_ * rho_ * nu_ * alpha_ * G * t);
  double nex =
      0.5 * eps_ * eps_ * alpha_ * alpha_ *
      (1 + 2 * eps_ * rho_ * nu_ * z + eps_ * eps_ * nu_ * nu_ * z * z);
  return nex * ex * l * l;
};

double pdeSABR::Mnormal(double F, double /*t*/) {
  return 0.5 * eps_ * eps_ *
         (alpha_ * alpha_ + 2 * rho_ * nu_ * alpha_ * (F - f_) +
          nu_ * nu_ * (F - f_) * (F - f_));
}

double pdeSABR::volNormal(double K) {
  double zeta = nu_ / alpha_ * (f_ - K);
  double x =
      log((sqrt(1 - rho_ * zeta + zeta * zeta) - rho_ + zeta) / (1 - rho_));
  return eps_ * alpha_ * zeta / x *
         (1 + (.25 * rho_ * nu_ * alpha_ +
               1 / 24.0 * (2 - 3 * rho_ * rho_) * nu_ * nu_) *
                  eps_ * eps_ * T_);
};

double pdeSABR::volShifted(double K) {
  double fbk = pow(f_ + shift_, beta_) - pow(K + shift_, beta_);
  double fmbk = pow(f_ + shift_, 1 - beta_) - pow(K + shift_, 1 - beta_);
  double lfk = log((f_ + shift_) / (K + shift_));
  double z = nu_ / alpha_ * fmbk / (1 - beta_);
  double x = log((sqrt(1 - rho_ * z + z * z) - rho_ + z) / (1 - rho_));
  double p1 = eps_ * alpha_ * (f_ - K) * (1 - beta_) / fmbk * z / x;
  double p2 = -1.0 / 24 *
              (beta_ * (2 - beta_) * (1 - beta_) * (1 - beta_) * alpha_ *
               alpha_ * lfk * lfk) /
              (fmbk * fmbk);
  double p3 = 0.25 * rho_ * nu_ * alpha_ * fbk / (f_ - K) +
              (2 - 3 * rho_ * rho_) / 24 * nu_ * nu_;
  return p1 * (1 + (p2 + p3) * eps_ * eps_ * T_);
};

double pdeSABR::volShifted1(double K) {
  double lfbk = log(f_ + shift_) + log(K + shift_);
  double z = nu_ / alpha_ * lfbk;
  double x = log((sqrt(1 - rho_ * z + z * z) - rho_ + z) / (1 - rho_));
  double p1 = eps_ * alpha_ * (f_ - K) / lfbk * z / x;
  double p3 = 0.25 * rho_ * nu_ * alpha_ +
              (2 - 3 * rho_ * rho_) / 24 * nu_ * nu_ -
              1.0 / 24 * alpha_ * alpha_;
  return p1 * (1.0 + p3 * eps_ * eps_ * T_);
};

void pdeSABR::calibrator(vector<volQuote> &quotes) {
  int n = 3;
  quotes_ = quotes;
  int m = quotes_.size();
  double *x = new double[n]; // initial estimate of parameters vector
  x[0] = alpha_;
  x[1] = rho_;
  x[2] = nu_;

  QL_ENSURE(m >= n, "too much freedom in Calibration  " << m - n);
  double *fvec = new double[m]; // no need to populate
  double ftol = 1e-10;          // tolerance
  double xtol = 1e-10;          // tolerance
  double gtol = 1e-10;          // tolerance
  int maxfev = 10000;           // maximum function evaluations
  double epsfcn = 1e-10;        // tolerance
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
  boost::function<void(pdeSABR *, int, int, double *, double *, int *)> obj =
      &pdeSABR::objFcnCalibration;
  lmfcn fcn = boost::bind(obj, this, _1, _2, _3, _4, _5);
  lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
        nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);

  QL_ENSURE(info != 4, "SABR Model Calibration Fails " << info);
  // the below is output result

  alpha_ = x[0];
  rho_ = sin(x[1]);
  nu_ = x[2];

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

void pdeSABR::objFcnCalibration(int m, int /*n*/, double *x, double *fvec,
                                int * /*iflag*/) {
  alpha_ = x[0];
  rho_ = sin(x[1]);
  nu_ = x[2];
  double penalty = ((alpha_ <= 0) ? 1e8 * alpha_ * alpha_ : 0) +
                   ((nu_ <= 0) ? 1e8 * nu_ * nu_ : 0);
  double error;
  for (int i = 0; i < m; i++) {
    error = quotes_[i].IV - sabrVol(quotes_[i].Strike);
    fvec[i] = error + penalty;
  }
};

void pdeSABR::oneStepForward(double t0, vector<double> &inV, pde_func &Mf) {
  vector<double> l(sizeX_), c(sizeX_), u(sizeX_), d(sizeX_);
  double t1 = (t0 + .5) * dT_;
  vector<double> aux(inV.begin() + 1, inV.end() - 1);
  double ndt = 2 * dT_ / h_ / h_;

#pragma omp parallel for
  for (int r = 1; r <= sizeX_; r++) {
    if (r == 1) {
      l[r - 1] = 0;
      u[r - 1] = -ndt * Mf((r + 1) * h_, t1);
      c[r - 1] = 3 * ndt * Mf((r)*h_, t1) + 1.0;
      d[r - 1] = inV[r] + ndt * (-2 * Mf((r - 0) * h_, t0) * inV[r - 0] +
                                 Mf((r + 1) * h_, t0) * inV[r + 1]);
    } else if (r == sizeX_ - 1) {
      u[r - 1] = 0;
      l[r - 1] = -ndt * Mf((r - 1) * h_, t1);
      c[r - 1] = 3 * ndt * Mf((r)*h_, t1) + 1.0;
      d[r - 1] = inV[r] + ndt * (Mf((r - 1) * h_, t0) * inV[r - 1] -
                                 3 * Mf((r - 0) * h_, t0) * inV[r - 0]);
    } else {
      d[r - 1] = inV[r] + ndt * (Mf((r - 1) * h_, t0) * inV[r - 1] -
                                 2 * Mf((r - 0) * h_, t0) * inV[r - 0] +
                                 Mf((r + 1) * h_, t0) * inV[r + 1]);
      l[r - 1] = -ndt * Mf((r - 1) * h_, t1);
      c[r - 1] = 2 * ndt * Mf((r)*h_, t1) + 1.0;
      u[r - 1] = -ndt * Mf((r + 1) * h_, t1);
    }
  }
  TriDiagonalSolve(sizeX_, l, c, u, d, aux);
  inV[0] = inV[0] + dT_ / h_ * (Mf(h_, t1) * inV[1] + Mf(h_, t0) * aux[0]);
  inV[sizeX_ - 1] =
      inV[sizeX_ - 1] - dT_ / h_ *
                            (Mf(h_ * sizeX_, t1) * inV[sizeX_ - 2] +
                             Mf(h_ * sizeX_, t0) * aux[sizeX_ - 1]);
#pragma omp parallel for
  for (int i = 0; i < sizeX_ - 1; i++)
    inV[i + 1] = aux[i];
}

void pdeSABR::oneStepForwardN(vector<double> &inV) {
  vector<double> l(sizeX_ + 2), c(sizeX_ + 2), u(sizeX_ + 2), d(sizeX_ + 2);
  vector<double> aux(inV.begin(), inV.end());
  double ndt = dT_ / (2 * h_);
#pragma omp parallel for
  for (int i = 0; i <= sizeX_ + 1; i++) {

    u[i] = (i == sizeX_ + 1)
               ? 0
               : -ndt * Cm_[i + 1] / (Fm_[i + 1] - Fm_[i]) * Em_[i + 1];
    l[i] =
        (i == 0) ? 0 : -ndt * Cm_[i - 1] / (Fm_[i] - Fm_[i - 1]) * Em_[i - 1];
    if (i == 0)
      c[i] = Cm_[0] / (Fm_[1] - Fm_[0]) * Em_[0];
    else if (i == sizeX_ + 1)
      c[i] = Cm_[i] / (Fm_[i] - Fm_[i - 1]) * Em_[i];
    else
      c[i] = ndt * Em_[i] * Cm_[i] *
             (1.0 / (Fm_[i + 1] - Fm_[i]) + 1.0 / (Fm_[i] - Fm_[i - 1]));
  }
  TriDiagonalSolve(sizeX_, l, c, u, d, inV);
  inV[0] = inV[sizeX_ + 1] = 0;
  QL_ += dT_ * Cm_[1] / (Fm_[1] - Fm_[0]) * Em_[1] * inV[1];
  QR_ += dT_ * Cm_[sizeX_] / (Fm_[sizeX_ + 1] - Fm_[sizeX_]) * Em_[sizeX_] *
         inV[sizeX_];
};

//
// void pdeSABR::oneStepForward(double t0, vector<double>&inV)
//{
//	vector<double> l(sizeX_),c(sizeX_),u(sizeX_),d(sizeX_);
//	double t1 = (t0+.5)*dT_;
//	vector<double>aux(inV.begin()+1,inV.end()-1);
//	double ndt = 2*dT_/h_/h_;
//
// #pragma omp parallel for
//	for(int r=1; r<=sizeX_; r++)
//	{
//		if (r==1)
//		{
//			l[r-1]=0;
//			u[r-1]=-ndt*Mgrid_[r];
//			c[r-1] = 3*ndt*Mgrid_[r-1]+1.0;
//			d[r-1] =
// inV[r]+ndt*(-3*Mgrid_[r-1]*inV[r]+Mgrid_[r]*inV[r+1]);
//		}
//		else if (r==sizeX_)
//		{
//			u[r-1]=0;
//			l[r-1] =-ndt*Mgrid_[r-2];
//			c[r-1] = 3*ndt*Mgrid_[r-1]+1.0;
//			d[r-1] =
// inV[r]+ndt*(Mgrid_[r-2]*inV[r-1]-3*Mgrid_[r-1]*inV[r]);
//		}
//		else
//		{
//			d[r-1] =
// inV[r]+ndt*(Mgrid_[r-2]*inV[r-1]-2*Mgrid_[r-1]*inV[r]+Mgrid_[r]*inV[r+1]);
//			l[r-1] = -ndt*Mgrid_[r-2];
//			c[r-1] = 2*ndt*Mgrid_[r-1]+1.0;
//			u[r-1] = -ndt*Mgrid_[r];
//
//		}
//	}
//	TriDiagonalSolve(sizeX_,l,c,u,d,aux);
//	inV[0]+=		dT_/h_*Mgrid_[0]*(inV[1]+aux[0]);
//	inV[sizeX_+1]+=	dT_/h_*Mgrid_[sizeX_-1]*(inV[sizeX_]+aux[sizeX_-1]);
// #pragma omp parallel for
//	for(int i = 0; i < sizeX_;i++) inV[i+1] = aux[i];
//}

vector<double> pdeSABR::calculateDensity() {
  dT_ = T_ / sizeT_;
  setFbounds();
  pde_func Mf;
  if (type_ == "Normal")
    Mf = boost::bind(&pdeSABR::Mnormal, this, _1, _2);
  else
    Mf = boost::bind(&pdeSABR::Mshifted, this, _1, _2);
  vector<double> inV(sizeX_ + 2, 0.0);
  inV[int((f_ - Fmin_) / h_ + 0.5)] = 1 / h_;
  for (int t = 0; t < sizeT_; t++)
    oneStepForward(t * dT_, inV, Mf);
  return inV;
};
void pdeSABR::setDensity() { Density_ = calculateDensity(); }
double pdeSABR::sabrVol(double K) {
  if (beta_ == 0.0)
    if (K == f_)
      return (volNormal(K * .99999) + volNormal(K * 1.00001)) / 2;
    else
      return volNormal(K);
  else if (beta_ == 1.0)
    if (K == f_)
      return (volShifted1(K * .99999) + volShifted1(K * 1.00001)) / 2;
    else
      return volShifted1(K);
  else if (K == f_)
    return (volShifted(K * .99999) + volShifted(K * 1.00001)) / 2;
  else
    return volShifted(K);
}

void pdeSABR::setFbounds() {
  int N = 6;
  double en = eps_ * nu_;
  double theta = en / 2 * N * sqrt(T_);
  double fup = 2 / en * sinh(theta) * (cosh(theta) + rho_ * sinh(theta));
  double fdn = 2 / en * sinh(theta) * (cosh(theta) - rho_ * sinh(theta));
  Fmax_ = pow(fup * eps_ * alpha_ * (1 - beta_) + pow(f_ + shift_, 1 - beta_),
              1.0 / (1 - beta_)) -
          shift_;
  Fmin_ = pow(-fdn * eps_ * alpha_ * (1 - beta_) + pow(f_ + shift_, 1 - beta_),
              1.0 / (1 - beta_)) -
          shift_;
  if (Fmin_ >= Fmax_)
    throw("Error in grid Bounds!");
  if (f_ < Fmin_)
    throw std::runtime_error("Bad grid bounds");
  if (f_ > Fmax_)
    throw std::runtime_error("Bad grid bounds");
  if ((Fmin_ != Fmin_) || (Fmax_ != Fmax_))
    throw("Illegal bounds, ajust shift and recalibrate");
  h_ = (Fmax_ - Fmin_) / sizeX_;
  Fgrid_.resize(sizeX_ + 2);
  Fgrid_[0] = Fmin_;
  Fgrid_[sizeX_ + 1] = Fmax_;
  Mgrid_.resize(sizeX_, 0.0);
  pde_func Mf;
  // if (type_=="Normal") Mf = boost::bind(&pdeSABR::Mnormal, this, _1,0);
  // else Mf = boost::bind(&pdeSABR::Mshifted, this, _1,0);
#pragma omp parallel for

  for (int i = 0; i < sizeX_; i++) {
    double aux = (Fmin_ + h_ / 2) + i * h_;
    Fgrid_[i + 1] = aux;
    //		Mgrid_[i]=Mf(aux);
  }
};

double pdeSABR::sabr_option(double strike, const string &type) {
  double payoff = 0.0;
  if (Density_.size() <= 1)
    setDensity();
  if (type == "Call") {
    if (strike < Fmin_)
      payoff = f_ - strike;
    else if ((strike <= Fmax_) && (strike >= Fmin_)) {
      int i = 0;
      for (i = 1; i <= sizeX_; i++)
        if ((Fgrid_[i] >= strike) && (Fgrid_[i + 1] < strike))
          break;
      do {
        payoff += (Fgrid_[i] - strike) * Density_[i] * h_;
      } while (i <= sizeX_);
      payoff += (Fmax_ - strike) * Density_[sizeX_ + 1];
    }
  } else {
    if (strike > Fmax_)
      payoff = f_ - strike;
    for (int i = 1; ((Fmin_ + h_ / 2) + i * h_) < strike; i++) {
      payoff += (i >= sizeX_) ? 0 : (strike - Fgrid_[i]) * Density_[i] * h_;
    }
    payoff += (strike - Fmin_) * Density_[0];
  }
  return payoff;
};
} // namespace velesquant