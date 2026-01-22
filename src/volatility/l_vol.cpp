#include <algorithm>
#include <assert.h>
#include <cmath>

#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <ql/quantlib.hpp>
#include <velesquant/models/utility.h>
#include <velesquant/numerics/tri_diag_matrix.h>
#include <velesquant/volatility/l_vol.h>
#include <velesquant/volatility/sabr.h>

// using namespace std;

namespace velesquant {

lVol::lVol(std::vector<double> maturities, std::vector<double> forwards,
           std::vector<double> betas, std::vector<double> alphas,
           std::vector<double> nus, std::vector<double> rhos, double spot)
    : spot_(spot) {
  int N = maturities.size();
  sabrModels_.resize(N);
  for (int i = 0; i < N; i++) {
    Sabr sabrModel(maturities[i], forwards[i], betas[i], alphas[i], nus[i],
                   rhos[i]);
    sabrModels_[i] = sabrModel;
  }
};

static const int Mx = 137;
static const double Wx[Mx] = {
    0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.225, 0.250, 0.275, 0.300, 0.320,
    0.340, 0.360, 0.380, 0.400, 0.420, 0.440, 0.460, 0.480, 0.500, 0.520, 0.540,
    0.560, 0.580, 0.600, 0.620, 0.640, 0.660, 0.680, 0.700, 0.710, 0.720, 0.730,
    0.740, 0.750, 0.760, 0.770, 0.780, 0.790, 0.800, 0.810, 0.820, 0.830, 0.840,
    0.850, 0.860, 0.870, 0.880, 0.890, 0.900, 0.905, 0.910, 0.915, 0.920, 0.925,
    0.930, 0.935, 0.940, 0.945, 0.950, 0.955, 0.960, 0.965, 0.970, 0.975, 0.980,
    0.985, 0.990, 0.995, 1.000, 1.005, 1.010, 1.015, 1.020, 1.025, 1.030, 1.035,
    1.040, 1.045, 1.050, 1.055, 1.060, 1.065, 1.070, 1.075, 1.080, 1.085, 1.090,
    1.095, 1.100, 1.110, 1.120, 1.130, 1.140, 1.150, 1.160, 1.170, 1.180, 1.190,
    1.200, 1.210, 1.220, 1.230, 1.240, 1.250, 1.260, 1.270, 1.280, 1.290, 1.300,
    1.320, 1.340, 1.360, 1.380, 1.400, 1.420, 1.440, 1.460, 1.480, 1.500, 1.525,
    1.550, 1.575, 1.600, 1.650, 1.700, 1.750, 1.800, 1.850, 1.900, 1.950, 2.000,
    2.250, 2.500, 3.000, 4.000, 5.000};

std::vector<std::vector<double>>
lVol::exportLV(const std::vector<double> times) {
  buildGrid(times);
  return gridLV_;
};

void lVol::buildGrid(const std::vector<double> times) {
  int Nt = times.size();
  gridT_.resize(Nt);
  gridFwd_.resize(Nt);
  gridLV_.resize(Nt);
  for (int t = 0; t < Nt; t++) {
    gridT_[t] = times[t];
    gridFwd_[t].resize(Mx);
    gridLV_[t].resize(Mx);

    double currentTime = gridT_[t];
    double currentForward = getForward(currentTime);
    for (int x = 0; x < Mx; x++) {
      gridFwd_[t][x] = Wx[x] * currentForward;
      currentTime = std::max(0.002, currentTime); // one day 1/365=0.002
      gridLV_[t][x] = getLocVol(currentTime, gridFwd_[t][x], currentForward);
    }
  }
};

void lVol::buildGrid(double lastTime, int Nt) {
  std::vector<double> times(Nt);
  double delT = lastTime / Nt;
  times[Nt - 1] = lastTime;
  for (int i = 0; i < Nt - 1; i++)
    times[i] = (i + 1) * delT;
  buildGrid(times);
};

double lVol::callPDE(double maturity, double strike, const int Nt) {
  buildGrid(maturity, Nt);
  std::vector<double> payoff(Mx), pv(Mx);
  double fwd = getForward(maturity);
  for (int x = 0; x < Mx; x++)
    payoff[x] = std::max(0.0, Wx[x] * fwd - strike);
  for (int t = Nt - 2; t >= 0; t--) {
    oneStepBackward(t, payoff, pv);
    for (int x = 0; x < Mx; x++)
      payoff[x] = pv[x];
  }
  int x;
  for (x = 0; x < Mx; x++)
    if (Wx[x] * spot_ >= spot_)
      break;
  return pv[x];
};
double lVol::putPDE(double maturity, double strike, const int Nt) {
  buildGrid(maturity, Nt);
  std::vector<double> payoff(Mx), pv(Mx);
  double fwd = getForward(maturity);
  for (int x = 0; x < Mx; x++)
    payoff[x] = std::max(0.0, strike - Wx[x] * fwd);
  for (int t = Nt - 2; t >= 0; t--) {
    oneStepBackward(t, payoff, pv);
    for (int x = 0; x < Mx; x++)
      payoff[x] = pv[x];
  }
  int x;
  for (x = 0; x < Mx; x++)
    if (Wx[x] * spot_ >= spot_)
      break;
  return pv[x];
};

void lVol::oneStepBackward(const int t, const std::vector<double> &inV,
                           std::vector<double> &outV) {
  std::vector<double> l(Mx - 2), c(Mx - 2), u(Mx - 2), d(Mx - 2), V(Mx - 2);
  double delT = gridT_[t + 1] - gridT_[t];
  // x in the middle range
  for (int r = 1; r < Mx - 1; r++) {
    double R = gridFwd_[t][r];
    double Ru = gridFwd_[t][r + 1];
    double Rl = gridFwd_[t][r - 1];
    double difu =
        0.5 * gridLV_[t][r] * gridLV_[t][r] * gridFwd_[t][r] * gridFwd_[t][r];
    l[r - 1] = -difu / (R - Rl) / (Ru - Rl);
    c[r - 1] = 1.0 / delT + difu / (R - Rl) / (Ru - R);
    u[r - 1] = -difu / (Ru - R) / (Ru - Rl);
    R = gridFwd_[t + 1][r];
    Ru = gridFwd_[t + 1][r + 1];
    Rl = gridFwd_[t + 1][r - 1];
    difu = 0.5 * gridLV_[t + 1][r] * gridLV_[t + 1][r] * gridFwd_[t + 1][r] *
           gridFwd_[t + 1][r];
    d[r - 1] = difu / (R - Rl) / (Ru - Rl) * inV[r - 1] +
               (1.0 / delT - difu / (R - Rl) / (Ru - R)) * inV[r] +
               difu / (Ru - R) / (Ru - Rl) * inV[r + 1];
  }
  // r=0 lower boundary condition
  double R0 = gridFwd_[t][0];
  double R1 = gridFwd_[t][1];
  double R2 = gridFwd_[t][2];
  double difu =
      0.5 * gridLV_[t][0] * gridLV_[t][0] * gridFwd_[t][0] * gridFwd_[t][0];
  double k0 = 1.0 / delT - difu / (R1 - R0) / (R2 - R0);
  double k1 = difu / (R1 - R0) / (R2 - R1);
  double k2 = -difu / (R2 - R1) / (R2 - R0);
  R0 = gridFwd_[t + 1][0];
  R1 = gridFwd_[t + 1][1];
  R2 = gridFwd_[t + 1][2];
  difu = 0.5 * gridLV_[t + 1][0] * gridLV_[t + 1][0] * gridFwd_[t + 1][0] *
         gridFwd_[t + 1][0];
  double d0 = (1.0 / delT + difu / (R1 - R0) / (R2 - R0)) * inV[0] -
              difu / (R1 - R0) / (R2 - R1) * inV[1] +
              difu / (R2 - R1) / (R2 - R0) * inV[2];
  c[0] = c[0] - l[0] * k1 / k0;
  u[0] = u[0] - l[0] * k2 / k0;
  d[0] = d[0] - l[0] * d0 / k0;
  l[0] = 0.0;
  // r=MR-1 upper boundary condition
  double R3 = gridFwd_[t][Mx - 3];
  R2 = gridFwd_[t][Mx - 2];
  R1 = gridFwd_[t][Mx - 1];
  difu = 0.5 * gridLV_[t][Mx - 1] * gridLV_[t][Mx - 1] * gridFwd_[t][Mx - 1] *
         gridFwd_[t][Mx - 1];
  double kl3 = -difu / (R2 - R3) / (R1 - R3);
  double kl2 = difu / (R2 - R3) / (R1 - R2);
  double kl1 = 1.0 / delT - difu / (R1 - R2) / (R1 - R3);
  R3 = gridFwd_[t + 1][Mx - 3];
  R2 = gridFwd_[t + 1][Mx - 2];
  R1 = gridFwd_[t + 1][Mx - 1];
  difu = 0.5 * gridLV_[t + 1][Mx - 1] * gridLV_[t + 1][Mx - 1] *
         gridFwd_[t + 1][Mx - 1] * gridFwd_[t + 1][Mx - 1];
  double dl1 = difu / (R2 - R3) / (R1 - R3) * inV[Mx - 3] -
               difu / (R2 - R3) / (R1 - R2) * inV[Mx - 2] +
               (1.0 / delT + difu / (R1 - R2) / (R1 - R3)) * inV[Mx - 1];
  l[Mx - 3] = l[Mx - 3] - u[Mx - 3] * kl3 / kl1;
  c[Mx - 3] = c[Mx - 3] - u[Mx - 3] * kl2 / kl1;
  d[Mx - 3] = d[Mx - 3] - u[Mx - 3] * dl1 / kl1;
  u[Mx - 3] = 0.0;
  TriDiagonalSolve(Mx - 2, l, c, u, d, V);
  for (int r = 1; r < Mx - 1; r++)
    outV[r] = V[r - 1];
  // update r=0 lower boundary
  outV[0] = (d0 - k1 * outV[1] - k2 * outV[2]) / k0;
  // update r=MR-1 upper boundary
  outV[Mx - 1] = (dl1 - kl3 * outV[Mx - 3] - kl2 * outV[Mx - 2]) / kl1;
};

void lVol::backwardOneStep(const int t, std::vector<double> &y,
                           const std::vector<double> &z) {
  std::vector<double> a(Mx), b(Mx), c(Mx), d(Mx);
  double delT = gridT_[t + 1] - gridT_[t];
  // lower boundary condition
  d[0] = z[0];
  c[0] = 0.0;
  b[0] = 1.0;
  a[0] = 0.0;
  // middle range
  for (int x = 1; x < Mx - 1; x++) {
    double delFu = gridFwd_[t][x + 1] - gridFwd_[t][x];
    double delFm = 0.5 * (gridFwd_[t][x + 1] - gridFwd_[t][x - 1]);
    double delFl = gridFwd_[t][x] - gridFwd_[t][x - 1];
    double var = 0.5 * delT * gridLV_[t][x] * gridLV_[t][x] * gridFwd_[t][x] *
                 gridFwd_[t][x];
    double su = 0.5 * var / delFu / delFm;
    double sl = 0.5 * var / delFl / delFm;
    c[x] = -su;
    b[x] = 1.0 + su + sl;
    a[x] = -sl;
    double delFu2 = gridFwd_[t + 1][x + 1] - gridFwd_[t + 1][x];
    double delFm2 = 0.5 * (gridFwd_[t + 1][x + 1] - gridFwd_[t + 1][x - 1]);
    double delFl2 = gridFwd_[t + 1][x] - gridFwd_[t + 1][x - 1];
    double var2 = 0.5 * delT * gridLV_[t + 1][x] * gridLV_[t + 1][x] *
                  gridFwd_[t + 1][x] * gridFwd_[t + 1][x];
    double su2 = 0.5 * var2 / delFu2 / delFm2;
    double sl2 = 0.5 * var2 / delFl2 / delFm2;
    d[x] = su2 * z[x + 1] + (1.0 - su2 - sl2) * z[x] + sl2 * z[x - 1];
  }
  // upper boundary condition
  d[Mx - 1] = z[Mx - 1];
  c[Mx - 1] = 0.0;
  b[Mx - 1] = 1.0;
  a[Mx - 1] = 0.0;
  TriDiagonalSolve(Mx, a, b, c, d, y);
};

double lVol::dntPDE(double maturity, double upperBarrier, double lowerBarrier,
                    const int Nt) {
  gridT_.resize(Nt);
  gridFwd_.resize(Nt);
  gridLV_.resize(Nt);
  double delT = maturity / Nt;
  for (int t = 0; t < Nt; t++) {
    double currentTime = (t + 1) * delT;
    gridT_[t] = currentTime;
    gridFwd_[t].push_back(lowerBarrier);
    double currentForward = getForward(currentTime);
    double localV = getLocVol(currentTime, lowerBarrier, currentForward);
    gridLV_[t].push_back(localV);
    for (int x = 0; x < Mx; x++) {
      double spot = Wx[x] * currentForward;
      if (spot > lowerBarrier && spot < upperBarrier) {
        gridFwd_[t].push_back(spot);
        localV = getLocVol(currentTime, spot, currentForward);
        gridLV_[t].push_back(localV);
      }
    }
    gridFwd_[t].push_back(upperBarrier);
    localV = getLocVol(currentTime, upperBarrier, currentForward);
    gridLV_[t].push_back(localV);
  };
  int M = gridFwd_[Nt - 1].size();
  std::vector<double> payoff(M), pv(M);
  for (int x = 0; x < M; x++) {
    payoff[x] = 1.0;
    double spot = gridFwd_[Nt - 1][x];
    if (spot <= lowerBarrier || spot >= upperBarrier)
      payoff[x] = 0.0;
  }
  for (int t = Nt - 2; t >= 0; t--) {
    oneStepBackwardDirichlet(t, payoff, pv);
    payoff = pv;
  }
  int x;
  for (x = 0; x < M; x++)
    if (gridFwd_[0][x] >= spot_)
      break;
  return pv[x];
};

void lVol::oneStepBackwardDirichlet(const int t, const std::vector<double> &inV,
                                    std::vector<double> &outV) {
  int M = inV.size();
  std::vector<double> l(M - 2), c(M - 2), u(M - 2), d(M - 2), V(M - 2);
  double delT = gridT_[t + 1] - gridT_[t];
  // x in the middle range
  for (int r = 1; r < M - 1; r++) {
    double R = gridFwd_[t][r];
    double Ru = gridFwd_[t][r + 1];
    double Rl = gridFwd_[t][r - 1];
    double difu =
        0.5 * gridLV_[t][r] * gridLV_[t][r] * gridFwd_[t][r] * gridFwd_[t][r];
    l[r - 1] = -difu / (R - Rl) / (Ru - Rl);
    c[r - 1] = 1.0 / delT + difu / (R - Rl) / (Ru - R);
    u[r - 1] = -difu / (Ru - R) / (Ru - Rl);
    R = gridFwd_[t + 1][r];
    Ru = gridFwd_[t + 1][r + 1];
    Rl = gridFwd_[t + 1][r - 1];
    difu = 0.5 * gridLV_[t + 1][r] * gridLV_[t + 1][r] * gridFwd_[t + 1][r] *
           gridFwd_[t + 1][r];
    d[r - 1] = difu / (R - Rl) / (Ru - Rl) * inV[r - 1] +
               (1.0 / delT - difu / (R - Rl) / (Ru - R)) * inV[r] +
               difu / (Ru - R) / (Ru - Rl) * inV[r + 1];
  }
  // r=0 lower boundary condition
  d[0] = d[0] + l[0] * inV[0];
  l[0] = 0.0;
  // r=MR-1 upper boundary condition
  d[M - 3] = d[M - 3] + u[M - 3] * inV[M - 1];
  u[M - 3] = 0.0;
  TriDiagonalSolve(M - 2, l, c, u, d, V);
  for (int r = 1; r < M - 1; r++)
    outV[r] = V[r - 1];
  // update r=0 lower boundary
  outV[0] = 0.0;
  // update r=MR-1 upper boundary
  outV[M - 1] = 0.0;
};

void lVol::forwardOneStep(const int t, const std::vector<double> &y,
                          std::vector<double> &z) {
  std::vector<double> a(Mx), b(Mx), c(Mx), d(Mx);
  double delT = gridT_[t] - gridT_[t - 1];
  // lower boundary condition
  d[0] = y[0];
  c[0] = 0.0;
  b[0] = 1.0;
  a[0] = 0.0;
  // middle range
  for (int x = 1; x < Mx - 1; x++) {
    double delFu = gridFwd_[t][x + 1] - gridFwd_[t][x];
    double delFm = 0.5 * (gridFwd_[t][x + 1] - gridFwd_[t][x - 1]);
    double delFl = gridFwd_[t][x] - gridFwd_[t][x - 1];
    double varu = 0.5 * delT * gridLV_[t][x + 1] * gridLV_[t][x + 1] *
                  gridFwd_[t][x + 1] * gridFwd_[t][x + 1];
    double varm = 0.5 * delT * gridLV_[t][x] * gridLV_[t][x] * gridFwd_[t][x] *
                  gridFwd_[t][x];
    double varl = 0.5 * delT * gridLV_[t][x - 1] * gridLV_[t][x - 1] *
                  gridFwd_[t][x - 1] * gridFwd_[t][x - 1];
    double su = 0.5 * varu / delFu / delFm;
    double sm = 0.5 * varm * (1.0 / delFu + 1.0 / delFl) / delFm;
    double sl = 0.5 * varl / delFl / delFm;
    c[x] = -su;
    b[x] = 1.0 + sm;
    a[x] = -sl;
    double delFu2 = gridFwd_[t - 1][x + 1] - gridFwd_[t - 1][x];
    double delFm2 = 0.5 * (gridFwd_[t - 1][x + 1] - gridFwd_[t - 1][x - 1]);
    double delFl2 = gridFwd_[t - 1][x] - gridFwd_[t - 1][x - 1];
    double varu2 = 0.5 * delT * gridLV_[t - 1][x + 1] * gridLV_[t - 1][x + 1] *
                   gridFwd_[t - 1][x + 1] * gridFwd_[t - 1][x + 1];
    double varm2 = 0.5 * delT * gridLV_[t - 1][x] * gridLV_[t - 1][x] *
                   gridFwd_[t - 1][x] * gridFwd_[t - 1][x];
    double varl2 = 0.5 * delT * gridLV_[t - 1][x - 1] * gridLV_[t - 1][x - 1] *
                   gridFwd_[t - 1][x - 1] * gridFwd_[t - 1][x - 1];
    double su2 = 0.5 * varu2 / delFu2 / delFm2;
    double sm2 = 0.5 * varm2 * (1.0 / delFu2 + 1.0 / delFl2) / delFm2;
    double sl2 = 0.5 * varl2 / delFl2 / delFm2;
    d[x] = su2 * y[x + 1] + (1.0 - sm2) * y[x] + sl2 * y[x - 1];
  }
  // upper boundary condition
  d[Mx - 1] = y[Mx - 1];
  c[Mx - 1] = 0.0;
  b[Mx - 1] = 1.0;
  a[Mx - 1] = 0.0;
  TriDiagonalSolve(Mx, a, b, c, d, z);
};

std::vector<double> lVol::density(double maturity, const int Nt) {
  buildGrid(maturity, Nt);
  std::vector<double> idensity(Mx), tdensity(Mx);
  for (int x = 0; x < Mx; x++) {
    idensity[x] = 0.0;
    if (Wx[x] * spot_ >= spot_ && Wx[x - 1] * spot_ < spot_)
      idensity[x] = 1.0 / ((Wx[x + 1] - Wx[x - 1]) / 2.0 * spot_);
  }
  for (int t = 1; t < Nt; t++) {
    forwardOneStep(t, idensity, tdensity);
    for (int x = 0; x < Mx; x++)
      idensity[x] = tdensity[x];
  }

  double pr = 0.0;
  double fwd = getForward(maturity);
  for (int x = 1; x < Mx - 1; x++)
    pr += tdensity[x] * (0.5 * (Wx[x + 1] - Wx[x - 1]) * fwd);
  tdensity.push_back(pr);

  return tdensity;
};

double lVol::premiumCALL(double strike, double forward, double maturity,
                         double vol) const {
  double deviation = vol * std::sqrt(maturity);
  double d1 = std::log(forward / strike) / deviation + 0.5 * deviation;
  double d2 = d1 - deviation;
  double premiumCALL;
  if (forward > strike) {
    if ((d1 != d1) || (d2 != d2))
      premiumCALL = forward - strike;
    else
      premiumCALL = forward * cdf_normal(d1) - strike * cdf_normal(d2);
  } else {
    double premiumPUT;
    if ((d1 != d1) || (d2 != d2))
      premiumPUT = 0.0;
    else
      premiumPUT = strike * cdf_normal(-d2) - forward * cdf_normal(-d1);
    premiumCALL = premiumPUT + (forward - strike);
  }
  return premiumCALL;
};

double lVol::interpolatedCall(const Sabr &nextModel, double strike,
                              double forward, double maturity) const {
  double nextMat = nextModel.getMaturity();
  double nextFwd = nextModel.getForward();
  // double d = log(forward/strike)/sqrt(maturity);
  // double nextStrike = nextFwd * exp(-d*sqrt(nextMat));   // same delta
  double nextStrike = nextFwd * strike / forward; // same ratio - return
  double nextImpvol = nextModel.impliedVol(nextStrike);
  return premiumCALL(strike, forward, maturity, nextImpvol);
};

double lVol::interpolatedCall(const Sabr &preModel, const Sabr &nextModel,
                              double strike, double forward,
                              double maturity) const {
  double preMat = preModel.getMaturity();
  double preFwd = preModel.getForward();
  // double d = log(forward/strike)/sqrt(maturity);
  // double preStrike = preFwd * exp(-d*sqrt(preMat));  // same delta
  double preStrike = preFwd * strike / forward; // same ratio - return
  double preImpvol = preModel.impliedVol(preStrike);
  double preVar = preImpvol * preImpvol * preMat;
  double nextMat = nextModel.getMaturity();
  double nextFwd = nextModel.getForward();
  // double nextStrike = nextFwd * exp(-d*sqrt(nextMat));  // same delta
  double nextStrike = nextFwd * strike / forward; // same ratio - return
  double nextImpvol = nextModel.impliedVol(nextStrike);
  double nextVar = nextImpvol * nextImpvol * nextMat;
  double vars = preVar + std::max(0.0, nextVar - preVar) * (maturity - preMat) /
                             (nextMat - preMat); // linear var
  double vol = std::sqrt(vars / maturity);
  return premiumCALL(strike, forward, maturity, vol);
};

const double RFLO = 0.01;
const double RTOP = 50.0;
const double DAMP = 0.95;
void lVol::getLocVol(double time, double preSpot, double preFwd, double &locvol,
                     double &forward) const {
  double fixRatio =
      1.0 + DAMP * (std::min(std::max(RFLO, preSpot / preFwd), RTOP) - 1.0);
  double firstMaturity = sabrModels_[0].getMaturity();
  if (time < firstMaturity) {
    double afterForward = sabrModels_[0].getForward();
    forward = spot_ + (afterForward - spot_) * time /
                          firstMaturity; // linear interpolation
    double currSpot = forward * fixRatio;
    double cVal = interpolatedCall(sabrModels_[0], currSpot, forward, time);
    double cLft =
        interpolatedCall(sabrModels_[0], currSpot * 0.999, forward, time);
    double cRgt =
        interpolatedCall(sabrModels_[0], currSpot * 1.001, forward, time);
    double dS2val =
        (cRgt + cLft - 2.0 * cVal) / (0.001 * currSpot * 0.001 * currSpot);
    // assert(dS2val > 0.0);
    dS2val = std::max(1.0E-20, dS2val);
    double cBigT =
        interpolatedCall(sabrModels_[0], currSpot, forward, time * 1.001);
    double dTval = (cBigT - cVal) / (0.001 * time);
    // assert(dTval > 0.0);
    dTval = std::max(0.0, dTval);
    locvol = std::sqrt(2.0 * dTval / dS2val) / currSpot;
    // locvol = sabrModels_[0].localVol( afterForward * fixRatio );
    return;
  }
  if (time == firstMaturity) {
    forward = sabrModels_[0].getForward();
    locvol = sabrModels_[0].localVol(forward * fixRatio);
    return;
  }
  int N = sabrModels_.size();
  double lastMaturity = sabrModels_[N - 1].getMaturity();
  if (time > lastMaturity) {
    forward = sabrModels_[N - 1].getForward(); // flat extrapolation
    locvol = sabrModels_[N - 1].localVol(forward * fixRatio);
    return;
  }
  for (int i = 1; i < N; i++) {
    double afterMaturity = sabrModels_[i].getMaturity();
    double beforeMaturity = sabrModels_[i - 1].getMaturity();
    if (time > beforeMaturity && time < afterMaturity) {
      double beforeForward = sabrModels_[i - 1].getForward();
      double afterForward = sabrModels_[i].getForward();
      forward = beforeForward + (afterForward - beforeForward) *
                                    (time - beforeMaturity) /
                                    (afterMaturity - beforeMaturity);
      double currSpot = forward * fixRatio;
      double cVal = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                     currSpot, forward, time);
      double cLft = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                     currSpot * 0.999, forward, time);
      double cRgt = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                     currSpot * 1.001, forward, time);
      double dS2val =
          (cRgt + cLft - 2.0 * cVal) / (0.001 * currSpot * 0.001 * currSpot);
      // assert(dS2val > 0.0);
      dS2val = std::max(1.0E-20, dS2val);
      double cBigT = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                      currSpot, forward, time * 1.001);
      double dTval = (cBigT - cVal) / (0.001 * time);
      // assert(dTval > 0.0);
      dTval = std::max(0.0, dTval);
      locvol = std::sqrt(2.0 * dTval / dS2val) / currSpot;
      // if(locvol != locvol) locvol = 0.0;	//if(locvol > 111.111) locvol =
      // 111.111;
      return;
    }
    if (time == afterMaturity) {
      forward = sabrModels_[i].getForward();
      locvol = sabrModels_[i].localVol(forward * fixRatio);
      return;
    }
  }
};

double lVol::getForward(double time) const {
  double firstMaturity = sabrModels_[0].getMaturity();
  if (time <= firstMaturity) {
    double afterForward = sabrModels_[0].getForward();
    return spot_ + (afterForward - spot_) * time /
                       firstMaturity; // linear interpolation
  }
  int N = sabrModels_.size();
  double lastMaturity = sabrModels_[N - 1].getMaturity();
  if (time > lastMaturity)
    return sabrModels_[N - 1].getForward(); // flat extrapolation
  for (int i = 1; i < N; i++) {
    double afterMaturity = sabrModels_[i].getMaturity();
    double beforeMaturity = sabrModels_[i - 1].getMaturity();
    if (time > beforeMaturity && time <= afterMaturity) {
      double beforeForward = sabrModels_[i - 1].getForward();
      double afterForward = sabrModels_[i].getForward();
      return beforeForward + (afterForward - beforeForward) *
                                 (time - beforeMaturity) /
                                 (afterMaturity - beforeMaturity);
    }
  }
  return 0.0;
};

double lVol::getLocVol(double preTime, double preSpot, double preFwd) const {
  double firstMaturity = sabrModels_[0].getMaturity();
  if (preTime < firstMaturity) {
    double cVal = interpolatedCall(sabrModels_[0], preSpot, preFwd, preTime);
    double cLft =
        interpolatedCall(sabrModels_[0], preSpot * 0.999, preFwd, preTime);
    double cRgt =
        interpolatedCall(sabrModels_[0], preSpot * 1.001, preFwd, preTime);
    double dS2val =
        (cRgt + cLft - 2.0 * cVal) / (0.001 * preSpot * 0.001 * preSpot);
    // assert(dS2val>0.0);
    dS2val = std::max(1.0E-20, dS2val);
    double cBigT =
        interpolatedCall(sabrModels_[0], preSpot, preFwd, preTime * 1.001);
    double dTval = (cBigT - cVal) / (0.001 * preTime);
    // assert(dTval>0.0);
    dTval = std::max(0.0, dTval); // testing
    return std::sqrt(2.0 * dTval / dS2val) / preSpot;
  }
  if (preTime == firstMaturity)
    return sabrModels_[0].localVol(preSpot);
  int N = sabrModels_.size();
  double lastMaturity = sabrModels_[N - 1].getMaturity();
  if (preTime > lastMaturity)
    return sabrModels_[N - 1].localVol(preSpot);
  for (int i = 1; i < N; i++) {
    double afterMaturity = sabrModels_[i].getMaturity();
    double beforeMaturity = sabrModels_[i - 1].getMaturity();
    if (preTime > beforeMaturity && preTime < afterMaturity) {
      double cVal = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                     preSpot, preFwd, preTime);
      double cLft = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                     preSpot * 0.999, preFwd, preTime);
      double cRgt = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                     preSpot * 1.001, preFwd, preTime);
      double dS2val =
          (cRgt + cLft - 2.0 * cVal) / (0.001 * preSpot * 0.001 * preSpot);
      // assert(dS2val>0.0);
      dS2val = std::max(1.0E-20, dS2val);
      double cBigT = interpolatedCall(sabrModels_[i - 1], sabrModels_[i],
                                      preSpot, preFwd, preTime * 1.001);
      double dTval = (cBigT - cVal) / (0.001 * preTime);
      // assert(dTval>0.0);
      dTval = std::max(0.0, dTval); // testing
      return std::sqrt(2.0 * dTval / dS2val) / preSpot;
    }
    if (preTime == afterMaturity)
      return sabrModels_[i].localVol(preSpot);
  }
  return 0.0;
};

std::vector<double> lVol::simulation(std::vector<double> timesPath,
                                     std::vector<double> randsPath) const {
  double preSpot = spot_;
  double preFwd = spot_;
  double preTime = 0.0;
  double driftPath = 0.0;
  int N = timesPath.size();
  std::vector<double> spotsPath(N);
  for (int i = 0; i < N; i++) {
    double currentTime = timesPath[i];
    double currentLocVol, currentForward;
    getLocVol(currentTime, preSpot, preFwd, currentLocVol, currentForward);
    // getLocVolOLD(currentTime,preSpot,preFwd,currentLocVol,currentForward);
    double currentVol = currentLocVol * std::sqrt(currentTime - preTime);
    double currentDrift = -0.5 * currentVol * currentVol +
                          currentVol * randsPath[i]; // Euler scheme
    driftPath = driftPath + currentDrift;
    double currentSpot = currentForward * std::exp(driftPath); // lognormal
    spotsPath[i] = currentSpot;
    preSpot = currentSpot;
    preFwd = currentForward;
    preTime = currentTime;
  }
  return spotsPath;
};

std::vector<double> lVol::simulation(std::vector<double> timesPath) const {
  double preSpot = spot_;
  double preFwd = spot_;
  double preTime = 0.0;
  int N = timesPath.size();
  std::vector<double> spotsPath(N);
  for (int i = 0; i < N; i++) {
    double currentTime = timesPath[i];
    if (preTime == 0.0)
      preTime = 0.75 * currentTime;
    // preTime = max(0.002,preTime);  // one day 1/365=0.002
    double preLocVol = getLocVol(preTime, preSpot, preFwd);
    double currentVol = preLocVol * std::sqrt(currentTime - preTime);
    double z = random_normal();
    double currentDrift =
        currentVol * z - 0.5 * currentVol * currentVol; // Euler scheme
    double bumpLocVol = getLocVol(preTime, 1.001 * preSpot, preFwd);
    double skewLocVol = (bumpLocVol - preLocVol) / (0.001 * preSpot);
    currentDrift = currentDrift + 0.5 * (z * z - 1.0) * preLocVol *
                                      (currentTime - preTime) *
                                      skewLocVol; // Milstein scheme
    double currentForward = getForward(currentTime);
    double currentSpot =
        preSpot * std::exp(currentDrift) * (currentForward / preFwd);
    currentSpot = currentForward *
                  std::min(std::max(0.01, currentSpot / currentForward), 50.0);
    spotsPath[i] = currentSpot;
    preSpot = currentSpot;
    preFwd = currentForward;
    preTime = currentTime;
  }
  return spotsPath;
};

// std::vector<double> lVol::simulation(std::vector<double> timesPath) const
//{
//	double preSpot = spot_;
//	double preFwd = spot_;
//	double preTime = 0.0;
//	double pathDrift = 0.0;
//	int N = timesPath.size();
//	std::vector<double> spotsPath(N);
//	for(int i=0; i<N; i++)
//	{
//		preTime = max(0.004,preTime);
//		double preLocVol = getLocVol(preTime,preSpot,preFwd);
//		double currentTime = timesPath[i];
//		double currentVol = preLocVol*sqrt(currentTime-preTime);
//		double z = random_normal();
//		double currentDrift = currentVol*z -0.5*currentVol*currentVol;
////Euler scheme 		double bumpLocVol =
/// getLocVol(preTime, 1.001*preSpot, preFwd);
//		double skewLocVol = (bumpLocVol-preLocVol)/(0.001*preSpot);
//		currentDrift = currentDrift
//+0.5*(z*z-1.0)*preLocVol*(currentTime-preTime)*skewLocVol; //Milstein scheme
//		pathDrift = pathDrift + currentDrift;
//		double currentForward = getForward(currentTime);
//		double currentSpot = currentForward * min(max(0.001,
// exp(pathDrift)), 500.0); 		spotsPath[i] = currentSpot;
// preSpot = currentSpot; 		preFwd = currentForward;
// preTime = currentTime;
//	}
//	return spotsPath;
// };

std::vector<double> lVol::simulationTEST(std::vector<double> timesPath) const {
  double preSpot = spot_;
  double preFwd = spot_;
  double preTime = 0.0;
  double driftPath = 0.0;
  int N = timesPath.size();
  std::vector<double> spotsPath(N);
  for (int i = 0; i < N; i++) {
    double currentTime = timesPath[i];
    double currentSpot, currentForward, currentLocVol;
    getLocVol(currentTime, preSpot, preFwd, currentLocVol, currentForward);
    double currentVol = currentLocVol * std::sqrt(currentTime - preTime);
    double currentDrift = -0.5 * currentVol * currentVol +
                          currentVol * random_normal(); // Euler scheme
    // double bumpLocVol, bumpForward;
    // getLocVol(currentTime, 1.001*preSpot, preFwd,bumpLocVol,bumpForward);
    // double skewLocVol = (bumpLocVol-currentLocVol)/(0.001*preSpot);
    ////if(skewLocVol < -9.9) skewLocVol = -9.9;	if(skewLocVol > 9.9)
    /// skewLocVol = 9.9;
    // double z = random_normal();
    // double currentDrift = -0.5*currentVol*currentVol + currentVol*z
    //						+0.5*(z*z-1.0)*currentLocVol*(currentTime-preTime)*skewLocVol;
    ////Milstein scheme
    driftPath = driftPath + currentDrift;
    currentSpot = currentForward * std::min(std::max(RFLO, std::exp(driftPath)),
                                            RTOP); // lognormal
    spotsPath[i] = currentSpot;
    preSpot = currentSpot;
    preFwd = currentForward;
    preTime = currentTime;
  }
  return spotsPath;
};

void lVol::getLocVolOLD(double time, double preSpot, double preFwd,
                        double &locvol, double &forward) const {
  double firstMaturity = sabrModels_[0].getMaturity();
  if (time <= firstMaturity) {
    double nextForward = sabrModels_[0].getForward();
    forward = spot_ + time * (nextForward - spot_) /
                          firstMaturity;              // linear interpolation
    double nextSpot = nextForward * preSpot / preFwd; // same ratio
    locvol = sabrModels_[0].localVol(nextSpot);
    return;
  }
  int N = sabrModels_.size();
  double lastMaturity = sabrModels_[N - 1].getMaturity();
  if (time > lastMaturity) {
    double nextForward = sabrModels_[N - 1].getForward();
    forward = nextForward;
    double nextSpot = nextForward * preSpot / preFwd;
    locvol = sabrModels_[N - 1].localVol(nextSpot);
    return;
  }
  for (int i = 1; i < N; i++) {
    double preMaturity = sabrModels_[i - 1].getMaturity();
    double nextMaturity = sabrModels_[i].getMaturity();
    if (time > preMaturity && time <= nextMaturity) {
      double preForward = sabrModels_[i - 1].getForward();
      double nextForward = sabrModels_[i].getForward();
      forward = preForward + (time - preMaturity) * (nextForward - preForward) /
                                 (nextMaturity - preMaturity);
      double preLocvol = sabrModels_[i - 1].localVolzabr(preSpot);
      double nextSpot = nextForward * preSpot / preFwd;
      double nextLocvol = sabrModels_[i].localVolzabr(nextSpot);
      locvol = preLocvol + (time - preMaturity) * (nextLocvol - preLocvol) /
                               (nextMaturity - preMaturity); // linear locvol
      return;
    }
  }
};
} // namespace velesquant