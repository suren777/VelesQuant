
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <ql/errors.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <velesquant/volatility/lm.h>
#include <velesquant/numerics/tri_diag_matrix.h>
#include <velesquant/pde_solvers/short_rate_2f_pde.h>
using namespace std;
#pragma warning(disable : 4715)
#pragma warning(disable : 4018)

namespace velesquant {

static const int Mxy = 101;
static const double Wxy[Mxy] = {
    // equidistant (uniform) discretisation
    -4.0000, -3.6480, -3.3269, -3.0341, -2.7670, -2.5234, -2.3013, -2.0986,
    -1.9138, -1.7453, -1.5916, -1.4513, -1.3234, -1.2067, -1.1003, -1.0032,
    -0.9147, -0.8339, -0.7601, -0.6929, -0.6315, -0.5754, -0.5243, -0.4776,
    -0.4350, -0.3960, -0.3604, -0.3279, -0.2982, -0.2710, -0.2461, -0.2232,
    -0.2023, -0.1831, -0.1654, -0.1492, -0.1342, -0.1204, -0.1076, -0.0956,
    -0.0846, -0.0742, -0.0644, -0.0552, -0.0465, -0.0382, -0.0301, -0.0224,
    -0.0148, -0.0074, 0.0000,  0.0074,  0.0148,  0.0224,  0.0301,  0.0382,
    0.0465,  0.0552,  0.0644,  0.0742,  0.0846,  0.0956,  0.1076,  0.1204,
    0.1342,  0.1492,  0.1654,  0.1831,  0.2023,  0.2232,  0.2461,  0.2710,
    0.2982,  0.3279,  0.3604,  0.3960,  0.4350,  0.4776,  0.5243,  0.5754,
    0.6315,  0.6929,  0.7601,  0.8339,  0.9147,  1.0032,  1.1003,  1.2067,
    1.3234,  1.4513,  1.5916,  1.7453,  1.9138,  2.0986,  2.3013,  2.5234,
    2.7670,  3.0341,  3.3269,  3.6480,  4.0000};

void ShortRate2FPDE::buildGrid(double Time, int Nt) {
  gridT_.resize(Nt);
  double delT = Time / (Nt - 1);
  gridT_[0] = 0.0;
  for (int t = 1; t < Nt - 1; t++)
    gridT_[t] = t * delT;
  gridT_[Nt - 1] = Time;
  gridX_.resize(Mxy);
  double avgSigma = 0.0;
  int nP = timeSigma1s_.size();
  for (int n = 0; n < nP; n++) {
    if (n == 0)
      avgSigma += sigma1s_[n] * sigma1s_[n] * timeSigma1s_[n];
    else
      avgSigma +=
          sigma1s_[n] * sigma1s_[n] * (timeSigma1s_[n] - timeSigma1s_[n - 1]);
  }
  avgSigma = sqrt(avgSigma / timeSigma1s_[nP - 1]) *
             sqrt((1 - exp(-2 * kappa1_ * Time)) / (2 * kappa1_));
  for (int x = 0; x < Mxy; x++)
    gridX_[x] = Wxy[x] * avgSigma;
  gridY_.resize(Mxy);
  avgSigma = 0.0;
  nP = timeSigma2s_.size();
  for (int n = 0; n < nP; n++) {
    if (n == 0)
      avgSigma += sigma2s_[n] * sigma2s_[n] * timeSigma2s_[n];
    else
      avgSigma +=
          sigma2s_[n] * sigma2s_[n] * (timeSigma2s_[n] - timeSigma2s_[n - 1]);
  }
  avgSigma = sqrt(avgSigma / timeSigma2s_[nP - 1]) *
             sqrt((1 - exp(-2 * kappa2_ * Time)) / (2 * kappa2_));
  for (int x = 0; x < Mxy; x++)
    gridY_[x] = Wxy[x] * avgSigma;
};

double ShortRate2FPDE::pricingZB(double Maturity) {
  int Nt = int(Maturity * 2520 + 0.5) / 10; // daily time grid?
  buildGrid(Maturity, Nt);
  vector<vector<double>> payoff;
  payoff.resize(Mxy);
  for (int x = 0; x < Mxy; x++) {
    payoff[x].resize(Mxy);
    for (int y = 0; y < Mxy; y++)
      payoff[x][y] = 1.0;
  }
  vector<vector<double>> pv(payoff);
  for (int t = Nt - 2; t >= 0; t--) {
    oneStepBackwardADIDouglas(t, payoff, pv);
    payoff = pv;
  }
  double zeroBond = 0.0;
  for (int x = 0; x < Mxy; x++)
    if (gridX_[x] >= 0.0) {
      for (int y = 0; y < Mxy; y++)
        if (gridY_[y] >= 0.0) {
          zeroBond = payoff[x][y];
          break;
        }
      break;
    }
  QL_ENSURE(zeroBond <= 1.001, "pricingZB Functor Fails " << zeroBond);
  return zeroBond;
};
double ShortRate2FPDE::pricingSwaption(double Expiry, double Tenor,
                                       double Strike, double PayFrequency) {
  int Nt = int((Expiry + Tenor) * 2520 + 0.5) / 10; // working daily time grid
  buildGrid((Expiry + Tenor), Nt);
  int iExpiry = int(Expiry * 2520 + 0.5) / 10;
  int Ncoupon = int(Tenor / PayFrequency + 0.5);
  vector<vector<double>> payoff;
  payoff.resize(Mxy);
  for (int x = 0; x < Mxy; x++) {
    payoff[x].resize(Mxy);
    for (int y = 0; y < Mxy; y++)
      payoff[x][y] = -1.0;
  }
  vector<vector<double>> pv(payoff);
  for (int t = Nt - 2; t >= iExpiry; t--) // swap
  {
    double couponTime = Expiry + Ncoupon * PayFrequency;
    if (couponTime > gridT_[t] && couponTime <= gridT_[t + 1]) {
      for (int x = 0; x < Mxy; x++)
        for (int y = 0; y < Mxy; y++)
          payoff[x][y] -= PayFrequency * Strike;
      Ncoupon--;
    }
    oneStepBackwardADIDouglas(t, payoff, pv);
    payoff = pv;
  }
  for (int x = 0; x < Mxy; x++)
    for (int y = 0; y < Mxy; y++)
      payoff[x][y] = max(0.0, payoff[x][y] + 1.0); // swaption payoff
  for (int t = iExpiry - 1; t >= 0; t--) {
    oneStepBackwardADIDouglas(t, payoff, pv);
    payoff = pv;
  }
  double swaptionValue = 0.0;
  for (int x = 0; x < Mxy; x++)
    if (gridX_[x] >= 0.0) {
      for (int y = 0; y < Mxy; y++)
        if (gridY_[y] >= 0.0) {
          swaptionValue = payoff[x][y];
          break;
        }
      break;
    }
  QL_ENSURE(swaptionValue >= -0.0001,
            "pricingSwaption Functor Fails " << swaptionValue);
  return swaptionValue;
};

void ShortRate2FPDE::oneStepBackwardADIDouglas(
    int t, const vector<vector<double>> &inM, vector<vector<double>> &outM) {
  oneStepBackwardExplicit(t, inM, outM);
  vector<vector<double>> midM(outM);
  oneStepBackwardDouglasX(t, inM, midM, outM);
  midM = outM;
  oneStepBackwardDouglasY(t, inM, midM, outM);
};
void ShortRate2FPDE::oneStepBackwardExplicit(int t,
                                             const vector<vector<double>> &inM,
                                             vector<vector<double>> &outM) {
  double te = gridT_[t + 1];
  double ti = gridT_[t];
  double delT = te - ti;
  double tm = 0.5 * (te + ti);
  double sigma1 = whichValue(tm, timeSigma1s_, sigma1s_);
  double difu1 = sigma1 * sigma1;
  double sigma2 = whichValue(tm, timeSigma2s_, sigma2s_);
  double difu2 = sigma2 * sigma2;
  double alpha = whichValue(tm, timeAlphas_, alphas_);

  // j=0
  double X0 = gridX_[0];
  double X1 = gridX_[1];
  double X2 = gridX_[2];
  double conv1 = kappa1_ * X0;
  // j=0, k=0 lower boundary condition
  double Y0 = gridY_[0];
  double Y1 = gridY_[1];
  double Y2 = gridY_[2];
  double conv2 = lambda_ * X0 + kappa2_ * Y0;
  double rate = alpha + X0 + Y0;
  outM[0][0] =
      (1.0 / delT - rate) * inM[0][0] +
      (conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) * inM[0][0] +
      (-conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][0] +
      (conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][0] +
      (conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) * inM[0][0] +
      (-conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[0][1] +
      (conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[0][2];
  outM[0][0] = delT * outM[0][0];
  // j=0, k\=0 middle range
  double Y, Yu, Yl;
  for (int k = 1; k < Mxy - 1; k++) {
    Y = gridY_[k];
    Yu = gridY_[k + 1];
    Yl = gridY_[k - 1];
    conv2 = lambda_ * X0 + kappa2_ * Y;
    rate = alpha + X0 + Y;
    outM[0][k] =
        (1.0 / delT - rate) * inM[0][k] +
        (conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][k] +
        (-conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
        (conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
        (conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[0][k - 1] +
        (-conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) * inM[0][k] +
        (-conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[0][k + 1];
    outM[0][k] = delT * outM[0][k];
  }
  // j=0, k=M-1 upper boundary condition
  double Y3l = gridY_[Mxy - 3];
  double Y2l = gridY_[Mxy - 2];
  double Y1l = gridY_[Mxy - 1];
  conv2 = lambda_ * X0 + kappa2_ * Y1l;
  rate = alpha + X0 + Y1l;
  outM[0][Mxy - 1] =
      (1.0 / delT - rate) * inM[0][Mxy - 1] +
      (conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
          inM[0][Mxy - 1] +
      (-conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][Mxy - 1] +
      (conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][Mxy - 1] +
      (-conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
          inM[0][Mxy - 3] +
      (conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
          inM[0][Mxy - 2] +
      (-conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
          inM[0][Mxy - 1];
  outM[0][Mxy - 1] = delT * outM[0][Mxy - 1];

  // j\=0
  double X, Xu, Xl;
  for (int j = 1; j < Mxy - 1; j++) {
    X = gridX_[j];
    Xu = gridX_[j + 1];
    Xl = gridX_[j - 1];
    conv1 = kappa1_ * X;
    // j\=0, k=0 lower boundary condition
    conv2 = lambda_ * X + kappa2_ * Y0;
    rate = alpha + X + Y0;
    outM[j][0] =
        (1.0 / delT - rate) * inM[j][0] +
        (conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][0] +
        (-conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) * inM[j][0] +
        (-conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][0] +
        (conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[j][0] +
        (-conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
        (conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2];
    outM[j][0] = delT * outM[j][0];
    // j\=0, k\=0 middle range
    for (int k = 1; k < Mxy - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = lambda_ * X + kappa2_ * Y;
      rate = alpha + X + Y;
      outM[j][k] =
          (1.0 / delT - rate) * inM[j][k] +
          (conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
          (-conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
              inM[j][k] +
          (-conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
          (conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
          (-conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
              inM[j][k] +
          (-conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1];
      outM[j][k] = delT * outM[j][k];
    }
    // j\=0, k=M-1 upper boundary condition
    conv2 = lambda_ * X + kappa2_ * Y1l;
    rate = alpha + X + Y1l;
    outM[j][Mxy - 1] = (1.0 / delT - rate) * inM[j][Mxy - 1] +
                       (conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) *
                           inM[j - 1][Mxy - 1] +
                       (-conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) /
                           (Xu - X) * inM[j][Mxy - 1] +
                       (-conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) *
                           inM[j + 1][Mxy - 1] +
                       (-conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) /
                           (Y1l - Y3l) * inM[j][Mxy - 3] +
                       (conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) /
                           (Y1l - Y2l) * inM[j][Mxy - 2] +
                       (-conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) /
                           (Y1l - Y3l) * inM[j][Mxy - 1];
    outM[j][Mxy - 1] = delT * outM[j][Mxy - 1];
  }

  // j=M-1
  double X3l = gridX_[Mxy - 3];
  double X2l = gridX_[Mxy - 2];
  double X1l = gridX_[Mxy - 1];
  conv1 = kappa1_ * X1l;
  // j=M-1, k=0 lower boundary condition
  conv2 = lambda_ * X1l + kappa2_ * Y0;
  rate = alpha + X1l + Y0;
  outM[Mxy - 1][0] =
      (1.0 / delT - rate) * inM[Mxy - 1][0] +
      (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
          inM[Mxy - 3][0] +
      (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
          inM[Mxy - 2][0] +
      (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
          inM[Mxy - 1][0] +
      (conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
          inM[Mxy - 1][0] +
      (-conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[Mxy - 1][1] +
      (conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[Mxy - 1][2];
  outM[Mxy - 1][0] = delT * outM[Mxy - 1][0];
  // j=M-1, k\=0 middle range
  for (int k = 1; k < Mxy - 1; k++) {
    Y = gridY_[k];
    Yu = gridY_[k + 1];
    Yl = gridY_[k - 1];
    conv2 = lambda_ * X1l + kappa2_ * Y;
    rate = alpha + X1l + Y;
    outM[Mxy - 1][k] = (1.0 / delT - rate) * inM[Mxy - 1][k] +
                       (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) /
                           (X1l - X3l) * inM[Mxy - 3][k] +
                       (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) /
                           (X1l - X2l) * inM[Mxy - 2][k] +
                       (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) /
                           (X1l - X3l) * inM[Mxy - 1][k] +
                       (conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) *
                           inM[Mxy - 1][k - 1] +
                       (-conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) /
                           (Yu - Y) * inM[Mxy - 1][k] +
                       (-conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) *
                           inM[Mxy - 1][k + 1];
    outM[Mxy - 1][k] = delT * outM[Mxy - 1][k];
  }
  // j=M-1, k=M-1 upper boundary condition
  conv2 = lambda_ * X1l + kappa2_ * Y1l;
  rate = alpha + X1l + Y1l;
  outM[Mxy - 1][Mxy - 1] =
      (1.0 / delT - rate) * inM[Mxy - 1][Mxy - 1] +
      (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
          inM[Mxy - 3][Mxy - 1] +
      (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
          inM[Mxy - 2][Mxy - 1] +
      (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
          inM[Mxy - 1][Mxy - 1] +
      (-conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
          inM[Mxy - 1][Mxy - 3] +
      (conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
          inM[Mxy - 1][Mxy - 2] +
      (-conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
          inM[Mxy - 1][Mxy - 1];
  outM[Mxy - 1][Mxy - 1] = delT * outM[Mxy - 1][Mxy - 1];
};
void ShortRate2FPDE::oneStepBackwardDouglasX(int t,
                                             const vector<vector<double>> &inM,
                                             const vector<vector<double>> &midM,
                                             vector<vector<double>> &outM) {
  vector<double> l(Mxy - 2), c(Mxy - 2), u(Mxy - 2), d(Mxy - 2), V(Mxy - 2);
  double te = gridT_[t + 1];
  double ti = gridT_[t];
  double delT = te - ti;
  double tm = 0.5 * (te + ti);
  double sigma1 = whichValue(tm, timeSigma1s_, sigma1s_);
  double difu1 = 0.5 * sigma1 * sigma1;
  double alpha = whichValue(tm, timeAlphas_, alphas_);
  double conv1, rate;

  // for all y (k=0,...,M-1)
  for (int k = 0; k < Mxy; k++) {
    double Y = gridY_[k];
    // x in the middle range
    for (int j = 1; j < Mxy - 1; j++) {
      double X = gridX_[j];
      double Xu = gridX_[j + 1];
      double Xl = gridX_[j - 1];
      conv1 = 0.5 * kappa1_ * X;
      rate = 0.25 * (alpha + X + Y);
      l[j - 1] = (-conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl);
      c[j - 1] = 1.0 / delT + rate +
                 (conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X);
      u[j - 1] = (conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl);
      d[j - 1] =
          (-conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
          (rate + (conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X)) *
              inM[j][k] +
          (conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
          midM[j][k] / delT;
    }
    // j=0 lower boundary condition
    double X0 = gridX_[0];
    double X1 = gridX_[1];
    double X2 = gridX_[2];
    conv1 = 0.5 * kappa1_ * X0;
    rate = 0.25 * (alpha + X0 + Y);
    double k0 = 1.0 / delT + rate +
                (-conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0);
    double k1 = (conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1);
    double k2 = (-conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0);
    double d0 =
        (rate +
         (-conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0)) *
            inM[0][k] +
        (conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
        (-conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
        midM[0][k] / delT;
    c[0] = c[0] - l[0] * k1 / k0;
    u[0] = u[0] - l[0] * k2 / k0;
    d[0] = d[0] - l[0] * d0 / k0;
    l[0] = 0.0;
    // j=M-1 upper boundary condition
    double X3 = gridX_[Mxy - 3];
    X2 = gridX_[Mxy - 2];
    X1 = gridX_[Mxy - 1];
    conv1 = 0.5 * kappa1_ * X1;
    rate = 0.25 * (alpha + X1 + Y);
    double kl3 = (conv1 * (X1 - X2) - difu1) / (X2 - X3) / (X1 - X3);
    double kl2 = (-conv1 * (X1 - X3) + difu1) / (X2 - X3) / (X1 - X2);
    double kl1 = 1.0 / delT + rate +
                 (conv1 * (2.0 * X1 - X2 - X3) - difu1) / (X1 - X2) / (X1 - X3);
    double dl1 =
        (conv1 * (X1 - X2) - difu1) / (X2 - X3) / (X1 - X3) * inM[Mxy - 3][k] +
        (-conv1 * (X1 - X3) + difu1) / (X2 - X3) / (X1 - X2) * inM[Mxy - 2][k] +
        (rate +
         (conv1 * (2.0 * X1 - X2 - X3) - difu1) / (X1 - X2) / (X1 - X3)) *
            inM[Mxy - 1][k] +
        midM[Mxy - 1][k] / delT;
    l[Mxy - 3] = l[Mxy - 3] - u[Mxy - 3] * kl3 / kl1;
    c[Mxy - 3] = c[Mxy - 3] - u[Mxy - 3] * kl2 / kl1;
    d[Mxy - 3] = d[Mxy - 3] - u[Mxy - 3] * dl1 / kl1;
    u[Mxy - 3] = 0.0;
    TriDiagonalSolve(Mxy - 2, l, c, u, d, V);
    for (int j = 1; j < Mxy - 1; j++)
      outM[j][k] = V[j - 1];
    // update j=0 lower boundary
    outM[0][k] = (d0 - k1 * outM[1][k] - k2 * outM[2][k]) / k0;
    // update j=M-1 upper boundary
    outM[Mxy - 1][k] =
        (dl1 - kl3 * outM[Mxy - 3][k] - kl2 * outM[Mxy - 2][k]) / kl1;
  }
};
void ShortRate2FPDE::oneStepBackwardDouglasY(int t,
                                             const vector<vector<double>> &inM,
                                             const vector<vector<double>> &midM,
                                             vector<vector<double>> &outM) {
  vector<double> l(Mxy - 2), c(Mxy - 2), u(Mxy - 2), d(Mxy - 2), V(Mxy - 2);
  double te = gridT_[t + 1];
  double ti = gridT_[t];
  double delT = te - ti;
  double tm = 0.5 * (te + ti);
  double sigma2 = whichValue(tm, timeSigma2s_, sigma2s_);
  double difu2 = 0.5 * sigma2 * sigma2;
  double alpha = whichValue(tm, timeAlphas_, alphas_);
  double conv2, rate;

  // for all x (j=0,...,M-1)
  for (int j = 0; j < Mxy; j++) {
    double X = gridX_[j];
    // y in the middle range
    for (int k = 1; k < Mxy - 1; k++) {
      double Y = gridY_[k];
      double Yu = gridY_[k + 1];
      double Yl = gridY_[k - 1];
      conv2 = 0.5 * (lambda_ * X + kappa2_ * Y);
      rate = 0.25 * (alpha + X + Y);
      l[k - 1] = (-conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl);
      c[k - 1] = 1.0 / delT + rate +
                 (conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y);
      u[k - 1] = (conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl);
      d[k - 1] =
          (-conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
          (rate + (conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y)) *
              inM[j][k] +
          (conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1] +
          midM[j][k] / delT;
    }
    // k=0 lower boundary condition
    double Y0 = gridY_[0];
    double Y1 = gridY_[1];
    double Y2 = gridY_[2];
    conv2 = 0.5 * (lambda_ * X + kappa2_ * Y0);
    rate = 0.25 * (alpha + X + Y0);
    double k0 = 1.0 / delT + rate +
                (-conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0);
    double k1 = (conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1);
    double k2 = (-conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0);
    double d0 =
        (rate +
         (-conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0)) *
            inM[j][0] +
        (conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
        (-conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2] +
        midM[j][0] / delT;
    c[0] = c[0] - l[0] * k1 / k0;
    u[0] = u[0] - l[0] * k2 / k0;
    d[0] = d[0] - l[0] * d0 / k0;
    l[0] = 0.0;
    // k=M-1 upper boundary condition
    double Y3 = gridY_[Mxy - 3];
    Y2 = gridY_[Mxy - 2];
    Y1 = gridY_[Mxy - 1];
    conv2 = 0.5 * (lambda_ * X + kappa2_ * Y1);
    rate = 0.25 * (alpha + X + Y1);
    double kl3 = (conv2 * (Y1 - Y2) - difu2) / (Y2 - Y3) / (Y1 - Y3);
    double kl2 = (-conv2 * (Y1 - Y3) + difu2) / (Y2 - Y3) / (Y1 - Y2);
    double kl1 = 1.0 / delT + rate +
                 (conv2 * (2.0 * Y1 - Y2 - Y3) - difu2) / (Y1 - Y2) / (Y1 - Y3);
    double dl1 =
        (conv2 * (Y1 - Y2) - difu2) / (Y2 - Y3) / (Y1 - Y3) * inM[j][Mxy - 3] +
        (-conv2 * (Y1 - Y3) + difu2) / (Y2 - Y3) / (Y1 - Y2) * inM[j][Mxy - 2] +
        (rate +
         (conv2 * (2.0 * Y1 - Y2 - Y3) - difu2) / (Y1 - Y2) / (Y1 - Y3)) *
            inM[j][Mxy - 1] +
        midM[j][Mxy - 1] / delT;
    l[Mxy - 3] = l[Mxy - 3] - u[Mxy - 3] * kl3 / kl1;
    c[Mxy - 3] = c[Mxy - 3] - u[Mxy - 3] * kl2 / kl1;
    d[Mxy - 3] = d[Mxy - 3] - u[Mxy - 3] * dl1 / kl1;
    u[Mxy - 3] = 0.0;
    TriDiagonalSolve(Mxy - 2, l, c, u, d, V);
    for (int k = 1; k < Mxy - 1; k++)
      outM[j][k] = V[k - 1];
    // update k=0 lower boundary
    outM[j][0] = (d0 - k1 * outM[j][1] - k2 * outM[j][2]) / k0;
    // update k=M-1 upper boundary
    outM[j][Mxy - 1] =
        (dl1 - kl3 * outM[j][Mxy - 3] - kl2 * outM[j][Mxy - 2]) / kl1;
  }
};

void ShortRate2FPDE::calibrator(vector<double> timeDFs, vector<double> DFs,
                                vector<defSwap> swapQuotes) {
  DFs_ = DFs;
  quoteSwap_ = swapQuotes;
  int nD = timeDFs.size();
  timeAlphas_.resize(nD);
  timeAlphas_ = timeDFs;
  alphas_.resize(nD, 0.015);
  int ns1 = sigma1s_.size();
  int ns2 = sigma2s_.size();
  int n = ns1 + ns2 + 3;     // no. of model volatility term paremeters
  double *x = new double[n]; // initial estimate of parameters vector
  for (int i = 0; i < ns1; i++)
    x[i] = sigma1s_[i]; // sigma1s initial value
  for (int i = 0; i < ns2; i++)
    x[ns1 + i] = sigma2s_[i]; // sigma2s initial value
  x[n - 1] = kappa1_;         // kappa1 initial value
  x[n - 2] = kappa2_;         // kappa2 initial value
  x[n - 3] = lambda_;         // lambda initial value
  int m = quoteSwap_.size();  // no. of observations
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
  boost::function<void(ShortRate2FPDE *, int, int, double *, double *, int *)>
      obj = &ShortRate2FPDE::objFcnCalibrator;
  lmfcn fcn = boost::bind(obj, this, _1, _2, _3, _4, _5);
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
  QL_ENSURE(info != 4, "Model Calibration Fails " << info);
  // the below is output result
  for (int i = 0; i < ns1; i++)
    sigma1s_[i] = fabs(x[i]); // sigmas1 final value
  for (int i = 0; i < ns2; i++)
    sigma2s_[i] = fabs(x[ns1 + i]); // sigmas2 final value
  kappa1_ = x[n - 1];               // kappa1 final value
  kappa2_ = x[n - 2];               // kappa2 final value
  lambda_ = x[n - 3];               // lambda final value
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
void ShortRate2FPDE::objFcnCalibrator(int m, int n, double *x, double *fvec,
                                      int *iflag) {
  int ns1 = sigma1s_.size();
  int ns2 = sigma2s_.size();
  for (int i = 0; i < ns1; i++)
    sigma1s_[i] = fabs(x[i]); // sigmas1 final value
  for (int i = 0; i < ns2; i++)
    sigma2s_[i] = fabs(x[ns1 + i]); // sigmas2 final value
  kappa1_ = x[n - 1];               // kappa1 final value
  kappa2_ = x[n - 2];               // kappa2 final value
  lambda_ = x[n - 3];               // lambda final value
  termStructureCalibrator();        // TEREM STRUCTURE CALIBRATION
  for (int i = 0; i < m; i++) {
    double model =
        pricingSwaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                        quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
    double market = quoteSwap_[i].Value;
    fvec[i] = model - market;
  }
};
void ShortRate2FPDE::termStructureCalibrator() {
  int n = timeAlphas_.size();
  vector<int> iTs(n);
  for (int i = 0; i < n; i++)
    iTs[i] = int(timeAlphas_[i] * 2520 + 0.5) / 10;
  buildGrid(timeAlphas_[n - 1], iTs[n - 1]);
  vector<vector<double>> inV, outV, lastV;
  inV.resize(Mxy);
  for (int x = 0; x < Mxy; x++) {
    inV[x].resize(Mxy);
    for (int y = 0; y < Mxy; y++)
      inV[x][y] = 0.0;
  }
  outV = inV;
  for (int x = 1; x < Mxy - 1; x++)
    if (gridX_[x] >= 0.0) {
      for (int y = 1; y < Mxy - 1; y++)
        if (gridY_[y] >= 0.0) {
          inV[x][y] =
              4.0 / ((gridX_[x + 1] - gridX_[x - 1]) *
                     (gridY_[y + 1] -
                      gridY_[y - 1])); // Dirac delta as initial distribution
          break;
        }
      break;
    }
  lastV = inV; // remember the starting density
  for (int i = 0; i < n; i++) {
    int sT = 0;
    if (i > 0)
      sT = iTs[i - 1] - 1;
    double df = 0.0;
    int niter = 0;
    do {
      niter++;
      inV = lastV; // reset to the starting density
      if (alphas_[i] == 0.0)
        alphas_[i] = 0.0005; // case of 0 value handler
      alphas_[i] *= 1.001;   // 0.1% up for theta parameter
      for (int t = sT; t < iTs[i] - 1; t++) {
        oneStepForwardADIDouglas(t, inV, outV);
        inV = outV;
      }
      double dfUP = trapezoidal2D(inV);
      inV = lastV;         // reset to the starting density
      alphas_[i] /= 1.001; // back to theta
      for (int t = sT; t < iTs[i] - 1; t++) {
        oneStepForwardADIDouglas(t, inV, outV);
        inV = outV;
      }
      df = trapezoidal2D(inV);
      alphas_[i] *=
          (1 + 0.001 * (DFs_[i] - df) / (dfUP - df)); // Newton iteration
      alphas_[i] =
          max(-0.0050, alphas_[i]); // Cutoff the nagative value (BETTER)
    } while (fabs(1.0 - df / DFs_[i]) >= 1.0E-6 && niter < 7);
    inV = lastV; // reset to the starting density
    for (int t = sT; t < iTs[i] - 1; t++) {
      oneStepForwardADIDouglas(t, inV, outV);
      inV = outV;
    }
    lastV = inV; // remember the starting density
  };
  return;
};

vector<double> ShortRate2FPDE::calculateDFs(vector<double> &timePoints) {
  int n = timePoints.size();
  vector<int> iTs(n);
  for (int i = 0; i < n; i++)
    iTs[i] = int(timePoints[i] * 2520 + 0.5) / 10;
  buildGrid(timePoints[n - 1], iTs[n - 1]);
  vector<double> DFs(n, 0.0);
  vector<vector<double>> inV;
  inV.resize(Mxy);
  for (int x = 0; x < Mxy; x++) {
    inV[x].resize(Mxy);
    for (int y = 0; y < Mxy; y++)
      inV[x][y] = 0.0;
  }
  for (int x = 1; x < Mxy - 1; x++)
    if (gridX_[x] >= 0.0) {
      for (int y = 1; y < Mxy - 1; y++)
        if (gridY_[y] >= 0.0) {
          // inV[x][y]
          // = 4.0/((gridX_[x+1]-gridX_[x])*(gridY_[y+1]-gridY_[y])+(gridX_[x+1]-gridX_[x])*(gridY_[y]-gridY_[y-1])
          //	            +(gridX_[x]-gridX_[x-1])*(gridY_[y+1]-gridY_[y])+(gridX_[x]-gridX_[x-1])*(gridY_[y]-gridY_[y-1]));
          //// Dirac delta as initial distribution
          inV[x][y] =
              4.0 / ((gridX_[x + 1] - gridX_[x - 1]) *
                     (gridY_[y + 1] -
                      gridY_[y - 1])); // Dirac delta as initial distribution
          break;
        }
      break;
    }
  vector<vector<double>> outV(inV);
  for (int i = 0; i < n; i++) {
    int sT = 0;
    if (i > 0)
      sT = iTs[i - 1] - 1;
    for (int t = sT; t < iTs[i] - 1; t++) {
      oneStepForwardADIDouglas(t, inV, outV);
      inV = outV;
    }
    DFs[i] = trapezoidal2D(inV);
    QL_ENSURE(DFs[i] <= 1.001, "pricingDF Functor Fails " << DFs[i]);
  }
  return DFs;
};

double ShortRate2FPDE::trapezoidal2D(vector<vector<double>> &inV) {
  double value = 0.0;
  for (int x = 1; x < Mxy; x++)
    for (int y = 1; y < Mxy; y++)
      value += 0.25 *
               (inV[x][y] + inV[x][y - 1] + inV[x - 1][y] + inV[x - 1][y - 1]) *
               (gridX_[x] - gridX_[x - 1]) * (gridY_[y] - gridY_[y - 1]);
  return value;
};

void ShortRate2FPDE::oneStepForwardADIDouglas(int T,
                                              const vector<vector<double>> &inM,
                                              vector<vector<double>> &outM) {
  oneStepForwardExplicit(T, inM, outM);
  vector<vector<double>> midM(outM);
  oneStepForwardDouglasX(T, inM, midM, outM);
  midM = outM;
  oneStepForwardDouglasY(T, inM, midM, outM);
};
void ShortRate2FPDE::oneStepForwardExplicit(int T,
                                            const vector<vector<double>> &inM,
                                            vector<vector<double>> &outM) {
  double Te = gridT_[T + 1];
  double Ti = gridT_[T];
  double delT = Te - Ti;
  double Tm = 0.5 * (Te + Ti);
  double sigma1 = whichValue(Tm, timeSigma1s_, sigma1s_);
  double difu1 = sigma1 * sigma1;
  double sigma2 = whichValue(Tm, timeSigma2s_, sigma2s_);
  double difu2 = sigma2 * sigma2;
  double alpha = whichValue(Tm, timeAlphas_, alphas_);

  // j=0
  double X0 = gridX_[0];
  double X1 = gridX_[1];
  double X2 = gridX_[2];
  double conv1 = kappa1_ * X0;
  // j=0, k=0 lower boundary condition
  double Y0 = gridY_[0];
  double Y1 = gridY_[1];
  double Y2 = gridY_[2];
  double conv2 = lambda_ * X0 + kappa2_ * Y0;
  double rate = alpha + X0 + Y0 - kappa1_ - kappa2_;
  outM[0][0] =
      (1.0 / delT - rate) * inM[0][0] +
      (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
          inM[0][0] +
      (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][0] +
      (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][0] +
      (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
          inM[0][0] +
      (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[0][1] +
      (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[0][2];
  outM[0][0] = delT * outM[0][0];
  // j=0, k\=0 middle range
  double Y, Yu, Yl;
  for (int k = 1; k < Mxy - 1; k++) {
    Y = gridY_[k];
    Yu = gridY_[k + 1];
    Yl = gridY_[k - 1];
    conv2 = lambda_ * X0 + kappa2_ * Y;
    rate = alpha + X0 + Y - kappa1_ - kappa2_;
    outM[0][k] =
        (1.0 / delT - rate) * inM[0][k] +
        (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][k] +
        (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
        (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
        (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[0][k - 1] +
        (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) * inM[0][k] +
        (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[0][k + 1];
    outM[0][k] = delT * outM[0][k];
  }
  // j=0, k=M-1 upper boundary condition
  double Y3l = gridY_[Mxy - 3];
  double Y2l = gridY_[Mxy - 2];
  double Y1l = gridY_[Mxy - 1];
  conv2 = lambda_ * X0 + kappa2_ * Y1l;
  rate = alpha + X0 + Y1l - kappa1_ - kappa2_;
  outM[0][Mxy - 1] =
      (1.0 / delT - rate) * inM[0][Mxy - 1] +
      (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
          inM[0][Mxy - 1] +
      (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][Mxy - 1] +
      (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][Mxy - 1] +
      (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
          inM[0][Mxy - 3] +
      (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
          inM[0][Mxy - 2] +
      (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
          inM[0][Mxy - 1];
  outM[0][Mxy - 1] = delT * outM[0][Mxy - 1];

  // j\=0
  double X, Xu, Xl;
  for (int j = 1; j < Mxy - 1; j++) {
    X = gridX_[j];
    Xu = gridX_[j + 1];
    Xl = gridX_[j - 1];
    conv1 = kappa1_ * X;
    // j\=0, k=0 lower boundary condition
    conv2 = lambda_ * X + kappa2_ * Y0;
    rate = alpha + X + Y0 - kappa1_ - kappa2_;
    outM[j][0] =
        (1.0 / delT - rate) * inM[j][0] +
        (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][0] +
        (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) * inM[j][0] +
        (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][0] +
        (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[j][0] +
        (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
        (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2];
    outM[j][0] = delT * outM[j][0];
    // j\=0, k\=0 middle range
    for (int k = 1; k < Mxy - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = lambda_ * X + kappa2_ * Y;
      rate = alpha + X + Y - kappa1_ - kappa2_;
      outM[j][k] =
          (1.0 / delT - rate) * inM[j][k] +
          (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
          (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
              inM[j][k] +
          (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
          (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
          (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
              inM[j][k] +
          (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1];
      outM[j][k] = delT * outM[j][k];
    }
    // j\=0, k=M-1 upper boundary condition
    conv2 = lambda_ * X + kappa2_ * Y1l;
    rate = alpha + X + Y1l - kappa1_ - kappa2_;
    outM[j][Mxy - 1] = (1.0 / delT - rate) * inM[j][Mxy - 1] +
                       (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) *
                           inM[j - 1][Mxy - 1] +
                       (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) /
                           (Xu - X) * inM[j][Mxy - 1] +
                       (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) *
                           inM[j + 1][Mxy - 1] +
                       (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) /
                           (Y1l - Y3l) * inM[j][Mxy - 3] +
                       (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) /
                           (Y1l - Y2l) * inM[j][Mxy - 2] +
                       (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) /
                           (Y1l - Y3l) * inM[j][Mxy - 1];
    outM[j][Mxy - 1] = delT * outM[j][Mxy - 1];
  }

  // j=M-1
  double X3l = gridX_[Mxy - 3];
  double X2l = gridX_[Mxy - 2];
  double X1l = gridX_[Mxy - 1];
  conv1 = kappa1_ * X1l;
  // j=M-1, k=0 lower boundary condition
  conv2 = lambda_ * X1l + kappa2_ * Y0;
  rate = alpha + X1l + Y0 - kappa1_ - kappa2_;
  outM[Mxy - 1][0] =
      (1.0 / delT - rate) * inM[Mxy - 1][0] +
      (conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
          inM[Mxy - 3][0] +
      (-conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
          inM[Mxy - 2][0] +
      (conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
          inM[Mxy - 1][0] +
      (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
          inM[Mxy - 1][0] +
      (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[Mxy - 1][1] +
      (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[Mxy - 1][2];
  outM[Mxy - 1][0] = delT * outM[Mxy - 1][0];
  // j=M-1, k\=0 middle range
  for (int k = 1; k < Mxy - 1; k++) {
    Y = gridY_[k];
    Yu = gridY_[k + 1];
    Yl = gridY_[k - 1];
    conv2 = lambda_ * X1l + kappa2_ * Y;
    rate = alpha + X1l + Y - kappa1_ - kappa2_;
    outM[Mxy - 1][k] =
        (1.0 / delT - rate) * inM[Mxy - 1][k] +
        (conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
            inM[Mxy - 3][k] +
        (-conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
            inM[Mxy - 2][k] +
        (conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
            inM[Mxy - 1][k] +
        (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) *
            inM[Mxy - 1][k - 1] +
        (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
            inM[Mxy - 1][k] +
        (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[Mxy - 1][k + 1];
    outM[Mxy - 1][k] = delT * outM[Mxy - 1][k];
  }
  // j=M-1, k=M-1 upper boundary condition
  conv2 = lambda_ * X1l + kappa2_ * Y1l;
  rate = alpha + X1l + Y1l - kappa1_ - kappa2_;
  outM[Mxy - 1][Mxy - 1] =
      (1.0 / delT - rate) * inM[Mxy - 1][Mxy - 1] +
      (conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
          inM[Mxy - 3][Mxy - 1] +
      (-conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
          inM[Mxy - 2][Mxy - 1] +
      (conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
          inM[Mxy - 1][Mxy - 1] +
      (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
          inM[Mxy - 1][Mxy - 3] +
      (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
          inM[Mxy - 1][Mxy - 2] +
      (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
          inM[Mxy - 1][Mxy - 1];
  outM[Mxy - 1][Mxy - 1] = delT * outM[Mxy - 1][Mxy - 1];
};
void ShortRate2FPDE::oneStepForwardDouglasX(int T,
                                            const vector<vector<double>> &inM,
                                            const vector<vector<double>> &midM,
                                            vector<vector<double>> &outM) {
  vector<double> l(Mxy - 2), c(Mxy - 2), u(Mxy - 2), d(Mxy - 2), V(Mxy - 2);
  double Te = gridT_[T + 1];
  double Ti = gridT_[T];
  double delT = Te - Ti;
  double Tm = 0.5 * (Te + Ti);
  double sigma1 = whichValue(Tm, timeSigma1s_, sigma1s_);
  double difu1 = 0.5 * sigma1 * sigma1;
  double alpha = whichValue(Tm, timeAlphas_, alphas_);
  double conv1, rate;

  // for all y (k=0,...,M-1)
  for (int k = 0; k < Mxy; k++) {
    double Y = gridY_[k];
    // x in the middle range
    for (int j = 1; j < Mxy - 1; j++) {
      double X = gridX_[j];
      double Xu = gridX_[j + 1];
      double Xl = gridX_[j - 1];
      conv1 = 0.5 * kappa1_ * X;
      rate = 0.25 * (alpha + X + Y - kappa1_ - kappa2_);
      l[j - 1] = (conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl);
      c[j - 1] = 1.0 / delT + rate +
                 (-conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X);
      u[j - 1] = (-conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl);
      d[j - 1] =
          (conv1 * (Xu - X) - difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
          (rate +
           (-conv1 * (Xu - 2.0 * X + Xl) + difu1) / (X - Xl) / (Xu - X)) *
              inM[j][k] +
          (-conv1 * (X - Xl) - difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
          midM[j][k] / delT;
    }
    // j=0 lower boundary condition
    double X0 = gridX_[0];
    double X1 = gridX_[1];
    double X2 = gridX_[2];
    conv1 = 0.5 * kappa1_ * X0;
    rate = 0.25 * (alpha + X0 + Y - kappa1_ - kappa2_);
    double k0 = 1.0 / delT + rate +
                (conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0);
    double k1 = (-conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1);
    double k2 = (conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0);
    double d0 =
        (rate +
         (conv1 * (X2 + X1 - 2.0 * X0) - difu1) / (X1 - X0) / (X2 - X0)) *
            inM[0][k] +
        (-conv1 * (X2 - X0) + difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
        (conv1 * (X1 - X0) - difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
        midM[0][k] / delT;
    c[0] = c[0] - l[0] * k1 / k0;
    u[0] = u[0] - l[0] * k2 / k0;
    d[0] = d[0] - l[0] * d0 / k0;
    l[0] = 0.0;
    // j=M-1 upper boundary condition
    double X3 = gridX_[Mxy - 3];
    X2 = gridX_[Mxy - 2];
    X1 = gridX_[Mxy - 1];
    conv1 = 0.5 * kappa1_ * X1;
    rate = 0.25 * (alpha + X1 + Y - kappa1_ - kappa2_);
    double kl3 = (-conv1 * (X1 - X2) - difu1) / (X2 - X3) / (X1 - X3);
    double kl2 = (conv1 * (X1 - X3) + difu1) / (X2 - X3) / (X1 - X2);
    double kl1 =
        1.0 / delT + rate +
        (-conv1 * (2.0 * X1 - X2 - X3) - difu1) / (X1 - X2) / (X1 - X3);
    double dl1 =
        (-conv1 * (X1 - X2) - difu1) / (X2 - X3) / (X1 - X3) * inM[Mxy - 3][k] +
        (conv1 * (X1 - X3) + difu1) / (X2 - X3) / (X1 - X2) * inM[Mxy - 2][k] +
        (rate +
         (-conv1 * (2.0 * X1 - X2 - X3) - difu1) / (X1 - X2) / (X1 - X3)) *
            inM[Mxy - 1][k] +
        midM[Mxy - 1][k] / delT;
    l[Mxy - 3] = l[Mxy - 3] - u[Mxy - 3] * kl3 / kl1;
    c[Mxy - 3] = c[Mxy - 3] - u[Mxy - 3] * kl2 / kl1;
    d[Mxy - 3] = d[Mxy - 3] - u[Mxy - 3] * dl1 / kl1;
    u[Mxy - 3] = 0.0;
    TriDiagonalSolve(Mxy - 2, l, c, u, d, V);
    for (int j = 1; j < Mxy - 1; j++)
      outM[j][k] = V[j - 1];
    // update j=0 lower boundary
    outM[0][k] = (d0 - k1 * outM[1][k] - k2 * outM[2][k]) / k0;
    // update j=M-1 upper boundary
    outM[Mxy - 1][k] =
        (dl1 - kl3 * outM[Mxy - 3][k] - kl2 * outM[Mxy - 2][k]) / kl1;
  }
};
void ShortRate2FPDE::oneStepForwardDouglasY(int T,
                                            const vector<vector<double>> &inM,
                                            const vector<vector<double>> &midM,
                                            vector<vector<double>> &outM) {
  vector<double> l(Mxy - 2), c(Mxy - 2), u(Mxy - 2), d(Mxy - 2), V(Mxy - 2);
  double Te = gridT_[T + 1];
  double Ti = gridT_[T];
  double delT = Te - Ti;
  double Tm = 0.5 * (Te + Ti);
  double sigma2 = whichValue(Tm, timeSigma2s_, sigma2s_);
  double difu2 = 0.5 * sigma2 * sigma2;
  double alpha = whichValue(Tm, timeAlphas_, alphas_);
  double conv2, rate;

  // for all x (j=0,...,M-1)
  for (int j = 0; j < Mxy; j++) {
    double X = gridX_[j];
    // y in the middle range
    for (int k = 1; k < Mxy - 1; k++) {
      double Y = gridY_[k];
      double Yu = gridY_[k + 1];
      double Yl = gridY_[k - 1];
      conv2 = 0.5 * (lambda_ * X + kappa2_ * Y);
      rate = 0.25 * (alpha + X + Y - kappa1_ - kappa2_);
      l[k - 1] = (conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl);
      c[k - 1] = 1.0 / delT + rate +
                 (-conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y);
      u[k - 1] = (-conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl);
      d[k - 1] =
          (conv2 * (Yu - Y) - difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
          (rate +
           (-conv2 * (Yu - 2.0 * Y + Yl) + difu2) / (Y - Yl) / (Yu - Y)) *
              inM[j][k] +
          (-conv2 * (Y - Yl) - difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1] +
          midM[j][k] / delT;
    }
    // k=0 lower boundary condition
    double Y0 = gridY_[0];
    double Y1 = gridY_[1];
    double Y2 = gridY_[2];
    conv2 = 0.5 * (lambda_ * X + kappa2_ * Y0);
    rate = 0.25 * (alpha + X + Y0 - kappa1_ - kappa2_);
    double k0 = 1.0 / delT + rate +
                (conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0);
    double k1 = (-conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1);
    double k2 = (conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0);
    double d0 =
        (rate +
         (conv2 * (Y2 + Y1 - 2.0 * Y0) - difu2) / (Y1 - Y0) / (Y2 - Y0)) *
            inM[j][0] +
        (-conv2 * (Y2 - Y0) + difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
        (conv2 * (Y1 - Y0) - difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2] +
        midM[j][0] / delT;
    c[0] = c[0] - l[0] * k1 / k0;
    u[0] = u[0] - l[0] * k2 / k0;
    d[0] = d[0] - l[0] * d0 / k0;
    l[0] = 0.0;
    // k=M-1 upper boundary condition
    double Y3 = gridY_[Mxy - 3];
    Y2 = gridY_[Mxy - 2];
    Y1 = gridY_[Mxy - 1];
    conv2 = 0.5 * (lambda_ * X + kappa2_ * Y1);
    rate = 0.25 * (alpha + X + Y1 - kappa1_ - kappa2_);
    double kl3 = (-conv2 * (Y1 - Y2) - difu2) / (Y2 - Y3) / (Y1 - Y3);
    double kl2 = (conv2 * (Y1 - Y3) + difu2) / (Y2 - Y3) / (Y1 - Y2);
    double kl1 =
        1.0 / delT + rate +
        (-conv2 * (2.0 * Y1 - Y2 - Y3) - difu2) / (Y1 - Y2) / (Y1 - Y3);
    double dl1 =
        (-conv2 * (Y1 - Y2) - difu2) / (Y2 - Y3) / (Y1 - Y3) * inM[j][Mxy - 3] +
        (conv2 * (Y1 - Y3) + difu2) / (Y2 - Y3) / (Y1 - Y2) * inM[j][Mxy - 2] +
        (rate +
         (-conv2 * (2.0 * Y1 - Y2 - Y3) - difu2) / (Y1 - Y2) / (Y1 - Y3)) *
            inM[j][Mxy - 1] +
        midM[j][Mxy - 1] / delT;
    l[Mxy - 3] = l[Mxy - 3] - u[Mxy - 3] * kl3 / kl1;
    c[Mxy - 3] = c[Mxy - 3] - u[Mxy - 3] * kl2 / kl1;
    d[Mxy - 3] = d[Mxy - 3] - u[Mxy - 3] * dl1 / kl1;
    u[Mxy - 3] = 0.0;
    TriDiagonalSolve(Mxy - 2, l, c, u, d, V);
    for (int k = 1; k < Mxy - 1; k++)
      outM[j][k] = V[k - 1];
    // update k=0 lower boundary
    outM[j][0] = (d0 - k1 * outM[j][1] - k2 * outM[j][2]) / k0;
    // update k=M-1 upper boundary
    outM[j][Mxy - 1] =
        (dl1 - kl3 * outM[j][Mxy - 3] - kl2 * outM[j][Mxy - 2]) / kl1;
  }
};

double ShortRate2FPDE::whichValue(double t, const vector<double> times,
                                  const vector<double> values) {
  int n = times.size();
  if (t < times[0])
    return values[0];
  if (t >= times[n - 1])
    return values[n - 1];
  for (int i = 1; i < n; i++)
    if ((t >= times[i - 1]) && (t < times[i]))
      return values[i];
};
} // namespace velesquant