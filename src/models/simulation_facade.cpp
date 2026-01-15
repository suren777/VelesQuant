#include <memory>
#include <vector>
#include <velesquant/volatility/l_vol.h>
#include <velesquant/volatility/s_vol.h>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/simulation_facade.h>
#include <velesquant/models/utility.h>

namespace velesquant {

using namespace std;

// Helper to convert Matrix to std::vector (duplicated from local_vol_facade,
// ideally common util)
std::vector<double> matrixToVectorSim(const Matrix &m) {
  std::vector<double> v;
  int rows = m.rows();
  int cols = m.cols();
  v.reserve(rows * cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      v.push_back(m(i, j));
    }
  }
  return v;
}

// Struct just to hold vectors now, or just locals? Locals are easier.
// Removing SABRParams struct entirely to avoid confusion and use locals.

Matrix simulationMLV(Matrix u1Sabr, double u1Spot, Matrix u2Sabr, double u2Spot,
                     Matrix u3Sabr, double u3Spot, Matrix corrMatrix,
                     Matrix theTimes, int NP) {
  // u1Sabr is likely 6 rows x N columns ?? Or 1 x 6?
  // model.cpp: SABRParams p1(u1Sabr); -> usage implies u1Sabr contains all
  // params. The SABRParams ctor in model.cpp iterated `n = theParams.cols()`
  // and filled params. This implies u1Sabr is a Matrix where rows are: 0:Mat,
  // 1:Fwd, 2:Beta, 3:Alpha, 4:Nu, 5:Rho. I must unpack this manually.

  auto unpackSABR = [](const Matrix &m, std::vector<double> &mat,
                       std::vector<double> &fwd, std::vector<double> &beta,
                       std::vector<double> &alpha, std::vector<double> &nu,
                       std::vector<double> &rho) {
    int n = m.cols();
    if (m.rows() < 6)
      throw std::runtime_error("SABR Params Matrix must have 6 rows");
    mat.resize(n);
    fwd.resize(n);
    beta.resize(n);
    alpha.resize(n);
    nu.resize(n);
    rho.resize(n);
    for (int i = 0; i < n; ++i) {
      mat[i] = m(0, i);
      fwd[i] = m(1, i);
      beta[i] = m(2, i);
      alpha[i] = m(3, i);
      nu[i] = m(4, i);
      rho[i] = m(5, i);
    }
  };

  std::vector<double> m1, f1, b1, a1, n1, r1;
  unpackSABR(u1Sabr, m1, f1, b1, a1, n1, r1);
  auto model1LV = std::make_unique<lVol>(m1, f1, b1, a1, n1, r1, u1Spot);

  std::vector<double> m2, f2, b2, a2, n2, r2;
  unpackSABR(u2Sabr, m2, f2, b2, a2, n2, r2);
  auto model2LV = std::make_unique<lVol>(m2, f2, b2, a2, n2, r2, u2Spot);

  std::vector<double> m3, f3, b3, a3, n3, r3;
  unpackSABR(u3Sabr, m3, f3, b3, a3, n3, r3);
  auto model3LV = std::make_unique<lVol>(m3, f3, b3, a3, n3, r3, u3Spot);

  int m = corrMatrix.rows();
  vector<vector<double>> matrixCorr(m, vector<double>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++)
      matrixCorr[i][j] = corrMatrix(i, j);
  }
  vector<vector<double>> matrixCholesky = cholesky(matrixCorr);

  std::vector<double> times;
  if (theTimes.rows() > 0) {
    times.resize(theTimes.rows());
    for (int i = 0; i < theTimes.rows(); ++i)
      times[i] = theTimes(i, 0);
  }
  int n = times.size();

  Matrix theSpots(3 * n, NP);
  for (int j = 0; j < NP; j++) {
    vector<double> rands1(n), rands2(n), rands3(n);
    for (int i = 0; i < n; i++) {
      double ra1 = random_normal();
      double ra2 = random_normal();
      double ra3 = random_normal();
      rands1[i] = ra1;
      rands2[i] = matrixCholesky[1][0] * ra1 + matrixCholesky[1][1] * ra2;
      rands3[i] = matrixCholesky[2][0] * ra1 + matrixCholesky[2][1] * ra2 +
                  matrixCholesky[2][2] * ra3;
    }
    std::vector<double> spot1s(model1LV->simulation(times, rands1));
    for (int i = 0; i < n; i++)
      theSpots(i, j) = spot1s[i];
    std::vector<double> spot2s(model2LV->simulation(times, rands2));
    for (int i = 0; i < n; i++)
      theSpots(i + n, j) = spot2s[i];
    std::vector<double> spot3s(model3LV->simulation(times, rands3));
    for (int i = 0; i < n; i++)
      theSpots(i + n + n, j) = spot3s[i];
  }
  return theSpots;
}

Matrix simulation2LV(Matrix the1Maturities, Matrix the1Forwards,
                     Matrix the1Betas, Matrix the1Alphas, Matrix the1Nus,
                     Matrix the1Rhos, double the1Spot, Matrix the2Maturities,
                     Matrix the2Forwards, Matrix the2Betas, Matrix the2Alphas,
                     Matrix the2Nus, Matrix the2Rhos, double the2Spot,
                     double theCorr, Matrix theTimes, int NP) {
  std::vector<double> m1 = matrixToVectorSim(the1Maturities);
  std::vector<double> f1 = matrixToVectorSim(the1Forwards);
  std::vector<double> b1 = matrixToVectorSim(the1Betas);
  std::vector<double> a1 = matrixToVectorSim(the1Alphas);
  std::vector<double> n1 = matrixToVectorSim(the1Nus);
  std::vector<double> r1 = matrixToVectorSim(the1Rhos);

  auto model1LV = std::make_unique<lVol>(m1, f1, b1, a1, n1, r1, the1Spot);

  std::vector<double> m2 = matrixToVectorSim(the2Maturities);
  std::vector<double> f2 = matrixToVectorSim(the2Forwards);
  std::vector<double> b2 = matrixToVectorSim(the2Betas);
  std::vector<double> a2 = matrixToVectorSim(the2Alphas);
  std::vector<double> n2 = matrixToVectorSim(the2Nus);
  std::vector<double> r2 = matrixToVectorSim(the2Rhos);

  auto model2LV = std::make_unique<lVol>(m2, f2, b2, a2, n2, r2, the2Spot);

  std::vector<double> times;
  if (theTimes.rows() > 0) {
    times.resize(theTimes.rows());
    for (int i = 0; i < theTimes.rows(); ++i)
      times[i] = theTimes(i, 0);
  }
  int n = times.size();

  Matrix theSpots(2 * n, NP);
  for (int j = 0; j < NP; j++) {
    std::vector<double> rands1(n), rands2(n);
    double cor2 = sqrt(1.0 - theCorr * theCorr);
    for (int i = 0; i < n; i++) {
      rands1[i] = random_normal();
      rands2[i] = theCorr * rands1[i] + cor2 * random_normal();
    }
    std::vector<double> spot1s(model1LV->simulation(times, rands1));
    for (int i = 0; i < n; i++)
      theSpots(i, j) = spot1s[i];
    std::vector<double> spot2s(model2LV->simulation(times, rands2));
    for (int i = 0; i < n; i++)
      theSpots(i + n, j) = spot2s[i];
  }
  return theSpots;
}

Matrix simoreLV(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos,
                Matrix theTimes, double spot, int NP) {
  std::vector<double> m = matrixToVectorSim(theMaturities);
  std::vector<double> f = matrixToVectorSim(theForwards);
  std::vector<double> b = matrixToVectorSim(theBetas);
  std::vector<double> a = matrixToVectorSim(theAlphas);
  std::vector<double> n = matrixToVectorSim(theNus);
  std::vector<double> r = matrixToVectorSim(theRhos);

  auto modelLV = std::make_unique<lVol>(m, f, b, a, n, r, spot);

  std::vector<double> times;
  if (theTimes.rows() > 0) {
    times.resize(theTimes.rows());
    for (int i = 0; i < theTimes.rows(); ++i)
      times[i] = theTimes(i, 0);
  }
  int sz = times.size();

  Matrix theSpots(sz, NP);
  for (int j = 0; j < NP; j++) {
    std::vector<double> spots(modelLV->simulation(times));
    for (int i = 0; i < sz; i++)
      theSpots(i, j) = spots[i];
  }
  return theSpots;
}

} // namespace velesquant
