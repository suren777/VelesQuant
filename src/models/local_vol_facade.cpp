#include <memory>
#include <vector>
#include <velesquant/volatility/l_vol.h>
#include <velesquant/volatility/sabr.h>
#include <velesquant/models/local_vol_facade.h>

namespace velesquant {

using namespace std;

// Helper to convert Matrix to std::vector
std::vector<double> matrixToVector(const Matrix &m) {
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

Matrix theDensity(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                  Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                  double maturity) {
  std::vector<double> vMaturities = matrixToVector(theMaturities);
  std::vector<double> vForwards = matrixToVector(theForwards);
  std::vector<double> vBetas = matrixToVector(theBetas);
  std::vector<double> vAlphas = matrixToVector(theAlphas);
  std::vector<double> vNus = matrixToVector(theNus);
  std::vector<double> vRhos = matrixToVector(theRhos);

  auto modelLV = std::make_unique<lVol>(vMaturities, vForwards, vBetas, vAlphas,
                                        vNus, vRhos, spot);
  std::vector<double> aDensity = modelLV->density(maturity, 100);

  int n = aDensity.size();
  Matrix density(n, 1);
  for (int i = 0; i < n; i++)
    density(i, 0) = aDensity[i];
  return density;
}

Matrix lvExport(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                Matrix theTimes) {
  std::vector<double> vMaturities = matrixToVector(theMaturities);
  std::vector<double> vForwards = matrixToVector(theForwards);
  std::vector<double> vBetas = matrixToVector(theBetas);
  std::vector<double> vAlphas = matrixToVector(theAlphas);
  std::vector<double> vNus = matrixToVector(theNus);
  std::vector<double> vRhos = matrixToVector(theRhos);

  auto modelLV = std::make_unique<lVol>(vMaturities, vForwards, vBetas, vAlphas,
                                        vNus, vRhos, spot);

  std::vector<double> times;
  if (theTimes.rows() > 0) {
    times.resize(theTimes.rows());
    for (int i = 0; i < theTimes.rows(); ++i)
      times[i] = theTimes(i, 0);
  } else if (theTimes.cols() > 0) {
    times.resize(theTimes.cols());
    for (int i = 0; i < theTimes.cols(); ++i)
      times[i] = theTimes(0, i);
  }

  vector<vector<double>> LV = modelLV->exportLV(times);

  if (LV.empty())
    return Matrix(0, 0);

  int n = LV.size();
  int m = LV[0].size();
  Matrix theLV(n, m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      theLV(i, j) = LV[i][j];
  return theLV;
}

double dntPDElv(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                double maturity, double upperBarrier, double lowerBarrier) {
  std::vector<double> vMaturities = matrixToVector(theMaturities);
  std::vector<double> vForwards = matrixToVector(theForwards);
  std::vector<double> vBetas = matrixToVector(theBetas);
  std::vector<double> vAlphas = matrixToVector(theAlphas);
  std::vector<double> vNus = matrixToVector(theNus);
  std::vector<double> vRhos = matrixToVector(theRhos);

  auto modelLV = std::make_unique<lVol>(vMaturities, vForwards, vBetas, vAlphas,
                                        vNus, vRhos, spot);
  return modelLV->dntPDE(maturity, upperBarrier, lowerBarrier);
}

double putPDElv(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                double maturity, double strike) {
  std::vector<double> vMaturities = matrixToVector(theMaturities);
  std::vector<double> vForwards = matrixToVector(theForwards);
  std::vector<double> vBetas = matrixToVector(theBetas);
  std::vector<double> vAlphas = matrixToVector(theAlphas);
  std::vector<double> vNus = matrixToVector(theNus);
  std::vector<double> vRhos = matrixToVector(theRhos);

  auto modelLV = std::make_unique<lVol>(vMaturities, vForwards, vBetas, vAlphas,
                                        vNus, vRhos, spot);
  return modelLV->putPDE(maturity, strike);
}

Matrix callPDElv(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                 Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                 double maturity, double strike) {
  std::vector<double> vMaturities = matrixToVector(theMaturities);
  std::vector<double> vForwards = matrixToVector(theForwards);
  std::vector<double> vBetas = matrixToVector(theBetas);
  std::vector<double> vAlphas = matrixToVector(theAlphas);
  std::vector<double> vNus = matrixToVector(theNus);
  std::vector<double> vRhos = matrixToVector(theRhos);

  auto modelLV = std::make_unique<lVol>(vMaturities, vForwards, vBetas, vAlphas,
                                        vNus, vRhos, spot);
  double val = modelLV->callPDE(maturity, strike);
  Matrix res(1, 1);
  res(0, 0) = val;
  return res;
}

} // namespace velesquant
