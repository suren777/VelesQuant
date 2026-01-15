#pragma once

#include <vector>
#include <velesquant/errors.h>
#include <velesquant/types.h>

namespace velesquant {

// Templated helper to extract a column from CellMatrix into std::vector<double>
// Templated helper to extract a column from CellMatrix into std::vector<double>
template <typename T = double>
std::vector<T> extractColumn(const CellMatrix &matrix, int columnIndex) {
  int rows = matrix.RowsInStructure();
  std::vector<T> result(rows);
  for (int i = 0; i < rows; ++i) {
    result[i] = static_cast<T>(matrix(i, columnIndex).NumericValue());
  }
  return result;
}

// Overload for Eigen::MatrixXd
template <typename T = double>
std::vector<T> extractColumn(const Matrix &matrix, int columnIndex) {
  int rows = matrix.rows();
  std::vector<T> result(rows);
  for (int i = 0; i < rows; ++i) {
    result[i] = static_cast<T>(matrix(i, columnIndex));
  }
  return result;
}

// Struct to hold SABR parameters to reduce repetition
struct SABRParams {
  std::vector<double> maturities;
  std::vector<double> forwards;
  std::vector<double> betas;
  std::vector<double> alphas;
  std::vector<double> nus;
  std::vector<double> rhos;

  // Constructor extracting from standard CellMatrix inputs
  SABRParams(const CellMatrix &m, const CellMatrix &f, const CellMatrix &b,
             const CellMatrix &a, const CellMatrix &n, const CellMatrix &r) {
    maturities = extractColumn(m, 0);
    forwards = extractColumn(f, 0);
    betas = extractColumn(b, 0);
    alphas = extractColumn(a, 0);
    nus = extractColumn(n, 0);
    rhos = extractColumn(r, 0);
  }

  // Constructor extracting from Eigen::MatrixXd inputs
  SABRParams(const Matrix &m, const Matrix &f, const Matrix &b, const Matrix &a,
             const Matrix &n, const Matrix &r) {
    maturities = extractColumn(m, 0);
    forwards = extractColumn(f, 0);
    betas = extractColumn(b, 0);
    alphas = extractColumn(a, 0);
    nus = extractColumn(n, 0);
    rhos = extractColumn(r, 0);
  }

  // Constructor extracting from single combined matrix (if applicable)
  SABRParams(const CellMatrix &combined) {
    int m = combined.RowsInStructure();
    int cols = combined.ColumnsInStructure();
    VEL_CHECK(cols >= 6, "Combined SABR matrix must have at least 6 columns");

    maturities.resize(m);
    forwards.resize(m);
    betas.resize(m);
    alphas.resize(m);
    nus.resize(m);
    rhos.resize(m);

    for (int i = 0; i < m; ++i) {
      maturities[i] = combined(i, 0).NumericValue();
      forwards[i] = combined(i, 1).NumericValue();
      betas[i] = combined(i, 2).NumericValue();
      alphas[i] = combined(i, 3).NumericValue();
      nus[i] = combined(i, 4).NumericValue();
      rhos[i] = combined(i, 5).NumericValue();
    }
  }

  // Constructor extracting from single combined Eigen::MatrixXd
  SABRParams(const Matrix &combined) {
    int m = combined.rows();
    int cols = combined.cols();
    VEL_CHECK(cols >= 6, "Combined SABR matrix must have at least 6 columns");

    maturities.resize(m);
    forwards.resize(m);
    betas.resize(m);
    alphas.resize(m);
    nus.resize(m);
    rhos.resize(m);

    for (int i = 0; i < m; ++i) {
      maturities[i] = combined(i, 0);
      forwards[i] = combined(i, 1);
      betas[i] = combined(i, 2);
      alphas[i] = combined(i, 3);
      nus[i] = combined(i, 4);
      rhos[i] = combined(i, 5);
    }
  }
};

} // namespace velesquant
