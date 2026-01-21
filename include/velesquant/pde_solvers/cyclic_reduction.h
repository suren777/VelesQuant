#ifndef CYCLIC_REDUCTION_H
#define CYCLIC_REDUCTION_H

#include <array>
#include <vector>

namespace velesquant {
// using namespace std; // REMOVED by refactor
template <typename T, size_t N>
void TriDiagonalSolve(const int n, const std::array<T, N> &a,
                      const std::array<T, N> &b, const std::array<T, N> &c,
                      const std::array<T, N> &d, std::vector<double> &X) {
  std::vector<double> TempC(n);
  std::vector<double> TempD(n);

  for (int i = 0; i < n; i++) {
    TempC[i] = c[i];
    TempD[i] = d[i];
  }

  /* Modify the coefficients. */
  TempC[0] /= b[0];
  TempD[0] /= b[0];
  for (int i = 1; i < n; i++) {
    double id = 1 / (b[i] - TempC[i - 1] * a[i]);
    TempC[i] *= id;
    TempD[i] = (TempD[i] - TempD[i - 1] * a[i]) * id;
  }
  /* Now back substitute. */
  X[n - 1] = TempD[n - 1];
  for (int i = n - 2; i >= 0; i--)
    X[i] = TempD[i] - TempC[i] * X[i + 1];
}

template <typename T>
void TriDiagonalSolve(const int n, const std::vector<T> &a,
                      const std::vector<T> &b, const std::vector<T> &c,
                      const std::vector<T> &d, std::vector<T> &X) {
  std::vector<double> TempC(n);
  std::vector<double> TempD(n);

  for (int i = 0; i < n; i++) {
    TempC[i] = c[i];
    TempD[i] = d[i];
  }

  /* Modify the coefficients. */
  TempC[0] /= b[0];
  TempD[0] /= b[0];
  for (int i = 1; i < n; i++) {
    double id = 1 / (b[i] - TempC[i - 1] * a[i]);
    TempC[i] *= id;
    TempD[i] = (TempD[i] - TempD[i - 1] * a[i]) * id;
  }
  /* Now back substitute. */
  X[n - 1] = TempD[n - 1];
  for (int i = n - 2; i >= 0; i--)
    X[i] = TempD[i] - TempC[i] * X[i + 1];
}

template <typename T>
void TriDiagonalSolve(const int n, const std::vector<T> &a,
                      const std::vector<T> &b, const std::vector<T> &c,
                      const std::vector<T> &d, std::vector<T> &X,
                      std::vector<T> &TempC, std::vector<T> &TempD) {
  // Optimization: Remove OMP and memory allocation
  for (int i = 0; i < n; i++) {
    TempC[i] = c[i];
    TempD[i] = d[i];
  }

  /* Modify the coefficients. */
  TempC[0] /= b[0];
  TempD[0] /= b[0];
  for (int i = 1; i < n; i++) {
    double id = 1.0 / (b[i] - TempC[i - 1] * a[i]);
    TempC[i] *= id;
    TempD[i] = (TempD[i] - TempD[i - 1] * a[i]) * id;
  }
  /* Now back substitute. */
  X[n - 1] = TempD[n - 1];
  for (int i = n - 2; i >= 0; i--)
    X[i] = TempD[i] - TempC[i] * X[i + 1];
}

template <typename T, size_t N>
void CyclicReduction(std::vector<double> &a, std::vector<double> &b,
                     std::vector<double> &c, std::vector<double> &d,
                     std::vector<double> &x) {
  int n_size = x.size();
  std::vector<double> na(n_size, 0), nb(n_size - 2, 0), nc(n_size, 0),
      nd(n_size, 0);

  int i, id, auxN = n_size;
  int l = 0;

  // forward substitution
  do {
    id = pow(2, l - 1);
    i = 0;
    double gamma = -c[i] / b[i + id];
    nb[id] = b[i] + gamma * a[i + id];
    nc[id] = gamma * c[i + id];
    nd[id] = d[i] + gamma * d[i + id];

    // #pragma omp parallel for
    for (i = 1; i < n_size - 1; i++) {

      double alpha_cr = -a[i];
      double gamma_cr = -c[i] / b[i + id];
      na[id] = -a[i - id] * alpha_cr;
      nb[id] = b[i] + alpha_cr * c[i - id] + gamma_cr * a[i + id];
      nc[id] = gamma_cr * c[i + id];
      nd[id] = d[i] + alpha_cr * d[i - id] + gamma_cr * d[i + id];
    }
    i = n_size - 1;
    double alpha = -a[i] / b[i - id];
    na[id] = -a[i - id] * alpha;
    nb[id] = b[i] + alpha * c[i - id];
    nd[id] = d[i] + alpha * d[i - id];
    auxN /= 2;
    na.swap(a);
    nb.swap(b);
    nc.swap(c);
    nd.swap(d);
  } while (auxN > 1);
  // #pragma omp parallel for
  for (int k = 0; k < n_size; k++)
    x[k] = d[k] / b[k];
}
} // namespace velesquant

#endif // CYCLIC_REDUCTION_H