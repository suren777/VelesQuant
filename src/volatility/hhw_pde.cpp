#include <velesquant/volatility/hhw_pde.h>
// using namespace std;
#include <algorithm>

namespace velesquant {

double asinh(double a) { return std::log(a + std::sqrt(a * a + 1.0)); }

double HHWPDE::deltax(Vdoub aVec, int i) {
  if ((i - 1) < 0)
    throw("Memory Violation\n");
  return aVec[i] - aVec[i - 1];
};

double HHWPDE::alpha(Vdoub aVec, int i, int type) {
  double ret = 0;
  double a = deltax(aVec, i);
  double b = deltax(aVec, i - 1);
  if (type == -2)
    ret = a / (a * (a - b));
  else if (type == -1)
    ret = (-b - a) / (a * b);
  else if (type == 0)
    ret = (b + a + a) / (b * (a + b));
  else
    throw("Wrong type in Alpha\n");
  return ret;
};

double HHWPDE::beta(Vdoub aVec, int i, int type) {
  double ret = 0;
  double a = deltax(aVec, i + 1);
  double b = deltax(aVec, i);
  if (type == -1)
    ret = -a / (b * (a + b));
  else if (type == 0)
    ret = (a - b) / (a * b);
  else if (type == 1)
    ret = b / (a * (a + b));
  else
    throw("Wrong type in Beta\n");
  return ret;
};

double HHWPDE::gamma(Vdoub aVec, int i, int type) {
  double ret = 0;
  double a = deltax(aVec, i + 2);
  double b = deltax(aVec, i + 1);
  if (type == 0)
    ret = (-b - b - a) / (b * (a + b));
  else if (type == 1)
    ret = (a + b) / (a * b);
  else if (type == 2)
    ret = -b / (a * (a + b));
  else
    throw("Wrong type in Gamma\n");
  return ret;
};

double HHWPDE::delta(Vdoub aVec, int i, int type) {
  double ret = 0;
  double a = deltax(aVec, i + 1);
  double b = deltax(aVec, i);
  if (type == -1)
    ret = 2.0 / (b * (a + b));
  else if (type == 0)
    ret = -2.0 / (a * b);
  else if (type == 1)
    ret = 2.0 / (a * (a + b));
  else
    throw("Wrong type in Delta\n");
  return ret;
};

void HHWPDE::generateMeshes() {
  /*
  i = r index
  j = v index
  k = s index
  */
  N_ = (Mr_ + 1) * (Mv_ + 1) * (Ms_ + 1);

  sM_ = gsMesh();
  rM_ = grMesh();
  vM_ = gvMesh();

  for (int i = 0; i < Mr_; i++) {
    bRm1[i] = beta(rM_, i, -1);
    bR01[i] = beta(rM_, i, 0);
    bRp1[i] = beta(rM_, i, 1);
    dRm1[i] = beta(rM_, i, -1);
    dR01[i] = beta(rM_, i, 0);
    dRp1[i] = beta(rM_, i, 1);
  }
  for (int j = 0; j < Mv_; j++) {
    bVm1[j] = beta(vM_, j, -1);
    bV01[j] = beta(vM_, j, 0);
    bVp1[j] = beta(vM_, j, 1);
    dVm1[j] = beta(vM_, j, -1);
    dV01[j] = beta(vM_, j, 0);
    dVp1[j] = beta(vM_, j, 1);
  }
  for (int k = 0; k < Ms_; k++) {
    bSm1[k] = beta(rM_, k, -1);
    bS01[k] = beta(rM_, k, 0);
    bSp1[k] = beta(rM_, k, 1);
    dSm1[k] = beta(rM_, k, -1);
    dS01[k] = beta(rM_, k, 0);
    dSp1[k] = beta(rM_, k, 1);
  }
};

Vdoub HHWPDE::gsMesh() {
  Vdoub Mesh;
  Mesh.resize(Ms_ + 1);
  double Sleft = .5 * K_;
  double Sright = K_;
  double d1 = K_ / 20.0;
  double zmin = asinh(-Sleft / d1);
  double zint = (Sright - Sleft) / d1;
  double zmax = zint + asinh((Smax_ - Sright) / d1);
  double dz = (zmax - zint) / Ms_;

  for (int i = 0; i <= Ms_; i++) {
    double z = zmin + i * dz;
    if ((z >= zmin) && (z < 0))
      Mesh[i] = Sleft + d1 * std::sinh(z);
    else if ((z >= 0) && (z <= zint))
      Mesh[i] = Sleft + d1 * z;
    else
      Mesh[i] = Sright + d1 * std::sinh(z - zint);
  }
  return Mesh;
};

Vdoub HHWPDE::gvMesh() {
  Vdoub Mesh;
  Mesh.resize(Mv_ + 1);
  double d2 = Vmax_ / 500.0;
  double dn = 1.0 / Mv_ * asinh(Vmax_ / d2);
  for (int i = 0; i <= Mv_; i++)
    Mesh[i] = d2 * std::sinh(i * dn);
  return Mesh;
};

Vdoub HHWPDE::grMesh() {
  Vdoub Mesh;
  Mesh.resize(Mv_ + 1);
  double d3 = Rmax_ / 500.0;
  double start = asinh((-Rmax_ - r0_) / d3);
  double dx = 1 / Mr_ * (asinh((Rmax_ - r0_) / d3) - start);
  for (int i = 0; i <= Mv_; i++)
    Mesh[i] = r0_ + d3 * std::sinh(start + i * dx);
  return Mesh;
};

void HHWPDE::aAbB(double a, Vdoub &Ain, double b, Vdoub &Bin, Vdoub &Cout) {
  size_t n;
  if ((n = Ain.size()) == Bin.size())
    for (size_t i = 0; i < n; i++)
      Cout[i] = a * Ain[i] + b * Bin[i];
  else
    throw("Error: dim(A) not equal to dim(B)\n");
};

void HHWPDE::A0Ug(Vdoub &Uin, Vdoub &Uout) {
  for (int i = 0; i <= Mr_; i++)
    for (int j = 0; j <= Mv_; j++)
      for (int k = 0; k <= Ms_; k++)
        Uout[i + Mv_ * j + Mv_ * Ms_ * k] =
            rho12_ * sigma1_ * vM_[j] * sM_[k] *
                dudxdy(vM_, sM_, Uin, i, j, k, 1) +
            rho13_ * sigma2_ * sM_[k] * std::sqrt(vM_[j]) *
                dudxdy(rM_, sM_, Uin, i, j, k, 3) +
            rho23_ * sigma1_ * sigma2_ * std::sqrt(vM_[j]) *
                dudxdy(rM_, vM_, Uin, i, j, k, 2);
};
double HHWPDE::dudxdy(Vdoub &a, Vdoub &b, Vdoub &U, int i, int j, int k,
                      int type) {
  /*
  Types
  1: d^2/(dvds)
  2: d^2/(drdv)
  3: d^2/(drds)

  i = r index
  j = v index
  k = s index

  */

  double Sum = 0;
  //-----------------Boundary conditions:---------------

  if ((type == 3) && ((i == 0) || (i == Mr_) || (j == Mv_) || (j == 0) ||
                      (k == 0) || (k == Ms_)))
    return 0;

  if ((type == 1) && ((j == Mv_) || (j == 0) || (k == 0) || (k == Ms_)))
    return 0;

  if ((type == 2) &&
      ((i == 0) || (i == Mr_) || (j == Mv_) || (j == 0) || (k == 0)))
    return 0;
  //---------------------------------------------------------//
  if (type == 1)
    for (int l = -1; l <= 1; l++)
      for (int m = -1; m <= 1; m++) {
        if (((i + l) < 0) || ((i + l) > Mr_))
          continue;
        else
          Sum += beta(a, j, l) * beta(b, k, m) * fU(U, i, j + l, k + m);
      }
  else if (type == 2)
    for (int l = -1; l <= 1; l++)
      for (int m = -1; m <= 1; m++) {
        if (k > Ms_)
          continue;
        else
          Sum += beta(a, i, l) * beta(b, j, m) * fU(U, i + l, j + m, k);
      }
  else if (type == 3)
    for (int l = -1; l <= 1; l++)
      for (int m = -1; m <= 1; m++)
        Sum += beta(a, i, l) * beta(b, k, m) * fU(U, i + l, j, k + m);
  else
    throw("Wrong type\n");
  return Sum;
};

void HHWPDE::LU3dec(Vdoub &a, Vdoub &b, Vdoub &c, Vdoub &l, Vdoub &d) {

  d[0] = a[0];
  for (int i = 0; i < N_; i++) {
    l[i - 1] = b[i - 1] / d[i - 1];
    d[i] = a[i] - l[i - 1] * c[i - 1];
  }
};

void HHWPDE::LU3decSparce(Vdoub &c, Vdoub &a, Vdoub &b, Vdoub &l, Vdoub &d,
                          int /*Ns*/) {

  for (int i = 0; i < N_; i++) {
    if (i <= N_ / 2)
      d[i] = a[i];
    else {
      l[i - N_ / 2] = c[i - N_ / 2] / d[i - N_ / 2];
      d[i] = a[i] - l[i - N_ / 2] * b[i - N_ / 2];
    }
  }
};

void HHWPDE::LU5dec(Vdoub & /*l*/, Vdoub &m, Vdoub &n, Vdoub &o, Vdoub &p,
                    /*Outputs:*/ Vdoub &b, Vdoub &c, Vdoub &d, Vdoub &e) {

  d[0] = n[0];
  e[0] = o[0];
  c[0] = m[0] / d[0];
  d[1] = n[1] - c[0] * e[0];
  e[1] = o[1] - c[0] * p[0];
  for (int i = 2; i < N_; i++) {
    b[i - 2] = e[i - 2] / d[i - 2];
    c[i - 1] = (m[i - 1] - b[i - 2] * e[i - 2]) / d[i - 1];
    d[i] = n[i] - b[i - 2] * p[i - 2] + c[i - 1] * e[i - 1];
    if (i < N_ - 1)
      e[i] = o[i] - c[i - 1] * p[i - 1];
  }
};

void HHWPDE::A1Ug(Vdoub &Uin, Vdoub &Uout) {
  /*
  A_1 * U + g(t)
  Second derivatives of Spot

  i = r index
  j = v index
  k = s index

  */
  double sv, rs, L, M, R;

  for (int k = 0; k <= Ms_; k++)
    for (int j = 0; j <= Mv_; j++) {
      sv = 0.5 * sM_[k] * sM_[k] * vM_[j];
      for (int i = 0; i <= Mr_; i++) {
        rs = rM_[i] * sM_[k];
        if (i == 0)
          L = 0;
        else
          L = (sv * dSm1[k] + rs * bSm1[k]) * fU(Uin, i, j, k - 1);
        M = (sv * dS01[k] + rs * bS01[k] - 1.0 / 3.0 * rM_[i]) *
            fU(Uin, i, j, k);
        if (k == Ms_)
          R = 0;
        else
          R = (sv * dSp1[k] + rs * bSm1[k]) * fU(Uin, i, j, k + 1);

        Uout[i + Mv_ * j + Mv_ * Ms_ * k] = L + M + R;
      }
    }
};

void HHWPDE::A2Ug(Vdoub &Uin, Vdoub &Uout) {
  double L, M, R, sigv, knuv;
  for (int j = 0; j <= Mv_; j++) {
    sigv = 0.5 * sigma1_ * vM_[j];
    knuv = kappa_ * (eta_ - vM_[j]);
    for (int k = 0; k <= Ms_; k++)
      for (int i = 0; i <= Mr_; i++) {
        if (j == 0) {
          L = (knuv * gamma(vM_, j, 0) - 1.0 / 3.0 * rM_[i]) * fU(Uin, i, j, k);
          M = knuv * gamma(vM_, j, 1) * fU(Uin, i, j + 1, k);
          R = knuv * gamma(vM_, j, 2) * fU(Uin, i, j + 2, k);
        } else if ((vM_[j] <= eta_) && (j > 0)) {
          L = (sigv * dVm1[j] + knuv * bVm1[j]) * fU(Uin, i, j - 1, k);
          M = (sigv + knuv * bV01[j] - 1.0 / 3.0 * rM_[i]) * fU(Uin, i, j, k);
          R = (sigv + knuv * bVp1[j]) * fU(Uin, i, j + 1, k);
        } else if ((vM_[j] > eta_) && (j < Mv_)) {
          L = knuv * alpha(vM_, j, -2) * fU(Uin, i, j - 2, k);
          M = (sigv * dVm1[j] + knuv * alpha(vM_, j, -1)) *
              fU(Uin, i, j - 1, k);
          R = (sigv * dV01[j] + knuv * alpha(vM_, j, 0)) * fU(Uin, i, j, k) +
              sigv * dVp1[j] * fU(Uin, i, j + 1, k);
        } else {
          L = 0;
          M = sM_[k];
          R = 0;
        }

        Uout[i + Mv_ * j + Mv_ * Ms_ * k] = L + M + R;
      }
  }
};

void HHWPDE::A3Ug(double t, Vdoub &Uin, Vdoub &Uout) {
  /*Add boundary conditions*/
  double L, M, R;
  double s2 = .5 * sigma2_ * sigma2_;
  for (int i = 0; i <= Mr_; i++) {
    double abtr = a_ * (b(T_ - t) - rM_[i]);
    if (i == 0)
      L = 0;
    else
      L = s2 * dRm1[i] + abtr * bRm1[i];
    M = s2 * dR01[i] + abtr * bR01[i];
    if (i == Mr_)
      R = 0;
    else
      R = s2 * dRp1[i] + abtr * bRp1[i];
    for (int j = 0; j < Mv_; j++)
      for (int k = 0; k < Ms_; k++)
        Uout[i + Mv_ * j + Mv_ * Ms_ * k] =
            +((i == 0) ? 0 : L * fU(Uin, i - 1, j, k)) + M * fU(Uin, i, j, k) +
            ((i == Mr_) ? 0 : R * fU(Uin, i + 1, j, k));
  }
};

void HHWPDE::pentaSolve(Vdoub &E, Vdoub &A, Vdoub &D, Vdoub &C, Vdoub &F,
                        Vdoub &B, Vdoub &X, int /*t*/) {
  static double xmult;
  int N;
  N = D.size();
  for (int i = 1; i < N - 1; i++) {
    xmult = A[i - 1] / D[i - 1];
    D[i] = D[i] - xmult * C[i - 1];
    C[i] = C[i] - xmult * F[i - 1];
    B[i] = B[i] - xmult * B[i - 1];
    xmult = E[i - 1] / D[i - 1];
    D[i + 1] = D[i + 1] - xmult * F[i - 1];
    B[i + 1] = B[i + 1] - xmult * B[i - 1];

    xmult = A[N - 2] / D[N - 2];
    D[N - 1] = D[N - 1] - xmult * C[N - 2];
  }
  X[N - 1] = (B[N - 1] - xmult * B[N - 2]) / D[N - 1];
  X[N - 2] = (B[N - 2] - xmult * C[N - 2] * X[N - 1]) / D[N - 2];
  for (int i = N - 3; i >= 0; i--)
    X[i] = (B[i] - F[i] * X[i + 2] - C[i] * X[i + 1]) / D[i];
};

void HHWPDE::fillA1(Vdoub &A1l, Vdoub &A1d, Vdoub &A1u) {
  double sv, rs;

  int count = 0;
  for (int k = 0; k <= Ms_; k++)
    for (int j = 0; j <= Mv_; j++) {
      sv = 0.5 * sM_[k] * sM_[k] * vM_[j];
      for (int i = 0; i <= Mr_; i++) {
        count = i + Mv_ * j + Mv_ * Ms_ * k;
        rs = rM_[i] * sM_[k];
        if (i == 0)
          A1d[count] = (sv * dS01[k] + rs * bS01[k] - 1.0 / 3.0 * rM_[i]);
        else {
          A1l[count - Mv_ * Ms_] = (sv * dSm1[k] + rs * bSm1[k]);
          A1d[count] = (sv * dS01[k] + rs * bS01[k] - 1.0 / 3.0 * rM_[i]);
        }
        if (k < Ms_)
          A1u[count] = (sv * dSp1[k] + rs * bSm1[k]);
      }
    }
};

void HHWPDE::fillA2(Vdoub A2l1, Vdoub &A2l2, Vdoub &A2d, Vdoub &A2u1,
                    Vdoub &A2u2) {
  double sigv, knuv;
  int count;
  for (int j = 0; j <= Mv_; j++) {
    sigv = 0.5 * sigma1_ * vM_[j];
    knuv = kappa_ * (eta_ - vM_[j]);
    for (int k = 0; k <= Ms_; k++)
      for (int i = 0; i <= Mr_; i++) {
        count = i + Mv_ * j + Mv_ * Ms_ * k;
        if (j == 0) {
          A2d[count] = (knuv * gamma(vM_, j, 0) - 1.0 / 3.0 * rM_[i]);
          A2u1[count] = knuv * gamma(vM_, j, 1);
          A2u2[count] = knuv * gamma(vM_, j, 2);
        } else if ((vM_[j] <= eta_) && (j > 0)) {
          A2l2[count - Ms_ * Mr_] = (sigv * dVm1[j] + knuv * bVm1[j]);
          A2d[count] = (sigv + knuv * bV01[j] - 1.0 / 3.0 * rM_[i]);
          A2u1[count] = (sigv + knuv * bVp1[j]);
        } else if ((vM_[j] > eta_) && (j < Mv_)) {
          A2l1[count - 2 * Ms_ * Mr_] = knuv * alpha(vM_, j, -2);
          A2l2[count - Ms_ * Mr_] = (sigv * dVm1[j] + knuv * alpha(vM_, j, -1));
          A2d[count] = (sigv * dV01[j] + knuv * alpha(vM_, j, 0));
          A2u1[count] = sigv * dVp1[j];
        } else
          A2d[count] = 1;
      }
  }
};

void HHWPDE::fillA3(Vdoub &A3l, Vdoub &A3d, Vdoub &A3u, double t) {
  double s2 = .5 * sigma2_ * sigma2_;
  int count = 0;
  for (int i = 0; i <= Mr_; i++) {
    double abtr = a_ * (b(T_ - t) - rM_[i]);
    for (int j = 0; j < Mv_; j++)
      for (int k = 0; k < Ms_; k++) {
        count = i + Mv_ * j + Mv_ * Ms_ * k;
        if (i > 0)
          A3l[count - Ms_ * Mv_] = s2 * dRm1[i] + abtr * bRm1[i];
        A3d[count] = s2 * dR01[i] + abtr * bR01[i];
        if (i < Mr_)
          A3u[count] = s2 * dRp1[i] + abtr * bRp1[i];
      }
  }
};

void HHWPDE::solvePDE(Vdoub &Uin, Vdoub & /*Uout*/) {
  // Grid SetUp:
  setMr(30);
  setMv(30);
  setMs(30);
  setNt(100);
  double theta = 2.0 / 3.0;
  double dT = T_ / Nt_;

  // Mesh Generator:
  generateMeshes();

  Vdoub Y0(N_), Y1(N_), Y2(N_), Y3(N_);
  Vdoub Uold(N_, 0), Unew(N_, 0);
  Vdoub A0U(N_, 0), A1U(N_, 0), A2U(N_, 0), A3U(N_, 0);

  Vdoub A1l, A1d, A1u;
  Vdoub A2l1, A2l2, A2d, A2u1, A2u2;
  Vdoub A3lt, A3dt, A3ut;
  fillA1(A1l, A1d, A1u);
  fillA2(A2l1, A2l2, A2d, A2u1, A2u2);

  for (int t = 0; t <= Nt_; t++) {
    /*Step 1*/
    A0Ug(Uin, A0U);
    A1Ug(Uin, A1U);
    A2Ug(Uin, A2U);
    A3Ug(t * dT, Uin, A3U);
    aAbB(dT, A0U, dT, A1U, Y1);
    aAbB(dT, A2U, dT, A3U, Y2);
    aAbB(1, Y1, 1, Y2, Y3);
    aAbB(1, Uin, 1, Y3, Y0);

    /*Step 2*/

    aAbB(1.0, Y0, -theta * dT, A1U, Y0);
    fillA3(A3lt, A3dt, A3ut, t * dT);
  }
};
} // namespace velesquant