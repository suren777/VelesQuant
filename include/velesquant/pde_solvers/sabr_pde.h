#ifndef SABR_PDE_H
#define SABR_PDE_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/cyclic_reduction.h>

namespace velesquant {

template <typename ModelType> class SabrPDE {
public:
  SabrPDE(std::shared_ptr<ModelType> model, int sizeX, int sizeT)
      : model_(model), sizeX_(sizeX), sizeT_(sizeT) {
    if (model_) {
      alpha_ = model_->getParameterAlpha();
      beta_ = model_->getBeta();
      nu_ = model_->getParameterNu();
      rho_ = model_->getParameterRho();
      f_ = model_->getForward();
      T_ = model_->getMaturity();
    }
  };
  virtual ~SabrPDE() = default;

  double getAlpha() { return alpha_; };
  double getBeta() { return beta_; };
  double getNu() { return nu_; };
  double getRho() { return rho_; };

  std::vector<double> getDensity() {
    calculateDensity();
    std::vector<double> densityF(sizeX_, 0.0);
    if (sizeX_ < 2)
      return densityF;

    // Jacobian transformation: P(F) = P(z) * |dz/dF|
    // We approximate P(F_i) * deltaF_i = P(z_i) * h_
    // => P(F_i) = Q_[i] * h_ / deltaF_i
    // where deltaF_i approx (F_{i+1} - F_{i-1}) / 2

    for (int i = 1; i < sizeX_ - 1; i++) {
      double deltaF = (Fm_[i + 1] - Fm_[i - 1]) / 2.0;
      if (std::abs(deltaF) > 1e-12) {
        densityF[i] = Q_[i] * h_ / deltaF;
      }
    }

    // Boundaries (one-sided)
    double deltaF_start = Fm_[1] - Fm_[0];
    if (std::abs(deltaF_start) > 1e-12)
      densityF[0] = Q_[0] * h_ / deltaF_start;

    double deltaF_end = Fm_[sizeX_ - 1] - Fm_[sizeX_ - 2];
    if (std::abs(deltaF_end) > 1e-12)
      densityF[sizeX_ - 1] = Q_[sizeX_ - 1] * h_ / deltaF_end;

    return densityF;
  };
  std::vector<double> getFgrid() { return Fm_; }

protected:
  std::shared_ptr<ModelType> model_;
  double alpha_, beta_, rho_, nu_, f_, dT_, h_, T_, QL_, QR_, zmin_, zmax_;
  int sizeX_, sizeT_, j0_;
  std::vector<double> Q_, Cm_, Em_, Fm_;
  mutable std::vector<double> solveL_, solveC_, solveU_, solveTempC_,
      solveTempD_;

  virtual void setZbounds() = 0;
  virtual double F(double z) = 0;
  virtual double G(double z) = 0;
  virtual double C(double y, double f) = 0;
  virtual double L(double F) = 0;

  double Y(double z) {
    return alpha_ / nu_ *
           (std::sinh(nu_ * z) + rho_ * (std::cosh(nu_ * z) - 1));
  }

  std::vector<double> Y(std::vector<double> z) {
    std::vector<double> y;
    y.reserve(sizeX_);
    for (int i = 0; i < sizeX_; i++)
      y.push_back(Y(z[i]));
    return y;
  }

  std::vector<double> F(std::vector<double> z) {
    std::vector<double> y;
    y.reserve(sizeX_);
    for (int i = 0; i < sizeX_; i++)
      y.push_back(F(z[i]));
    return y;
  }

  std::vector<double> C(std::vector<double> y, std::vector<double> f) {
    std::vector<double> aux;
    aux.reserve(sizeX_);
    for (int i = 0; i < sizeX_; i++)
      aux.push_back(C(y[i], f[i]));
    return aux;
  }

  std::vector<double> G(std::vector<double> z) {
    std::vector<double> y;
    y.reserve(sizeX_);
    for (int i = 0; i < sizeX_; i++)
      y.push_back(G(z[i]));
    return y;
  }

  void calculateDensity() {
    solveL_.resize(sizeX_);
    solveC_.resize(sizeX_);
    solveU_.resize(sizeX_);
    solveTempC_.resize(sizeX_);
    solveTempD_.resize(sizeX_);

    setZbounds();
    int J = sizeX_ - 2;
    double h0 = (zmax_ - zmin_) / J;
    j0_ = int(-zmin_ / h0);
    h_ = -zmin_ / (j0_ - .5);
    std::vector<double> z(sizeX_, 0);
    zmax_ = (J + 1) * h_ + zmin_;
    std::vector<double> ym(sizeX_, 0);
    Fm_.resize(sizeX_);

    for (int i = 0; i < sizeX_; i++) {
      double zi = i * h_ + zmin_ - 0.5 * h_;
      z[i] = zi;
      ym[i] = Y(zi);
      Fm_[i] = F(zi);
    }

    Cm_ = C(ym, Fm_);
    Cm_[0] = Cm_[1];
    Cm_[J + 1] = Cm_[J];
    std::vector<double> Gamma = G(Fm_);
    dT_ = T_ / sizeT_;
    double b = 1.0 - std::sqrt(2.0) / 2.0;
    double dt1 = dT_ * b, dt2 = dT_ * (1 - 2 * b);
    Em_.resize(sizeX_);
    for (int i = 0; i < sizeX_; i++)
      Em_[i] = 1.0;
    std::vector<double> Emdt1 = Emd(dt1, Gamma);
    Emdt1[0] = Emdt1[1];
    Emdt1[sizeX_ - 1] = Emdt1[sizeX_ - 2];
    std::vector<double> Emdt2 = Emd(dt2, Gamma);
    Emdt2[0] = Emdt2[1];
    Emdt2[sizeX_ - 1] = Emdt2[sizeX_ - 2];
    QL_ = QR_ = 0;
    double PL1, PL2, PR1, PR2;
    PL1 = PL2 = PR1 = PR2 = 0;
    std::vector<double> inV(sizeX_, 0.0);
    std::vector<double> inV2;
    // Ensure j0_ is within bounds
    if (j0_ >= 0 && j0_ < sizeX_)
      inV[j0_] = 1.0 / h_;

    for (int i = 1; i <= sizeT_; i++) {
      upE(Em_, Emdt1);
      oneStepForward(inV, dt1, PL1, PR1);
      inV2 = inV;
      PL2 = PL1, PR2 = PR1;
      upE(Em_, Emdt1);
      oneStepForward(inV2, dt1, PL2, PR2);
      for (int j = 0; j < sizeX_; j++)
        inV[j] = (std::sqrt(2.0) + 1) * inV2[j] - std::sqrt(2.0) * inV[j];
      PL1 = (std::sqrt(2.0) + 1) * PL2 - std::sqrt(2.0) * PL1;
      PR1 = (std::sqrt(2.0) + 1) * PR2 - std::sqrt(2.0) * PR1;
      upE(Em_, Emdt2);
    }

    QL_ = PL1, QR_ = PR1;
    Q_ = inV;
  }

  void oneStepForward(std::vector<double> &inV, double dt, double &PL,
                      double &PR) {
    // std::vector<double> l(sizeX_), c(sizeX_), u(sizeX_); // Optimized to use
    // members

    double ndt = dt / (2 * h_);
    // #pragma omp parallel for

    for (int i = 0; i < sizeX_; i++) {
      solveU_[i] = (i == sizeX_ - 1)
                       ? 0
                       : -ndt * Cm_[i + 1] / (Fm_[i + 1] - Fm_[i]) * Em_[i + 1];
      solveL_[i] =
          (i == 0) ? 0 : -ndt * Cm_[i - 1] / (Fm_[i] - Fm_[i - 1]) * Em_[i - 1];
      if (i == 0)
        solveC_[i] = Cm_[0] / (Fm_[1] - Fm_[0]) * Em_[0];
      else if (i == sizeX_ - 1)
        solveC_[i] = Cm_[i] / (Fm_[i] - Fm_[i - 1]) * Em_[i];
      else
        solveC_[i] =
            1 + ndt * Em_[i] * Cm_[i] *
                    (1.0 / (Fm_[i + 1] - Fm_[i]) + 1.0 / (Fm_[i] - Fm_[i - 1]));
    }
    solveU_[0] = Cm_[1] / (Fm_[1] - Fm_[0]) * Em_[1];
    solveL_[sizeX_ - 1] =
        Cm_[sizeX_ - 2] / (Fm_[sizeX_ - 1] - Fm_[sizeX_ - 2]) * Em_[sizeX_ - 1];

    TriDiagonalSolve(sizeX_, solveL_, solveC_, solveU_, inV, inV, solveTempC_,
                     solveTempD_);

    inV[0] = inV[sizeX_ - 1] = 0;
    PL += dt * Cm_[1] / (Fm_[1] - Fm_[0]) * Em_[1] * inV[1];
    PR += dt * Cm_[sizeX_ - 2] / (Fm_[sizeX_ - 2] - Fm_[sizeX_ - 1]) *
          Em_[sizeX_ - 2] * inV[sizeX_ - 2];
  }

  void upE(std::vector<double> &inE, std::vector<double> &upE) {
    for (int i = 0; i < sizeX_; i++)
      inE[i] *= upE[i];
  }

  double Emd(double dt, double Gamma) {
    return std::exp(rho_ * nu_ * alpha_ * Gamma * dt);
  }

  std::vector<double> Emd(double dt, std::vector<double> Gamma) {
    std::vector<double> aux(sizeX_, 0);
    for (int i = 0; i < sizeX_; i++)
      aux[i] = Emd(dt, Gamma[i]);
    return aux;
  }
};

template <typename ModelType> class AfSABR : public SabrPDE<ModelType> {
public:
  AfSABR(std::shared_ptr<ModelType> model, int sizeX, int sizeT, double nd)
      : SabrPDE<ModelType>(model, sizeX, sizeT), nd_(nd) {
    if (this->model_) {
      shift_ = this->model_->getShift();
    } else {
      shift_ = 0.0;
    }
  };
  ~AfSABR() = default;

  double getShift() { return shift_; };

protected:
  double shift_;
  double nd_;

  double sgn(double x) { return (x > 0) - (x < 0); }

  double L(double F) override { return std::pow(F + shift_, this->beta_); }
  double Lm(double F) { return std::pow(F + shift_, 1.0 - this->beta_); }

  double F(double z) override {
    return std::pow(Lm(this->f_) + (1.0 - this->beta_) * z,
                    1.0 / (1.0 - this->beta_));
  }

  double C(double y, double f) override {
    return std::sqrt(this->alpha_ * this->alpha_ +
                     2 * this->rho_ * this->alpha_ * this->nu_ * y +
                     this->nu_ * this->nu_ * y * y) *
           L(f);
  }

  double G(double F) override {
    return (F == this->f_) ? this->beta_ / Lm(this->f_)
                           : (L(F) - L(this->f_)) / (F - this->f_);
  }

  void setZbounds() override {
    this->zmax_ = this->zmin_ = -nd_ * std::sqrt(this->T_);
    this->zmax_ *= -1;
    if (this->beta_ < 1) {
      double aux = -Lm(this->f_) / (1.0 - this->beta_);
      double val =
          (std::sqrt(1.0 - this->rho_ * this->rho_ +
                     std::pow(this->rho_ + this->nu_ * aux / this->alpha_, 2)) -
           this->rho_ - this->nu_ * aux / this->alpha_) /
          (1.0 - this->rho_);
      // Ensure val > 0 for log
      double zb = -1.0 / this->nu_ * std::log(val);
      this->zmin_ = (zb > this->zmin_) ? zb : this->zmin_;
    }
  }

  double yStrike(double z) {
    return (sgn(z) * Lm(z) - sgn(this->f_) * Lm(this->f_)) /
           (1.0 - this->beta_);
  }
};

template <typename ModelType> class AntonovSABR : public SabrPDE<ModelType> {
public:
  AntonovSABR(std::shared_ptr<ModelType> model, int sizeX, int sizeT, double nd)
      : SabrPDE<ModelType>(model, sizeX, sizeT), nd_(nd) {};
  ~AntonovSABR() = default;

protected:
  double nd_;

  double sgn(double x) { return (x > 0) - (x < 0); }

  double L(double F) override { return std::pow(std::abs(F), this->beta_); };
  double Lm(double F) { return std::pow(std::abs(F), 1.0 - this->beta_); };

  void setZbounds() override {
    this->zmax_ = this->zmin_ = -nd_ * std::sqrt(this->T_);
    this->zmax_ *= -1;
  };

  double yStrike(double z) {
    return (sgn(z) * Lm(z) - sgn(this->f_) * Lm(this->f_)) /
           (1.0 - this->beta_);
  };

  double F(double z) override {
    if (std::abs(this->beta_ - 1.0) < 1e-12) {
      return this->f_ * std::exp(z);
    }
    double u = sgn(this->f_) * Lm(this->f_) + (1.0 - this->beta_) * z;
    return sgn(u) * std::pow(std::abs(u), 1.0 / (1.0 - this->beta_));
  };

  double C(double y, double f) override {
    return std::sqrt(this->alpha_ * this->alpha_ +
                     2 * this->rho_ * this->alpha_ * this->nu_ * y +
                     this->nu_ * this->nu_ * y * y) *
           L(f);
  };

  double G(double z) override {
    return (z == this->f_) ? sgn(this->f_) * this->beta_ / Lm(this->f_)
                           : (L(z) - L(this->f_)) / (z - this->f_);
  };
};

} // namespace velesquant

#endif // SABR_PDE_H