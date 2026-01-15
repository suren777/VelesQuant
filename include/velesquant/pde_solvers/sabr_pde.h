#include <vector>

namespace velesquant {
// using namespace std; // REMOVED by refactor

class sabr_pde {
public:
  virtual ~sabr_pde() = default;
  double getAlpha() { return alpha_; };
  double getBeta() { return beta_; };
  double getNu() { return nu_; };
  double getRho() { return rho_; };

  std::vector<double> getDensity() {
    calculateDensity();
    return Q_;
  };
  std::vector<double> getFgrid() { return Fm_; }

protected:
  double alpha_, beta_, rho_, nu_, f_, dT_, h_, T_, QL_, QR_, zmin_, zmax_;
  int sizeX_, sizeT_, j0_;
  virtual void setZbounds() = 0;
  virtual double F(double z) = 0;
  virtual double G(double z) = 0;
  virtual double C(double y, double f) = 0;
  virtual double L(double F) = 0;
  std::vector<double> Q_, Cm_, Em_, Fm_;

  double Y(double z);
  void calculateDensity();

  std::vector<double> Y(std::vector<double> z);
  std::vector<double> F(std::vector<double> z);
  std::vector<double> C(std::vector<double> y, std::vector<double> f);
  std::vector<double> G(std::vector<double> z);
  void oneStepForward(std::vector<double> &inV, double dt, double &PL,
                      double &PR);

  double Emd(double dt, double Gamma);
  std::vector<double> Emd(double dt, std::vector<double> Gamma);
  void upE(std::vector<double> &inE, std::vector<double> &upE);
};

class afSABR : public sabr_pde {
public:
  afSABR(double alpha, double beta, double nu, double rho, double shift,
         double maturity, double F, int sizeX, int sizeT, double nd)
      : shift_(shift), nd_(nd) {
    alpha_ = alpha, beta_ = beta, rho_ = rho, nu_ = nu, T_ = maturity, f_ = F,
    sizeX_ = sizeX, sizeT_ = sizeT;
  };
  ~afSABR();
  double getShift() { return shift_; };

protected:
  void setZbounds();
  double yStrike(double z);
  double Y(double z);
  double F(double z);
  double C(double y, double f);
  double G(double z);
  double L(double F);
  double Lm(double F);

private:
  double shift_;
  double nd_;
};

class antonovSABR : public sabr_pde {
public:
  antonovSABR(double alpha, double beta, double nu, double rho, double maturity,
              double F, int sizeX, int sizeT, double nd)
      : nd_(nd) {
    alpha_ = alpha, beta_ = beta, rho_ = rho, nu_ = nu, T_ = maturity, f_ = F,
    sizeX_ = sizeX, sizeT_ = sizeT;
  };
  ~antonovSABR();

protected:
  void setZbounds();
  double Y(double z);
  double F(double z);
  double C(double y, double f);
  double G(double z);
  double L(double F);
  double Lm(double F);
  double yStrike(double z);

private:
  double nd_;
};
} // namespace velesquant