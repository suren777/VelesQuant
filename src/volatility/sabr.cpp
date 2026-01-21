
#include <algorithm>
#include <boost/math/distributions.hpp>
#include <cassert>
#include <cmath>
#include <ql/quantlib.hpp>
#include <velesquant/constants.h>
#include <velesquant/models/solver.h>
#include <velesquant/models/utility.h>
#include <velesquant/volatility/lm.h>
#include <velesquant/volatility/sabr.h>

// using namespace std;
#include <memory>
#include <velesquant/errors.h>
#include <velesquant/models/utility.h>

// using namespace std;
// using namespace QuantLib;

namespace velesquant {

using namespace velesquant::constants::sabr;
Sabr::Sabr(double maturity, double forward, double beta, double alpha,
           double nu, double rho)
    : maturity_(maturity), forward_(forward), beta_(beta), alpha_(alpha),
      nu_(nu), rho_(rho) {
  parameters();
  shift_ = 0;
  notQUTL_ = true;
  type_ = "BS";
};

Sabr::Sabr(double maturity, double forward, double beta, double alpha,
           double nu, double rho, double shift)
    : maturity_(maturity), forward_(forward), beta_(beta), shift_(shift),
      alpha_(alpha), nu_(nu), rho_(rho) {
  parameters();
  notQUTL_ = true;
  type_ = "Normal";
};

double Sabr::normalVol(double K) const {
  if (beta_ == 0.0)
    if (K == forward_)
      return (Nvolb0(K * ATM_BUMP_DOWN) + Nvolb0(K * ATM_BUMP_UP)) / 2;
    else
      return Nvolb0(K);
  else if (beta_ == 1.0)
    if (K == forward_)
      return (Nvolb1(K * ATM_BUMP_DOWN) + Nvolb1(K * ATM_BUMP_UP)) / 2;
    else
      return Nvolb1(K);
  else if (K == forward_)
    return (Nvolb(K * ATM_BUMP_DOWN) + Nvolb(K * ATM_BUMP_UP)) / 2;
  else
    return Nvolb(K);
}

double Sabr::Nvolb0(double K) const {
  double zeta = nu_ / alpha_ * (forward_ - K);
  double xOfZeta = std::log(
      (std::sqrt(1 - rho_ * zeta + zeta * zeta) - rho_ + zeta) / (1 - rho_));
  return alpha_ * zeta / xOfZeta *
         (1 + (COEFF_4 * rho_ * nu_ * alpha_ +
               COEFF_24 * (2 - 3 * rho_ * rho_) * nu_ * nu_) *
                  maturity_);
};

double Sabr::Nvolb(double K) const {
  double forwardBetaDiff =
      std::pow(forward_ + shift_, beta_) - std::pow(K + shift_, beta_);
  double forwardOneMinusBetaDiff =
      std::pow(forward_ + shift_, 1 - beta_) - std::pow(K + shift_, 1 - beta_);
  double logForwardOverStrike = std::log((forward_ + shift_) / (K + shift_));
  double zeta = nu_ / alpha_ * forwardOneMinusBetaDiff / (1 - beta_);
  double xOfZeta = std::log(
      (std::sqrt(1 - rho_ * zeta + zeta * zeta) - rho_ + zeta) / (1 - rho_));
  double baseTerm = alpha_ * (forward_ - K) * (1 - beta_) /
                    forwardOneMinusBetaDiff * zeta / xOfZeta;
  double betaCorrection =
      -1.0 / 24 *
      (beta_ * (2 - beta_) * (1 - beta_) * (1 - beta_) * alpha_ * alpha_ *
       logForwardOverStrike * logForwardOverStrike) /
      (forwardOneMinusBetaDiff * forwardOneMinusBetaDiff);
  double volOfVolCorrection =
      0.25 * rho_ * nu_ * alpha_ * forwardBetaDiff / (forward_ - K) +
      (2 - 3 * rho_ * rho_) / 24 * nu_ * nu_;
  return baseTerm * (1 + (betaCorrection + volOfVolCorrection) * maturity_);
};

double Sabr::Nvolb1(double K) const {
  double logForwardPlusLogStrike =
      std::log(forward_ + shift_) + std::log(K + shift_);
  double zeta = nu_ / alpha_ * logForwardPlusLogStrike;
  double xOfZeta = std::log(
      (std::sqrt(1 - rho_ * zeta + zeta * zeta) - rho_ + zeta) / (1 - rho_));
  double baseTerm =
      alpha_ * (forward_ - K) / logForwardPlusLogStrike * zeta / xOfZeta;
  double volOfVolCorrection = 0.25 * rho_ * nu_ * alpha_ +
                              (2 - 3 * rho_ * rho_) / 24 * nu_ * nu_ -
                              1.0 / 24 * alpha_ * alpha_;
  return baseTerm * (1.0 + volOfVolCorrection * maturity_);
};

double Sabr::getVol(double strike) const {
  if (type_ == "BS")
    return impliedVol(strike);
  else
    return normalVol(strike);
};
double Sabr::getPremium(double strike, OptionType callORput) const {
  if (type_ == "BS")
    return premiumBlackScholes(strike, callORput);
  else
    return premiumBachelier(strike, callORput);
};
double Sabr::impliedVol(double strike) const {
  if (!calibrated_)
    VEL_RAISE("SABR model is not calibrated yet!");
  double geometricMean = std::pow(forward_ * strike, 0.5 * (1.0 - beta_));
  double logMoneyness = std::log(forward_ / strike);
  double moneynessCorrectionSq =
      (1.0 - beta_) * (1.0 - beta_) * logMoneyness * logMoneyness;
  double denominator =
      geometricMean * (1.0 + moneynessCorrectionSq / 24.0 +
                       moneynessCorrectionSq * moneynessCorrectionSq / 1920.0);
  double timeCorrection = (1.0 - beta_) * alpha_ / geometricMean *
                              (1.0 - beta_) * alpha_ / geometricMean / 24.0 +
                          rho_ * beta_ * nu_ * alpha_ / geometricMean / 4.0 +
                          (2.0 - 3.0 * rho_ * rho_) * nu_ * nu_ / 24.0;
  double zeta = nu_ / alpha_ * geometricMean * logMoneyness;
  double multiplier;
  if (std::fabs(zeta * zeta) > 1.0E-20) {
    double xOfZeta = std::log(
        (std::sqrt(1.0 - 2.0 * rho_ * zeta + zeta * zeta) + zeta - rho_) /
        (1.0 - rho_));
    multiplier = zeta / xOfZeta;
  } else
    multiplier = 1.0 - 0.5 * rho_ * zeta +
                 (2.0 - 3.0 * rho_ * rho_) * zeta * zeta / 12.0;
  double impVol =
      alpha_ / denominator * multiplier * (1.0 + maturity_ * timeCorrection);
  if (impVol < 0.0)
    impVol = 1.0E-09;
  return impVol;
}

double Sabr::premiumBachelier(double strike, OptionType callORput) const {
  int I = (callORput == OptionType::Call) ? 1 : -1;
  double vol = normalVol(strike) * std::sqrt(maturity_);
  double q = (forward_ - strike) / vol;
  return I * (forward_ - strike) * cdf_normal(I * q) + vol * pdf_normal(I * q);
};

double Sabr::premiumBlackScholes(double strike, OptionType callORput) const {
  double vol = impliedVol(strike) * std::sqrt(maturity_);
  double d1 = std::log(forward_ / strike) / vol + 0.5 * vol;
  double d2 = d1 - vol;
  double premium;
  if (forward_ >= strike) {
    if ((d1 != d1) || (d2 != d2))
      premium = forward_ - strike;
    else
      premium = forward_ * cdf_normal(d1) - strike * cdf_normal(d2);
    if (callORput == OptionType::Put)
      premium -= (forward_ - strike);
  } else {
    if ((d1 != d1) || (d2 != d2))
      premium = 0.0;
    else
      premium = strike * cdf_normal(-d2) - forward_ * cdf_normal(-d1);
    // Uses OptionType check
    if (callORput != OptionType::Put)
      premium += (forward_ - strike);
  }
  return premium;
};

double Sabr::localVol(double spot) const { return localVolCall(spot); };

double Sabr::localVolCall(double spot) const // based on the call premium
{
  double cLft = premiumBlackScholes(spot * SPOT_BUMP_DOWN);
  double cVal = premiumBlackScholes(spot);
  double cRgt = premiumBlackScholes(spot * SPOT_BUMP_UP);
  double dS2val =
      (cRgt + cLft - 2.0 * cVal) / (constants::numerical::FD_BUMP_SIZE * spot *
                                    constants::numerical::FD_BUMP_SIZE * spot);
  if (dS2val < 1.0E-10)
    dS2val = 1.0E-10;
  assert(dS2val > 0.0);
  assert(dS2val > 0.0);
  auto incrT = std::make_unique<Sabr>(maturity_ * SPOT_BUMP_UP, forward_, beta_,
                                      alpha_, nu_, rho_);
  double cBigT = incrT->premiumBlackScholes(spot);
  double dTval =
      (cBigT - cVal) / (constants::numerical::FD_BUMP_SIZE * maturity_);
  if (dTval < 1.0E-10)
    dTval = 1.0E-10;
  assert(dTval > 0.0);
  double lv = std::sqrt(2.0 * dTval / dS2val) / spot;
  return lv;
};

double Sabr::localVolIV(double spot) const // based on the implied vol
{
  double volL = impliedVol(spot * SPOT_BUMP_DOWN);
  double vol = impliedVol(spot);
  double volR = impliedVol(spot * SPOT_BUMP_UP);
  double dSvol =
      (volR - volL) / (2.0 * constants::numerical::FD_BUMP_SIZE * spot);
  double dS2vol =
      (volR + volL - 2.0 * vol) / (constants::numerical::FD_BUMP_SIZE * spot *
                                   constants::numerical::FD_BUMP_SIZE * spot);
  auto incrT = std::make_unique<Sabr>(maturity_ * SPOT_BUMP_UP, forward_, beta_,
                                      alpha_, nu_, rho_);
  double volT = incrT->impliedVol(spot);
  double dTvol =
      (volT - vol) / (constants::numerical::FD_BUMP_SIZE * maturity_);
  double derivation = vol * sqrt(maturity_);
  double d1 = log(forward_ / spot) / derivation + 0.5 * derivation;
  double d2 = d1 - derivation;
  double fz = vol * vol + 2.0 * vol * maturity_ * dTvol;
  double fm =
      1.0 + 2.0 * d1 * spot * sqrt(maturity_) * dSvol +
      spot * spot * maturity_ * (d1 * d2 * dSvol * dSvol + vol * dS2vol);
  double lv = sqrt(fz / fm);
  return lv;
};

double Sabr::localVolzabr(double spot) const {
  double y =
      nu_ / alpha_ *
      (std::pow(forward_, (1.0 - beta_)) - std::pow(spot, (1.0 - beta_))) /
      (1.0 - beta_);
  double jy = std::sqrt(1.0 - 2.0 * rho_ * y + y * y);
  double lv = jy * alpha_ * std::pow(spot, beta_) / spot;
  return lv;
};

void Sabr::objFcn(int m, int /*n*/, double *x, double *fvec, int * /*iflag*/,
                  CalibrationTarget quoteType) {
  setParameterAlpha(x[0]);
  setParameterNu(x[1]);
  setParameterRho(x[2]);
  double penalty =
      ((alpha_ <= 0) ? constants::calibration::PENALTY_HUGE * alpha_ * alpha_
                     : 0) +
      ((nu_ <= 0) ? constants::calibration::PENALTY_HUGE * nu_ * nu_ : 0);
  for (int i = 0; i < m; i++) {
    if (quoteType == CalibrationTarget::Price)
      fvec[i] = (getPremium(strikes_[i]) - marketQuotes_[i]) + penalty;
    else
      fvec[i] = (getVol(strikes_[i]) - marketQuotes_[i]) + penalty;
  }
};

void Sabr::objFcnATM(int m, int /*n*/, double *x, double *fvec, int * /*iflag*/,
                     CalibrationTarget quoteType) {

  setParameterNu(x[0]);
  setParameterRho(x[1]);

  if (beta_ == 1) {
    setParameterAlpha(ATMvolRoots(0.0, 1.0, 0.0001));
  } else {
    setParameterAlpha(ATMvolRoots(0.0, 100.00, 0.0001));
  }
  double penalty =
      ((nu_ <= 0) ? constants::calibration::PENALTY_LARGE * nu_ * nu_ : 0);
  for (int i = 0; i < m; i++) {
    if (quoteType == CalibrationTarget::Price)
      fvec[i] = (premiumBlackScholes(strikes_[i]) - marketQuotes_[i]) + penalty;
    else
      fvec[i] = (impliedVol(strikes_[i]) - marketQuotes_[i]) + penalty;
  }
};

void Sabr::calibratorWithInitialATM(std::vector<double> strikes,
                                    std::vector<double> marketQuotes,
                                    CalibrationTarget quoteType) {
  int m = strikes.size();
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  int n = 2;
  std::vector<double> x(n);

  x[0] = nu_;
  x[1] = rho_;

  std::vector<double> fvec(m);
  double ftol = 1e-10;
  double xtol = 1e-10;
  double gtol = 1e-10;
  int maxfev = 5000;
  double epsfcn = 1e-10;
  std::vector<double> diag(n);
  int mode = 1;
  double factor = 1;
  int nprint = 0;
  int info = 0;
  int nfev = 0;
  std::vector<double> fjac(m * n);
  int ldfjac = m;
  std::vector<int> ipvt(n);
  std::vector<double> qtf(n);
  std::vector<double> wa1(n);
  std::vector<double> wa2(n);
  std::vector<double> wa3(n);
  std::vector<double> wa4(m);

  lmfcn fcn = [this, quoteType](int m, int n, double *x, double *fvec,
                                int *iflag) {
    this->objFcnATM(m, n, x, fvec, iflag, quoteType);
  };

  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcn);

  QL_ENSURE(info >= 1 && info <= 4,
            "SABR Calibration Fails: " << getLmdifMessage(info)
                                       << " (info=" << info << ")");

  // the below is output result
  setParameterNu(x[0]);
  setParameterRho(x[1]);
  if (beta_ == 1) {
    setParameterAlpha(ATMvolRoots(0.0, 1.0, 0.0001));
  } else {
    setParameterAlpha(ATMvolRoots(0.0, 100.00, 0.0001));
  }

  parameters();
};

void Sabr::calibratorWithInitial(std::vector<double> strikes,
                                 std::vector<double> marketQuotes,
                                 CalibrationTarget quoteType) {
  int m = strikes.size();
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  int n = 3;
  std::vector<double> x(n);

  x[0] = alpha_;
  x[1] = nu_;
  x[2] = rho_;

  std::vector<double> fvec(m);
  double ftol = 1e-10;
  double xtol = 1e-10;
  double gtol = 1e-10;
  int maxfev = 5000;
  double epsfcn = 1e-10;
  std::vector<double> diag(n);
  int mode = 1;
  double factor = 1;
  int nprint = 0;
  int info = 0;
  int nfev = 0;
  std::vector<double> fjac(m * n);
  int ldfjac = m;
  std::vector<int> ipvt(n);
  std::vector<double> qtf(n);
  std::vector<double> wa1(n);
  std::vector<double> wa2(n);
  std::vector<double> wa3(n);
  std::vector<double> wa4(m);

  lmfcn fcn = [this, quoteType](int m, int n, double *x, double *fvec,
                                int *iflag) {
    this->objFcn(m, n, x, fvec, iflag, quoteType);
  };

  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcn);

  QL_ENSURE(info != 4, "SABR caliration fails with the info = " << info);

  setParameterAlpha(x[0]);
  setParameterNu(x[1]);
  setParameterRho(x[2]);
  parameters();
};

void Sabr::calibrator(std::vector<double> strikes,
                      std::vector<double> marketQuotes,
                      CalibrationTarget quoteType) {
  int m = strikes.size();
  strikes_.resize(m);
  strikes_ = strikes;
  marketQuotes_.resize(m);
  marketQuotes_ = marketQuotes;
  int n = 3;
  std::vector<double> x(n);
  x[0] = alpha_;
  x[1] = nu_;
  x[2] = rho_;

  std::vector<double> fvec(m);
  double ftol = 1e-10;
  double xtol = 1e-10;
  double gtol = 1e-10;
  int maxfev = 5000;
  double epsfcn = 1e-10;
  std::vector<double> diag(n);
  int mode = 1;
  double factor = 1;
  int nprint = 0;
  int info = 0;
  int nfev = 0;
  std::vector<double> fjac(m * n);
  int ldfjac = m;
  std::vector<int> ipvt(n);
  std::vector<double> qtf(n);
  std::vector<double> wa1(n);
  std::vector<double> wa2(n);
  std::vector<double> wa3(n);
  std::vector<double> wa4(m);

  lmfcn fcn = [this, quoteType](int m, int n, double *x, double *fvec,
                                int *iflag) {
    this->objFcn(m, n, x, fvec, iflag, quoteType);
  };

  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcn);

  QL_ENSURE(info != 4, "SABR caliration fails with the info = " << info);

  setParameterAlpha(x[0]);
  setParameterNu(x[1]);
  setParameterRho(x[2]);
  parameters();
};

boost::math::normal_distribution<> stdnormal(0.0, 1.0);
void Sabr::qutlTable() {
  int N = 1200;
  std::vector<double> spots(N), qutls, cdfss;
  if (forward_ < 1.0)
    spots = grMesh(N);
  else
    spots = gsMesh(N);
  std::vector<double> isp;

  for (int i = 0; i < N; i++) {
    double cVal = getPremium(spots[i]);
    double cRgt = getPremium(spots[i] * 1.001);
    double cdf = 1.0 + 1000.0 * (cRgt - cVal) / spots[i];
    if ((cdf > 0.0) && (cdf < 1.0)) {
      qutls.push_back(boost::math::quantile(stdnormal, cdf));
      cdfss.push_back(cdf);
      isp.push_back(spots[i]);
    }
  }
  // sorting the cdf distribution in an increasing order
  std::vector<double> scdfs;
  asort(qutls, isp, cdfss);
  double current = qutls[0];
  int k = static_cast<int>(isp.size());
  for (int i = 0; i < k; i++)
    if ((qutls[i] > current)) {
      current = qutls[i];
      spots_.push_back(isp[i]);
      qutls_.push_back(qutls[i]);
      cdfs_.push_back(cdfss[i]);
    }
  // adjusting (rebiaseing) the cdf to meet the forward
  N = static_cast<int>(spots_.size());
  double mean = 0.0;
  for (int i = 0; i < N; i++) {
    if (i == 0)
      mean += 0.5 * spots_[i] * cdfs_[i];
    else
      mean += 0.5 * (spots_[i] + spots_[i - 1]) * (cdfs_[i] - cdfs_[i - 1]);
  }
  double bias = forward_ - mean;
  for (auto &spot : spots_)
    spot += bias;

  notQUTL_ = false;
};

double Sabr::simulation(double corrRN) {
  if (notQUTL_)
    qutlTable();
  if (qutls_.empty())
    VEL_RAISE("SABR calibration failed: qutls table is empty");
  QuantLib::LinearInterpolation interp(qutls_.begin(), qutls_.end(),
                                       spots_.begin());
  double qmin = qutls_[0];
  double qmax = qutls_[qutls_.size() - 1];
  bool allowExtrapolation = true;
  if (corrRN > qmax)
    corrRN = qmax;
  else if (corrRN < qmin)
    corrRN = qmin;
  return interp((qmin > corrRN) ? qmin : corrRN, allowExtrapolation);
};
std::vector<double> Sabr::simulations(std::vector<double> correlatedRNs) {
  if (notQUTL_)
    qutlTable();
  if (qutls_.empty())
    VEL_RAISE("SABR calibration failed: qutls table is empty");
  QuantLib::LogLinearInterpolation interp(qutls_.begin(), qutls_.end(),
                                          spots_.begin());
  int N = correlatedRNs.size();
  std::vector<double> spots(N);
  bool allowExtrapolation = true;
  double qmin = qutls_[0];
  double qmax = qutls_[qutls_.size() - 1];
  for (int i = 0; i < N; i++)
    spots[i] =
        interp((((qmin > correlatedRNs[i]) ? qmin : correlatedRNs[i]) > qmax)
                   ? qmax
                   : correlatedRNs[i],
               allowExtrapolation);
  return spots;
};

int Sabr::amin(const std::vector<double> &values, int position) {
  auto it = std::min_element(values.begin() + position, values.end());
  return static_cast<int>(std::distance(values.begin(), it));
};

int Sabr::amax(const std::vector<double> &values, int position) {
  auto it = std::max_element(values.begin() + position, values.end());
  return static_cast<int>(std::distance(values.begin(), it));
};

void Sabr::asort(std::vector<double> &a, std::vector<double> &b,
                 std::vector<double> &c) {
  int N = a.size();
  int imin = amin(a, 0);

  for (int i = 1; i < N; i++) {
    double auxa = a[imin];
    double auxb = b[imin];
    double auxc = c[imin];
    a[imin] = a[i - 1];
    b[imin] = b[i - 1];
    c[imin] = c[i - 1];
    a[i - 1] = auxa;
    b[i - 1] = auxb;
    c[i - 1] = auxc;
    imin = amin(a, i);
  }
};

std::vector<double> Sabr::gsMesh(int meshSize) {
  std::vector<double> Mesh(meshSize);
  double spotLeft = .7 * forward_;
  double spotRight = 1.4 * forward_;
  double meshSpacing = forward_ / 20.0;
  double spotMax = forward_ * 14.0;
  double zMin = std::asinh(-spotLeft / meshSpacing);
  double zIntermediate = (spotRight - spotLeft) / meshSpacing;
  double zMax = zIntermediate + std::asinh((spotMax - spotRight) / meshSpacing);
  double deltaZ = (zMax - zMin) / meshSize;
  for (int i = 0; i < meshSize; i++) {
    double z = zMin + i * deltaZ;
    if (z < 0)
      Mesh[i] = spotLeft + meshSpacing * std::sinh(z);
    else
      Mesh[i] = (z <= zIntermediate)
                    ? (spotLeft + meshSpacing * z)
                    : (spotRight + meshSpacing * std::sinh(z - zIntermediate));
  }
  return Mesh;
};

std::vector<double> Sabr::grMesh(int meshSize) {
  std::vector<double> Mesh(meshSize);
  double rateMax = 1.0;
  double meshSpacing = rateMax / (meshSize * .9);
  double start = std::asinh((-rateMax - forward_) / meshSpacing);
  double deltaX =
      1.0 / meshSize * (std::asinh((rateMax - forward_) / meshSpacing) - start);
  for (int i = 0; i < meshSize; i++)
    Mesh[i] = forward_ + meshSpacing * std::sinh(start + i * deltaX);
  return Mesh;
};

double Sabr::ATMvolPoly(double alpha) {
  double a = (std::pow(1 - beta_, 2) * maturity_) /
             (24 * std::pow(forward_, 2 - 2 * beta_));
  double b =
      (rho_ * beta_ * nu_ * maturity_) / (4 * std::pow(forward_, 1 - beta_));
  double c = (1 + (2 - 3 * std::pow(rho_, 2)) / (24) * nu_ * nu_ * maturity_);
  double d = atmVol_ * std::pow(forward_, 1 - beta_);

  return a * std::pow(alpha, 3) + b * std::pow(alpha, 2) + c * alpha - d; //=0
};

double Sabr::ATMvolRoots(double lBound, double uBound, double tol) {
  // Bisection algorithm
  double alpha = 0.5 * (lBound + uBound);
  double y = ATMvolPoly(alpha);
  int Step = 0;
  do {
    if (y < 0.0)
      lBound = alpha;
    if (y > 0.0)
      uBound = alpha;
    alpha = 0.5 * (lBound + uBound);
    y = ATMvolPoly(alpha);
    if ((uBound - lBound) < tol / 10)
      Step++;
  } while ((std::fabs(y - 0.0) > tol) && (Step < 50));
  return alpha;
}
} // namespace velesquant