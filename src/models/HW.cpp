#include <velesquant/models/hw.h>
#include <velesquant/models/utility.h>
#include <velesquant/volatility/lm.h>

#include <algorithm>
#include <boost/math/distributions.hpp>
#include <cmath>
#include <functional>
#include <ql/quantlib.hpp>

#include <velesquant/constants.h>

const double DT = velesquant::constants::numerical::HULL_WHITE_DT;
const double SDT = sqrt(DT);
using namespace std;
using namespace QuantLib;
#pragma warning(disable : 4715)

namespace velesquant {

double HullWhite::optionBond(double Expiry, double Maturity, double Strike,
                             OptionType type) {
  double v0 = totalVariance(Expiry);
  double impVol =
      (1 - exp(-kappa_ * (Maturity - Expiry))) / kappa_ * sqrt(v0 / Expiry);
  return formulaBlack(getDF(Expiry), getDF(Maturity), Strike, impVol, Expiry,
                      type);
};

double HullWhite::ZC(double Expiry) { return getDF(Expiry); }

double HullWhite::swaption(double Expiry, double Tenor, double Strike,
                           double PayFrequency) {
  double CP = criticalPoint(Expiry, Tenor, Strike, PayFrequency);
  double vO = totalVariance(Expiry);
  double dfT0 = getDF(Expiry);
  int nC = int(Tenor / PayFrequency + 0.5);
  double swaptionV = 0.0;
  for (int i = 1; i <= nC; i++) {
    double Ti = Expiry + i * PayFrequency;
    double dfTi = getDF(Ti);
    double Gi = (1 - exp(-kappa_ * (Ti - Expiry))) / kappa_;
    double Ki = dfTi / dfT0 * exp(-CP * Gi - 0.5 * vO * Gi * Gi);
    double impVol =
        (1 - exp(-kappa_ * (Ti - Expiry))) / kappa_ * sqrt(vO / Expiry);
    double Puti = formulaBlack(dfT0, dfTi, Ki, impVol, Expiry, OptionType::Put);
    swaptionV += Strike * PayFrequency * Puti;
    if (i == nC)
      swaptionV += Puti;
  }
  return swaptionV;
};

vector<double> HullWhite::simulation(vector<double> times) const {
  int N = times.size();
  vector<double> path(N);
  int timeIndex = 0;
  double variance = 0.0;
  double shortRate = 0.0;
  for (int i = 0; i <= int(times[N - 1] / DT + 0.5); i++) {
    if (i * DT <= times[timeIndex] && (i + 1) * DT > times[timeIndex]) {
      path[timeIndex] = shortRate + whichFwdRate(i * DT);
      timeIndex++;
    }
    double sigma = whichSigma(i * DT);
    variance += sigma * sigma *
                (exp(2 * kappa_ * i * DT) - exp(2 * kappa_ * (i - 1) * DT)) /
                (2 * kappa_);
    shortRate +=
        (exp(-2 * kappa_ * i * DT) * variance - kappa_ * shortRate) * DT +
        sigma * SDT * random_normal();
  }
  return path;
};

void HullWhite::calibrator(const std::vector<defSwap> &swapQuotes,
                           CalibrationTarget target) {
  quoteSwap_ = swapQuotes;
  double pos;
  timeSigmas_.resize(0);
  timeSigmas_.push_back(quoteSwap_[0].Expiry);
  pos = timeSigmas_[0];
  for (int i = 1; i < int(quoteSwap_.size()); i++) {
    if (quoteSwap_[i].Expiry > pos) {
      timeSigmas_.push_back(quoteSwap_[i].Expiry);
      pos = quoteSwap_[i].Expiry;
    }
  }
  sigmas_.resize(timeSigmas_.size());
  int n = static_cast<int>(sigmas_.size()) + 1;
  std::vector<double> x(n);
  std::vector<double> lb(n);
  std::vector<double> ub(n);
  for (int i = 0; i < n - 1; i++) {
    x[i] = sigmas_[i];
    lb[i] = 1e-6;
    ub[i] = 1;
  }
  x[n - 1] = kappa_;
  lb[n - 1] = 1e-6;
  ub[n - 1] = 1;
  int m = static_cast<int>(quoteSwap_.size());
  std::vector<double> fvec(m);
  double ftol = constants::calibration::TOLERANCE;
  double xtol = constants::calibration::TOLERANCE;
  double gtol = constants::calibration::TOLERANCE;
  int maxfev = constants::calibration::MAX_ITERATIONS;
  double epsfcn = constants::calibration::EPSILON;
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
  marketSwaption_.resize(m);
  if (target == CalibrationTarget::Price) {
    for (int i = 0; i < m; i++)
      marketSwaption_[i] = swaptionATM(
          quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].VolATM);
  } else {
    for (int i = 0; i < m; i++)
      marketSwaption_[i] = quoteSwap_[i].VolATM;
  }

  double *lb_ptr = lb.data();
  double *ub_ptr = ub.data();
  bool isSwap = (target == CalibrationTarget::Price);
  auto fcn = [this, lb_ptr, ub_ptr, isSwap](int m, int n, double *x,
                                            double *fvec, int *iflag) {
    if (isSwap)
      this->objFcnPrice(m, n, x, fvec, iflag, lb_ptr, ub_ptr);
    else
      this->objFcnIV(m, n, x, fvec, iflag, lb_ptr, ub_ptr);
  };

  lmdif(m, n, x.data(), fvec.data(), ftol, xtol, gtol, maxfev, epsfcn,
        diag.data(), mode, factor, nprint, &info, &nfev, fjac.data(), ldfjac,
        ipvt.data(), qtf.data(), wa1.data(), wa2.data(), wa3.data(), wa4.data(),
        fcn);
  // the below is output result
  for (int i = 0; i < n - 1; i++)
    sigmas_[i] = fabs(x[i]); // sigmas final value
  kappa_ = x[n - 1];         // kappa final value
  // RAII: vectors automatically cleaned up
};

double HullWhite::pen_fun(double *x, double *lb, double *ub, int n) {
  double penalty = 0.0;
  double lambda = constants::calibration::PENALTY_LAMBDA;
  for (int i = 0; i < n; i++)
    penalty += ((x[i] < lb[i]) ? lambda * pow(x[i], 2) : 0) + ((x[i] > ub[i]))
                   ? lambda * pow(x[i], 2)
                   : 0;
  return penalty;
}

void HullWhite::objFcnPrice(int m, int n, double *x, double *fvec, int *iflag,
                            double *lb, double *ub) {
  double penalty = pen_fun(x, lb, ub, n);
  for (int i = 0; i < n - 1; i++) {
    sigmas_[i] = x[i];
  }
  kappa_ = x[n - 1];
  for (int i = 0; i < m; i++) {
    double modelSwaption =
        swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                 quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
    fvec[i] = modelSwaption - marketSwaption_[i] + penalty;
  }
}

void HullWhite::objFcnIV(int m, int n, double *x, double *fvec, int *iflag,
                         double *lb, double *ub) {
  double penalty = pen_fun(x, lb, ub, n);
  for (int i = 0; i < n - 1; i++)
    sigmas_[i] = x[i];
  kappa_ = x[n - 1];
  for (int i = 0; i < m; i++) {
    double modelSwaption =
        swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                 quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
    fvec[i] = swaptionIVblack(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor,
                              modelSwaption) -
              marketSwaption_[i] + penalty;
  }
};

double HullWhite::totalVariance(double T0) {
  double var = 0.0;
  int nP = timeSigmas_.size();
  for (int n = 0; n < nP; n++) {
    double ti = 0.0;
    if (n > 0)
      ti = timeSigmas_[n - 1];
    double te = timeSigmas_[n];
    if (T0 > te)
      var += sigmas_[n] * sigmas_[n] *
             (exp(2 * kappa_ * te) - exp(2 * kappa_ * ti)) / (2 * kappa_);
    else {
      var += sigmas_[n] * sigmas_[n] *
             (exp(2 * kappa_ * T0) - exp(2 * kappa_ * ti)) / (2 * kappa_);
      break;
    }
  }

  if (T0 > timeSigmas_[nP - 1])
    var += sigmas_[nP - 1] * sigmas_[nP - 1] *
           (exp(2 * kappa_ * T0) - exp(2 * kappa_ * timeSigmas_[nP - 1])) /
           (2 * kappa_);

  return exp(-2 * kappa_ * T0) * var;
};

double HullWhite::equalityCP(double x, double vO, double Expiry, double Tenor,
                             double Strike, double PayFrequency) {
  double dfT0 = getDF(Expiry);
  int nC = int(Tenor / PayFrequency + 0.5);
  double CP = 1.0;
  for (int i = 1; i <= nC; i++) {
    double Ti = Expiry + i * PayFrequency;
    double Gi = (1 - exp(-kappa_ * (Ti - Expiry))) / kappa_;
    double Ki = getDF(Ti) / dfT0 * exp(-x * Gi - 0.5 * vO * Gi * Gi);
    CP -= Strike * PayFrequency * Ki;
    if (i == nC)
      CP -= Ki;
  }
  return CP;
};

double HullWhite::criticalPoint(double Expiry, double Tenor, double Strike,
                                double PayFrequency) {
  double totalVar = totalVariance(Expiry);
  double lowerBound = -0.5;
  double upperBound = 0.5;
  double shortRate = 0.5 * (lowerBound + upperBound);
  double equalityValue =
      equalityCP(shortRate, totalVar, Expiry, Tenor, Strike, PayFrequency);
  int Niter = 0;
  do {
    Niter++;
    if (equalityValue < 0.0)
      lowerBound = shortRate;
    if (equalityValue > 0.0)
      upperBound = shortRate;
    shortRate = 0.5 * (lowerBound + upperBound);
    equalityValue =
        equalityCP(shortRate, totalVar, Expiry, Tenor, Strike, PayFrequency);
  } while (fabs(equalityValue) >= 1.0E-10 && Niter < 20 &&
           (upperBound - lowerBound) >= 10e-6);
  return shortRate;
};

double HullWhite::formulaBlack(double dfT0, double dfTN, double Strike,
                               double impVol, double T0, OptionType type) {
  double sd = impVol * sqrt(T0);
  double d1 = log(dfTN / (dfT0 * Strike)) / sd + 0.5 * sd;
  double d2 = d1 - sd;
  double price;
  if (type == OptionType::Call)
    price = dfTN * cdf_normal(d1) - (dfT0 * Strike) * cdf_normal(d2);
  else
    price = (dfT0 * Strike) * cdf_normal(-d2) - dfTN * cdf_normal(-d1);
  return price;
};

double HullWhite::whichFwdRate(double t) const {
  if (t < timeDFs_.front())
    return -std::log(DFs_.front()) / timeDFs_.front();

  auto it = std::upper_bound(timeDFs_.begin(), timeDFs_.end(), t);

  if (it == timeDFs_.end()) {
    size_t n = timeDFs_.size();
    return -std::log(DFs_.back() / DFs_[n - 2]) /
           (timeDFs_.back() - timeDFs_[n - 2]);
  }

  size_t i = std::distance(timeDFs_.begin(), it);
  return -std::log(DFs_[i] / DFs_[i - 1]) / (timeDFs_[i] - timeDFs_[i - 1]);
}

double HullWhite::whichSigma(double t) const {
  if (t < timeSigmas_.front())
    return sigmas_.front();

  auto it = std::upper_bound(timeSigmas_.begin(), timeSigmas_.end(), t);

  if (it == timeSigmas_.end())
    return sigmas_.back();

  size_t i = std::distance(timeSigmas_.begin(), it);
  return sigmas_[i];
}

double HullWhite::getSwapRate(double Expiry, double Tenor,
                              double PayFrequency) {
  double dfT0 = getDF(Expiry);
  int nC = int(Tenor / PayFrequency + 0.5);
  double Level = 0.0;
  for (int i = 1; i <= nC; i++) {
    double dfTi = getDF(Expiry + i * PayFrequency);
    Level += PayFrequency * dfTi;
    if (i == nC)
      dfT0 -= dfTi;
  }
  return dfT0 / Level;
};

void HullWhite::calibratorBstrp(const std::vector<defSwap> &swapQuotes,
                                CalibrationTarget target) {
  quoteSwap_ = swapQuotes;

  int MAX_ITER = constants::calibration::MAX_ITER_LOW;
  double X_TOL = constants::calibration::TOL_COARSE;
  double F_TOL = constants::calibration::TOL_FINE;
  double dF, F;
  double x = sigmas0_[0];
  int k = 0;
  std::function<double(double)> objFcnB;
  double min = 1e-6;
  double max = 1;

  calibrateKappa();
  double pos;
  timeSigmas_.resize(0);
  timeSigmas_.push_back(quoteSwap_[0].Expiry);
  pos = timeSigmas_[0];
  for (int i = 1; i < int(quoteSwap_.size()); i++) {
    if (quoteSwap_[i].Expiry > pos) {
      timeSigmas_.push_back(quoteSwap_[i].Expiry);
      pos = quoteSwap_[i].Expiry;
    }
  }
  sigmas_.resize(timeSigmas_.size());
  for (int j = 0; j < int(timeSigmas_.size()); j++) {
    iter_ = j;

    qSB_.resize(0);
    for (int l = k; l < int(quoteSwap_.size()); l++, k++) {
      if (l <= int(quoteSwap_.size())) {
        if (!(quoteSwap_[l].Expiry <= timeSigmas_[j])) {
          break;
        } else
          qSB_.push_back(quoteSwap_[l]);
      }
    }
    int m = static_cast<int>(qSB_.size()); // no. of observations
    marketSwaption_.resize(m);
    if (target == CalibrationTarget::Price) {
      for (int i = 0; i < m; i++)
        marketSwaption_[i] =
            swaptionATM(qSB_[i].Expiry, qSB_[i].Tenor, qSB_[i].VolATM);

      objFcnB = [this, m](double x) { return this->objFcnBswap(x, m); };
    } else {
      for (int i = 0; i < m; i++)
        marketSwaption_[i] = qSB_[i].VolATM;

      objFcnB = [this, m](double x) { return this->objFcnBIV(x, m); };
    }
    double delta;
    for (int i = 0; i < MAX_ITER; i++) {
      F = objFcnB(x);
      double Fu = objFcnB(x + X_TOL);
      double Fd = objFcnB(x - X_TOL);
      dF = 0.5 * (Fu - Fd) / X_TOL;
      if ((fabs(F) < F_TOL) || (fabs(dF) < F_TOL))
        break;
      double x1 = x;
      delta = F / dF;
      x1 -= delta;
      if (x1 <= min) {
        delta = 0.5F * (x - min);
        x1 = x - delta;
        if ((x1 == min) || (x1 == max))
          break;
      } else if (x1 >= max) {
        delta = 0.5F * (x - max);
        x1 = x - delta;
        if ((x1 == min) || (x1 == max))
          break;
      }
      if (delta > 0)
        max = x;
      else
        min = x;
      if (abs(x1 - x) < X_TOL * fabs(x1) || (x1 != x1))
        break;
      x = x1;
      if (fabs(max - min) < 1e-7)
        break;
    }
    sigmas_[j] = x;
    min = 1e-6;
    max = 1.0;
  }
};

double HullWhite::objFcnBIV(double x, int m) {
  double sum = 0.0;
  sigmas_[iter_] = x; // sigmas
  for (int i = 0; i < m; i++) {
    double modelSwaption = swaption(qSB_[i].Expiry, qSB_[i].Tenor,
                                    qSB_[i].SwapRate, qSB_[i].Frequency);
    double a =
        1.0 - swaptionIVblack(qSB_[i].Expiry, qSB_[i].Tenor, modelSwaption) /
                  marketSwaption_[i];
    sum += fabs(a);
  }
  return sum / m * 1e4;
};
double HullWhite::objFcnBswap(double x, int m) {
  double sum = 0.0;
  sigmas_[iter_] = x; // sigmas
  for (int i = 0; i < m; i++) {
    double modelSwaption = swaption(qSB_[i].Expiry, qSB_[i].Tenor,
                                    qSB_[i].SwapRate, qSB_[i].Frequency);
    double a = modelSwaption / marketSwaption_[i] - 1.0;
    sum += a * a;
  }
  return sum / m;
};

void HullWhite::calibrateKappa() {
  int MAX_ITER = constants::calibration::MAX_ITER_HIGH;
  double X_TOL = 1e-7; // Specific tolerance
  double F_TOL = constants::calibration::TOL_FINE;
  double a0;
  double a1;
  double dF, F;
  double min = 1e-6;
  double max = 1.0;
  int N = quoteSwap_.size();
  vector<double> iv_ratio;
  vector<double> bond_ratio;
  vector<int> index;

  for (int i = 0; i < N - 1; i++) {
    if (quoteSwap_[i].Expiry == quoteSwap_[i + 1].Expiry) {
      double ex = quoteSwap_[i].Expiry;
      double t1 = quoteSwap_[i].Tenor + ex;
      double t2 = quoteSwap_[i + 1].Tenor + ex;
      index.push_back(i);
      iv_ratio.push_back(quoteSwap_[i + 1].VolATM / quoteSwap_[i].VolATM);
      bond_ratio.push_back((getDF(ex) - getDF(t1)) / (getDF(ex) - getDF(t2)));
    }
  }

  a0 = kappa0_;
  double delta;
  for (int i = 0; i < MAX_ITER; i++) {
    F = objKappa(a0, iv_ratio, bond_ratio, index);
    double Fu = objKappa(a0 + X_TOL, iv_ratio, bond_ratio, index);
    double Fd = objKappa(a0 - X_TOL, iv_ratio, bond_ratio, index);
    dF = 0.5 * (Fu - Fd) / X_TOL;
    if ((fabs(F) < F_TOL) || (fabs(dF) < F_TOL))
      break;
    a1 = a0;
    delta = F / dF;
    a1 -= delta;
    if (a1 <= min) {
      delta = 0.5F * (a0 - min);
      a1 = a0 - delta;
      if ((a1 == min) || (a1 == max))
        break;
    } else if (a1 >= max) {
      delta = 0.5F * (a0 - max);
      a1 = a0 - delta;
      if ((a1 == min) || (a1 == max))
        break;
    }
    if (delta > 0)
      max = a0;
    else
      min = a0;
    if (abs(a1 - a0) < X_TOL * fabs(a1) || (a1 != a1))
      break;
    a0 = a1;
    if (fabs(max - min) < 1e-5)
      break;
  }
  QL_ENSURE(a0 == a0, "Failed to calibrate Kappa \n");
  kappa_ = a0;
};

double HullWhite::Bratio(double a, double M, double t1, double t2) {
  return (1.0 - exp(-a * (t2 - M))) / (1.0 - exp(-a * (t1 - M)));
};

double HullWhite::objKappa(double x, const vector<double> &ivr,
                           const vector<double> &br, const vector<int> &ind) {
  int N = ivr.size();
  double sum = 0.0;
  if (x == 0)
    return 1000;
  for (int i = 0; i < N; i++) {
    double ex = quoteSwap_[ind[i]].Expiry;
    double t1 = quoteSwap_[ind[i]].Tenor + ex;
    double t2 = quoteSwap_[ind[i] + 1].Tenor + ex;
    double a = br[i] * Bratio(x, ex, t1, t2) - ivr[i];
    sum += a * a;
  }
  return sum / N;
};

double HullWhite::swaptionIVblack(double Expiry, double Tenor,
                                  double swaption_price) {
  double xl = 1e-6;
  double xu = 1;
  double x = 0.5 * (xl + xu);

  double DF = getDF(Expiry) - getDF(Expiry + Tenor);
  double b =
      DF * (2.0 * cdf_normal(0.5 * x * sqrt(Expiry)) - 1.0) - swaption_price;

  int Niter = 0;
  do {
    Niter++;
    if (b < 0.0)
      xl = x;
    if (b > 0.0)
      xu = x;
    x = 0.5 * (xl + xu);
    b = DF * (2.0 * cdf_normal(0.5 * x * sqrt(Expiry)) - 1.0) - swaption_price;
  } while ((fabs(b) >= 1.0E-12) && (Niter < 40) && (xu - xl > 1E-7));
  return x;
};

double HullWhite::BlackStrike(double T0, double TN, double impVol) {
  double xl = 1e-6;
  double xu = 1;
  double x = 0.5 * (xl + xu);
  double dfT0 = getDF(T0);
  double dfTN = getDF(T0 + TN);
  double DF = dfT0 - dfTN;
  double atm_swap = DF * (2.0 * cdf_normal(0.5 * impVol * sqrt(T0)) - 1.0);
  double b =
      formulaBlack(dfT0, dfTN, x, impVol, T0, OptionType::Call) - atm_swap;

  int Niter = 0;
  do {
    Niter++;
    if (b < 0.0)
      xl = x;
    if (b > 0.0)
      xu = x;
    x = 0.5 * (xl + xu);
    b = formulaBlack(dfT0, dfTN, x, impVol, T0, OptionType::Call) - atm_swap;
  } while ((fabs(b) >= 1.0E-12) && (Niter < 40) && (xu - xl > 1E-7));
  return x;
};
double HullWhite::BlackStrikePlain(double T0, double TN) {
  double dfT0 = getDF(T0);
  double dfTN = getDF(T0 + TN);
  return dfTN / dfT0;
};

double HullWhite::swaptionIVblackPub(double Expiry, double Tenor,
                                     double swap_price) {
  return swaptionIVblack(Expiry, Tenor, swap_price);
};
} // namespace velesquant