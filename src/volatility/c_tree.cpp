#include <algorithm>
#include <cmath>
#include <velesquant/errors.h>
#include <velesquant/volatility/c_tree.h>

using namespace std;

namespace velesquant {

CTree::CTree() : spot_(0.0) {}

CTree::~CTree() {}

CTree::CTree(double S, Vdoub T, Vdoub F, Vdoub IV, Vdoub r, Vdoub q)
    : spot_(S), T_(std::move(T)), F_(std::move(F)), IV_(std::move(IV)),
      r_(std::move(r)), q_(std::move(q)) {};

double CTree::calculateBinomial(double strike, double Maturity, int Nnodes,
                                exStyle style, OptionType pay, tType /*tree*/) {
  double dt = Maturity / double(Nnodes);
  if (r_.empty())
    r_.push_back(0.0);
  if (q_.empty())
    q_.push_back(0.0);

  Mdoub myTree;
  Vdoub optionValues(Nnodes + 1, 0);
  Vdoub upFactors(Nnodes), upProbabilities(Nnodes);

  for (int i = 0; i < Nnodes; i++) {
    upFactors[i] = exp(interp_vol(T_, IV_, i * dt, (i + 1) * dt) * sqrt(dt));
    double downFactor = 1.0 / upFactors[i];
    upProbabilities[i] =
        (exp((r_[0] - q_[0]) * dt) - downFactor) / (upFactors[i] - downFactor);
  }

  double type = (pay == OptionType::Call) ? -1.0 : 1.0;
  for (int i = 0; i <= Nnodes; i++) {
    optionValues[i] = max(
        (strike - spot_ * pow(upFactors[Nnodes - 1], 2 * i - Nnodes)) * type,
        0.0);
  }
  for (int j = Nnodes - 1; j >= 0; j--) {
    double discountFactor = exp(-r_[0] * dt);
    for (int i = 0; i <= j; i++) {
      if (style == European)
        optionValues[i] =
            discountFactor * ((1.0 - upProbabilities[j]) * optionValues[i] +
                              upProbabilities[j] * optionValues[i + 1]);
      else {
        double exerciseValue =
            max((strike - spot_ * pow(upFactors[j], double(2 * i - j))) * type,
                0.0);
        optionValues[i] =
            max(exerciseValue,
                discountFactor * ((1.0 - upProbabilities[j]) * optionValues[i] +
                                  upProbabilities[j] * optionValues[i + 1]));
      }
    }
  }
  return optionValues[0];
};

double CTree::calculateBinomial(double strike, double Maturity, int Nnodes,
                                exStyle style, OptionType pay, tType /*tree*/,
                                Vdoub /*schedule*/) {
  double dt = Maturity / double(Nnodes);
  if (r_.empty())
    r_.push_back(0.0);
  if (q_.empty())
    q_.push_back(0.0);
  Mdoub myTree;
  Vdoub optionValues(Nnodes + 1, 0);
  Vdoub upFactors(Nnodes), upProbabilities(Nnodes);

  for (int i = 0; i < Nnodes; i++) {
    upFactors[i] = exp(interp_vol(T_, IV_, i * dt, (i + 1) * dt) * sqrt(dt));
    double downFactor = 1.0 / upFactors[i];
    upProbabilities[i] =
        (exp((r_[0] - q_[0]) * dt) - downFactor) / (upFactors[i] - downFactor);
  }

  double type = (pay == OptionType::Call) ? -1.0 : 1.0;
  for (int i = 0; i <= Nnodes; i++) {
    optionValues[i] = max(
        (strike - spot_ * pow(upFactors[Nnodes - 1], 2 * i - Nnodes)) * type,
        0.0);
  }
  for (int j = Nnodes - 1; j >= 0; j--) {
    double discountFactor = exp(-r_[0] * dt);
    for (int i = 0; i <= j; i++) {
      if (style == European)
        optionValues[i] =
            discountFactor * ((1.0 - upProbabilities[j]) * optionValues[i] +
                              upProbabilities[j] * optionValues[i + 1]);
      else {
        double exerciseValue =
            max((strike - spot_ * pow(upFactors[j], double(2 * i - j))) * type,
                0.0);
        optionValues[i] =
            max(exerciseValue,
                discountFactor * ((1.0 - upProbabilities[j]) * optionValues[i] +
                                  upProbabilities[j] * optionValues[i + 1]));
      }
    }
  }
  return optionValues[0];
};

double CTree::interp_vol(vector<double> &T, vector<double> &F, double t1,
                         double t2) const {
  double sig0 = interpolate(T, F, t1);
  double sigT = interpolate(T, F, t2);
  double dt = t2 - t1;
  if (dt <= 0.0)
    VEL_RAISE("Incorrect dt");
  if (sig0 == sigT)
    return sig0;
  else
    return sqrt((t2 * sigT * sigT - t1 * sig0 * sig0) / dt);
};

double CTree::interpolate(vector<double> &T, vector<double> &F,
                          double t) const {
  int i, N = T.size();
  for (i = 0; i < N; i++)
    if (t < T[i])
      break;
  if (i == 0)
    return F[i] / T[i] * t;
  if ((i > 0) && (i < N))
    return F[i - 1] + (F[i] - F[i - 1]) / (T[i] - T[i - 1]) * (t - T[i - 1]);
  else
    return F[N - 2] +
           (t - T[N - 2]) / (T[N - 1] - T[N - 2]) * (F[N - 1] - F[N - 2]);
};

double CTree::calculateTrinomial(double strike, double Maturity, int Nnodes,
                                 exStyle style, OptionType pay,
                                 tType /*tree*/) {
  double dt = Maturity / double(Nnodes);
  if (r_.empty())
    r_.push_back(0.0);
  if (q_.empty())
    q_.push_back(0.0);
  int valuesSize = 2 * Nnodes + 1;
  Vdoub optionValues(valuesSize, 0);
  Vdoub upFactors(Nnodes), upProbabilities(Nnodes), downProbabilities(Nnodes);
  double midProbability;
  double continuationValue;
  double sqrtHalfDt = sqrt(0.5 * dt);
  for (int i = 0; i < Nnodes; i++) {
    double sig = interp_vol(T_, IV_, i * dt, (i + 1) * dt);
    double sigSqrtHalfDt = sig * sqrtHalfDt;
    double halfRateDt = (r_[0] - q_[0]) * dt * 0.5;
    double denominator = exp(sigSqrtHalfDt) - exp(-sigSqrtHalfDt);
    upFactors[i] = exp(sig * sqrt(2.0 * dt));
    upProbabilities[i] = (exp(halfRateDt) - exp(-sigSqrtHalfDt)) / denominator;
    upProbabilities[i] *= upProbabilities[i];
    downProbabilities[i] = (exp(sigSqrtHalfDt) - exp(halfRateDt)) / denominator;
    downProbabilities[i] *= downProbabilities[i];
  }

  double type = (pay == OptionType::Call) ? -1.0 : 1.0;
  for (int i = 0; i < valuesSize; i++) {
    double payoff =
        (strike - spot_ * pow(upFactors[Nnodes - 1], i - Nnodes)) * type;
    optionValues[i] = max(payoff, 0.0);
  }

  for (int j = Nnodes - 1; j >= 0; j--) {
    double discountFactor = exp(-r_[0] * dt);
    for (int i = 0; i <= 2 * j; i++) {
      midProbability = (1.0 - downProbabilities[j] - upProbabilities[j]);
      continuationValue =
          discountFactor * (downProbabilities[j] * optionValues[i] +
                            midProbability * optionValues[i + 1] +
                            upProbabilities[j] * optionValues[i + 2]);
      if (style == American) {
        double exerciseValue = max(
            (strike - spot_ * pow(upFactors[j], double(i - j))) * type, 0.0);
        optionValues[i] = max(exerciseValue, continuationValue);
      } else
        optionValues[i] = continuationValue;
    }
  }
  return optionValues[0];
};

} // namespace velesquant