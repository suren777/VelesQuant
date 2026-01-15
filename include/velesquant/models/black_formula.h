// BlackFormula.h

#ifndef BLACKFORMULA_H
#define BLACKFORMULA_H

#include <cmath>
#include <cstdlib>
#include <velesquant/models/utility.h>

namespace velesquant {

inline double BlackFormulaCall(double Forward, double Strike, double Vol,
                               double Expiry) {
  double standardDeviation = Vol * sqrt(Expiry);
  double d1 =
      log(Forward / Strike) / standardDeviation + 0.5 * standardDeviation;
  double d2 = d1 - standardDeviation;
  return Forward * cdf_normal(d1) - Strike * cdf_normal(d2);
};

inline double BlackFormulaCallVega(double Forward, double Strike, double Vol,
                                   double Expiry) {
  double standardDeviation = Vol * sqrt(Expiry);
  double d1 =
      log(Forward / Strike) / standardDeviation + 0.5 * standardDeviation;
  return Forward * sqrt(Expiry) * pdf_normal(d1);
};

class BlackCall {
public:
  BlackCall(double T_, double Forward_, double Strike_)
      : T(T_), Forward(Forward_), Strike(Strike_) {};

  double operator()(double Vol) const {
    return BlackFormulaCall(Forward, Strike, Vol, T);
  };
  double Value(double Vol) const { return operator()(Vol); };
  double Vega(double Vol) const {
    return BlackFormulaCallVega(Forward, Strike, Vol, T);
  };

private:
  double T, Forward, Strike;
};

} // namespace velesquant
#endif