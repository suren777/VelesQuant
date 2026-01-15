
namespace velesquant {
// solver.h

#ifndef SOLVER_H
#define SOLVER_H

template <typename T>
double myBisection(double Target, double Low, double High, double Tolerance,
                   const T &TheFunction) {
  double x = 0.5 * (Low + High);
  double y = TheFunction(x);
  int Step = 0;
  do {
    if (y < Target)
      Low = x;
    if (y > Target)
      High = x;
    x = 0.5 * (Low + High);
    y = TheFunction(x);
    if ((High - Low) < Tolerance / 10)
      Step++;
  } while ((fabs(y - Target) > Tolerance) && (Step < 50));
  return x;
}

template <typename T, double (T::*Price)(double) const,
          double (T::*Derivative)(double) const>
double myNewtonRaphson(double Target, double Start, double Tolerance,
                       const T &TheObject) {
  double x = Start;
  double y = (TheObject.*Price)(x);
  while (fabs(y - Target) > Tolerance) {
    double d = (TheObject.*Derivative)(x);
    x += (Target - y) / d;
    y = (TheObject.*Price)(x);
  }
  return x;
}
}
#endif