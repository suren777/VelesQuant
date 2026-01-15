
#include <ql/quantlib.hpp>
#include <velesquant/numerics/interpolation.h>
#include <velesquant/types.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace velesquant::xlw;
using namespace std;

namespace velesquant {

double doInterpolation(string method, CellMatrix &xVec, CellMatrix &yVec,
                       double x) {
  if ((xVec.ColumnsInStructure() != 1) || (yVec.ColumnsInStructure() != 1))
    throw("one columns required!");
  int Ix = xVec.RowsInStructure();
  std::vector<double> vecX(Ix);
  for (int i = 0; i < Ix; ++i)
    vecX[i] = xVec(i, 0).NumericValue();
  int Iy = yVec.RowsInStructure();
  std::vector<double> vecY(Iy);
  for (int i = 0; i < Iy; ++i)
    vecY[i] = yVec(i, 0).NumericValue();
  return interpolation(method, vecX, vecY, x);
}

double interpolation(string method, std::vector<double> vecX,
                     std::vector<double> vecY, double x) {
  if (method == "CubicNaturalSpline") {
    QuantLib::CubicNaturalSpline interp(vecX.begin(), vecX.end(), vecY.begin());
    return interp(x);
  } else if (method == "LinearInterpolation") {
    QuantLib::LinearInterpolation interp(vecX.begin(), vecX.end(),
                                         vecY.begin());
    return interp(x);
  } else {
    QuantLib::LinearInterpolation interp(vecX.begin(), vecX.end(),
                                         vecY.begin());
    return interp(x);
  }
}
} // namespace velesquant