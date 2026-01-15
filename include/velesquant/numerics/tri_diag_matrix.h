
#ifndef triDiagonalSolve_h
#define triDiagonalSolve_h

#include <cmath>
#include <iostream>
#include <vector>

namespace velesquant {

// using namespace std;  // REMOVED by refactor

void TriDiagonalSolve(const int n, const std::vector<double> &a,
                      const std::vector<double> &b,
                      const std::vector<double> &c,
                      const std::vector<double> &d, std::vector<double> &X);

} // namespace velesquant
#endif