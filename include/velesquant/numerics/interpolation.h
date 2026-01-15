//      interpolation.h

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <string>
#include <vector>

namespace velesquant {

// using namespace std; // REMOVED by refactor

double interpolation(std::string method, std::vector<double> vecX,
                     std::vector<double> vecY, double x);

} // namespace velesquant
#endif