#ifndef HESTON_FACADE_H
#define HESTON_FACADE_H

#include <string>
#include <velesquant/types.h>

namespace velesquant {

std::string CreateHestonObj(const std::string &theObjName, double spot,
                            double var0, double kappa, double theta, double xi,
                            double rho, int seed);

double hestonPrice(std::string theName, double maturity, double forward,
                   double strike);

Matrix hestonCalibrator(const std::string &theName, Matrix theMaturitys,
                        Matrix theForwards, Matrix theStrikes, Matrix theQuotes,
                        std::string quoteType);

Matrix hestonSimulation(const std::string &theName, Matrix theTimes,
                        Matrix theForwards, int NP);
} // namespace velesquant

#endif
