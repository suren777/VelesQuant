#ifndef SABR_FACADE_H
#define SABR_FACADE_H

#include <string>
#include <velesquant/types.h>

namespace velesquant {

std::string CreateSABRObj(const std::string &theObjName, double maturity,
                          double forward, double beta, double alpha, double nu,
                          double rho);

std::string dvlCreateSABRObjShifted(const std::string &theObjName,
                                    double maturity, double forward,
                                    double beta, double alpha, double nu,
                                    double rho, double shift);

double sabrVolatility(double strike, double forward, double maturity,
                      double alpha, double beta, double nu, double rho);

Matrix sabrCalibrator(const std::string &theName, Matrix theStrikes,
                      Matrix theQuotes, std::string quoteType);

Matrix sabrCalibratorWithInitial(const std::string &theName, Matrix theStrikes,
                                 Matrix theQuotes, std::string quoteType,
                                 std::string atmFlag);

} // namespace velesquant

#endif
