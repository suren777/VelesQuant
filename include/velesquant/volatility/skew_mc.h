//      skewMC.h

#ifndef SKEWMC_H
#define SKEWMC_H

#include <string>
#include <vector>
#include <velesquant/volatility/sabr.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class skewMC {
public:
  skewMC(std::vector<double> maturities, std::vector<double> forwards,
         std::vector<double> betas, std::vector<double> alphas,
         std::vector<double> nus, std::vector<double> rhos);
  skewMC(std::vector<double> maturities, std::vector<double> forwards,
         std::vector<double> betas, std::vector<double> alphas,
         std::vector<double> nus, std::vector<double> rhos,
         std::vector<double> shifts);
  skewMC(std::vector<Sabr> sabrModels) : sabrModels_(sabrModels) {};
  skewMC() {};
  ~skewMC() {};

  std::vector<double> simulation(std::vector<double> times, double spot,
                                 double kappa);
  double simulatedSpot(double time, double spot, double corrRN);

private:
  std::vector<Sabr> sabrModels_;

  std::vector<std::vector<double>> autoCorrMatrix(std::vector<double> times,
                                                  double kappa) const;
  std::vector<double> autoCorrRandoms(std::vector<double> times,
                                      double kappa) const;
  std::vector<std::vector<double>>
  cholesky(std::vector<std::vector<double>> corrMatrix) const;
};

} // namespace velesquant
#endif