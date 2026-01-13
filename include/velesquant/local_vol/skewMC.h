//      skewMC.h

#ifndef SKEWMC_H
#define SKEWMC_H

#include <velesquant/local_vol/sabr.h>
#include <string>
#include <vector>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class skewMC
{
public:
	skewMC(std::vector<double> maturities, std::vector<double> forwards, std::vector<double> betas, 
		   std::vector<double> alphas, std::vector<double> nus, std::vector<double> rhos);
	skewMC(std::vector<double> maturities, std::vector<double> forwards, std::vector<double> betas, 
		   std::vector<double> alphas, std::vector<double> nus, std::vector<double> rhos, std::vector<double> shifts);
	skewMC(std::vector<sabr> sabrModels) : sabrModels_(sabrModels) {};
	skewMC() {};
	~skewMC() {};
	
	std::vector<double> simulation(std::vector<double> times, double spot, double kappa);
	double simulatedSpot(double time, double spot, double corrRN);

private:
	std::vector<sabr> sabrModels_;
	
	std::vector<std::vector<double> > autoCorrMatrix(std::vector<double> times, double kappa) const;
	std::vector<double> autoCorrRandoms(std::vector<double> times, double kappa) const;
	std::vector<std::vector<double> > cholesky(std::vector<std::vector<double> > corrMatrix) const;
};

}
#endif