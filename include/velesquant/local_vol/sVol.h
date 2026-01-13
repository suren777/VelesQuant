//		sVol.h

#ifndef SVOL_H
#define SVOL_H

#include <vector>
#include <complex>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class sVol
{
public:
	sVol() {};
	sVol(double spot, double var0, double kappa, double theta, double xi, double rho, int seed);
	~sVol() {};
	
	double hestonPrice(double maturity, double forward, double strike, std::string optType) const;
	double hestonPriceCF(double maturity, double forward, double strike, std::string optType) const;
	double hestonPriceNET(double maturity, double forward, double strike, std::string optType) const;
	void calibrator(std::vector<double> maturitys, std::vector<double> forwards, std::vector<double> strikes, 
		            std::vector<double> marketQuotes, std::string quoteType);
	void IVcalibrator(std::vector<double> maturitys, std::vector<double> forwards, std::vector<double> strikes, 
		            std::vector<double> marketQuotes, std::string quoteType);
	void FXcalibrator(std::vector<double> maturitys, std::vector<double> forwards, std::vector<double> strikes, 
		            std::vector<double> marketQuotes, std::string quoteType);
	std::vector<double> simulationHeston(std::vector<double> times, std::vector<double> forwards) const;
	std::vector<double> simulationHestonDO(std::vector<double> times, std::vector<double> forwards, double barrier) const;
	std::vector<double> simulationHestonMax(std::vector<double> times, std::vector<double> forwards) const;
	std::vector<double> simulationHestonCliq(std::vector<double> times, std::vector<double> forwards, double gcap, double gfloor, double lcap, double lfloor, double alpha) const;
	std::vector<double> simulationHestonDNT(std::vector<double> times, std::vector<double> forwards, double UP, double DOWN) const;
	double simulationHestonDNTdt(std::vector<double> times, std::vector<double> forwards, double UP, double DOWN, double dt) const;
	int simulationHestonDNTdtS(std::vector<double> times, std::vector<double> forwards, double maturity, double UP, double DOWN, double dt) const;
	int simulationHestonDNTdtE(std::vector<double> times, std::vector<double> forwards, double maturity, double UP, double DOWN, double dt) const;
	std::vector<double> simulationHestonUpnOut(std::vector<double> times, std::vector<double> forwards, double BARRIER) const;
	std::vector<double> simulationHestonDownnOut(std::vector<double> times, std::vector<double> forwards, double BARRIER) const;
	double getParameterVar0() const {return var0_;};
	double getParameterKappa() const {return kappa_;};
	double getParameterTheta() const {return theta_;};
	double getParameterXi() const {return xi_;};
	double getParameterRho() const {return rho_;};
	void update_seed(int seed);

private:
	double spot_, var0_, kappa_, theta_, xi_, rho_;
	
	std::vector<double> maturitys_, forwards_, strikes_, marketQuotes_;
	double intCFfun(double u, double ki, double X, double maturity)const;
	std::vector<double> weights_;

	void objFcn(int m, int n, double* x, double* fvec, int* iflag);
	void objFcnIV(int m, int n, double* x, double* fvec, int* iflag);
	void objFcnFX(int m, int n, double* x, double* fvec, int* iflag);

	std::complex<double> hestonCF(std::complex<double> k, double maturity) const;
	double hestonIntegrand(std::complex<double> k, double maturity, double forward, double strike, int choice) const;
	double trapz(std::vector<double> x, std::vector<double> y) const;
	
	void setParameterVar0(double var0) 
		{
		var0_=fabs(var0);
		};
	void setParameterKappa(double kappa) {kappa_=fabs(kappa);};
	void setParameterTheta(double theta) {theta_=fabs(theta);};
	void setParameterXi(double xi) {xi_=fabs(xi);};
	void setParameterRho(double rho) {if (abs(rho)<=1.0) rho_=rho; else rho_=sin(rho);};

	double interpolateF(std::vector<double> T, std::vector<double> F, double t) const;

	double hestonIntegrandNET(std::complex<double> k, double maturity, double forward, double strike) const;
};

}
#endif