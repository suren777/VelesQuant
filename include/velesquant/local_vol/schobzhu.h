// schobzhu

//		schobzhu.h

#ifndef SCHOBZHU_H
#define SCHOBZHU_H

#include <vector>
#include <complex>

namespace velesquant {

typedef std::complex<double> Cdoub;
typedef std::vector<double> Vdoub;

class SchobZhu	
	{
	public:

		SchobZhu(double spot, double var0, double kappa, double theta,  double xi,  double rho);								// parameters constructor
		~SchobZhu(){};



		// Schobel&Zhu integrand_1 and integrand_2


		double SchobelPrice(double maturity, double forward, double strike) const;
		Vdoub simulationSchobZhu(Vdoub times, Vdoub forwards) const;


		double getParameterVar0() const{return sigma0_;};
		double getParameterKappa() const{return kappa_;};
		double getParameterTheta() const{return theta_; };
		double getParameterXi() const {return xi_;};
		double getParameterRho() const{return rho_;};


		void setParameterVar0(double var0) {sigma0_=var0;};
		void setParameterKappa(double kappa) {kappa_=fabs(kappa);};
		void setParameterTheta(double theta) {theta_=fabs(theta);};
		void setParameterXi(double xi) {xi_=fabs(xi);};
		void setParameterRho(double rho) {rho_=sin(rho);};


		void calibrator(Vdoub maturitys, Vdoub forwards, Vdoub strikes, Vdoub marketQuotes);

	private:
		double s0_;					// Spot price
		double sigma0_;					// Initial volatility
		double kappa_;					// Mean reversion rate for volatility
		double theta_;					// Long run volatility
		double xi_;						// Volatility of Volatility
		double rho_;						// Price-volatility correlation
		Vdoub maturitys_, forwards_, strikes_, marketQuotes_;
		double SchobelIntegrand(double k, double maturity, double forward, double strike, int type) const;
		void objFcn(int m, int n, double* x, double* fvec, int* iflag);
	};

}
#endif