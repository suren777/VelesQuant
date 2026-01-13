//!		Heston Time-Dependent Model Class
/*!
	Heston Time Dependent model allows time dependency for mean-reverting parameter for volatility process
	as well as for volatility of volatility.


*/


#ifndef SVOLT_H
#define SVOLT_H

#include <vector>
#include <complex>
#include <velesquant/local_vol/termstructure.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor


//! Model types
/*!
	This enumeration provides a choise of model types to use

*/
enum model_types{volvol, meanrev, full};

class sVolT
	{
	public:

		//! Class constructor
		sVolT(double spot, double var0, double kappa,double rho, std::vector<double>  &time, std::vector<double> &theta, std::vector<double> &xit, int seed); ///
		sVolT(double spot, double var0, double kappa,double rho, std::vector<double>  &time, std::vector<double> &theta,  double xit, int seed); 
		sVolT(double spot, double var0, double kappa,double rho, std::vector<double>  &time, double theta, std::vector<double> &xit, int seed); 
		~sVolT() {};

		double hestonPriceTd(double maturity, double forward, double strike, std::string optType) const;
		double hestonPriceTdCF(double maturity, double forward, double strike, std::string optType) const;

		double	getParameterVar0() const {return var0_;};
		double	getParameterKappa() const {return kappa_;};
		double	getParameterRho() const {return rho_;};
		std::vector<double>	getXi()  {return xit_; };
		std::vector<double>	getTheta() {return thetat_;};
		std::vector<double>	getT() {return t_;};
		double	getParameterXi(int i) const {if (i<nT_) return xit_[i]; else throw("Memory Violation");};
		double	getParameterTheta(int i) const {if (i<nT_) return thetat_[i]; else throw("Memory Violation");};
		int		getParameternT() const {return nT_;}
		int		getType() const {return model_;}

		void setParameterVar0(double var0) {var0_=fabs(var0);};
		void setParameterKappa(double kappa) {kappa_=fabs(kappa);};
		void setParameterRho(double rho) 
		{
			if (abs(rho)<=1) rho_=rho;
			else rho_=sin(rho);
		};
		void setParameterXi(double xit,int i)
			{
			if (i < nT_) xit_[i]=fabs(xit);
			}
		void setParameterTheta(double theta,int i)
			{
			if (i < nT_) thetat_[i]=fabs(theta);
			}

		
		void calibratorLM(std::vector<double> maturities, std::vector<double> forwards, std::vector<double> strikes, 
			std::vector<double> marketQuotes, std::string quoteType);

		//	void calibratorNM(vector<double> maturities, vector<double> forwards, vector<double> strikes, 
		//		vector<double> marketQuotes, string quoteType);


		std::vector<double> simulationHestonTd(std::vector<double> times, std::vector<double> forwards) const;
		std::vector<double> simulationHestonTdMax(std::vector<double> times, std::vector<double> forwards) const;
		std::vector<double> simulationHestonTdCliq(std::vector<double> times, std::vector<double> forwards, double gcap, double gfloor, double lcap, double lfloor, double alpha) const;
		std::vector<double> simulationHestonTdCliq(std::vector<double> times, std::vector<double> forwards, 
												double gcap, double gfloor, 
												double lcap, double lfloor, 
												double alpha,
												int Day,
												Termstructure *T
												) const;
		void update_seed(int seed);
	
	private:
		double spot_, var0_, kappa_, rho_;

		std::vector<double> t_;		// times
		std::vector<double> xit_;	// etas
		std::vector<double> thetat_; // thetas

		int nT_;

		void createT(std::vector<double> &time)
			{
			nT_ = time.size();
			for (int i=0 ; i < nT_ ; i++)
				t_.push_back(time[i]);
			}



		void createXitT(std::vector<double> &xin){
			for (int i=0 ; i < nT_ ; i++)				
				xit_.push_back(abs(xin[i]));
			};

		void createXitT(double xin){
			for (int i=0 ; i < nT_ ; i++)				
				xit_.push_back(abs(xin));
			};

		void createthetat(std::vector<double> &thetat)
			{
			for (int i=0 ; i < nT_ ; i++)
				thetat_.push_back(abs(thetat[i]));

			};
		void createthetat(double thetat)
			{
			for (int i=0 ; i < nT_ ; i++)
				thetat_.push_back(abs(thetat));

			};

		std::vector<double> maturities_, forwards_, strikes_, marketQuotes_;

		double trapz(std::vector<double> x, std::vector<double> y) const;
		double hestonIntegrandTd(std::complex<double> k, double maturity, 
			double forward, double strike, int choice) const;
		double intCFfun(double u, double ki, double X, double maturity) const;
		std::complex<double> hestonCFTd(std::complex<double> k, double maturity) const;
		void objFcnTd(int m, int n, double* x, double* fvec, int* iflag);
		//		double objFcnNM(vector<double>);
		double interpolateF(std::vector<double> T, std::vector<double> F, double t) const;
		model_types model_;

	};

}
#endif