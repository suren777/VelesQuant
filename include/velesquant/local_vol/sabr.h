//		sabr.h

#ifndef SABR_H
#define SABR_H

#include <string>
#include <vector>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class sabr
	{
	public:
		sabr() {};
		sabr(double maturity, double forward, double beta=0.85, 
			double alpha=0.5, double nu=0.25, double rho=-0.75);
		sabr(double maturity, double forward, double beta, 
			double alpha, double nu, double rho, double shift);
		~sabr() {};

		double impliedVol(double strike) const;
		double normalVol(double K) const; 
		double premiumBachelier(double strike, std::string callORput="call") const;
		double premiumBlackScholes(double strike, std::string callORput="call") const;
		double localVol(double spot) const;
		double localVolCall(double spot) const; //based on the call premium
		double localVolzabr(double spot) const;	
		double localVolIV(double spot) const;   //based on the implied vol

		void qutlTable();
		double simulation(double corrRN);
		std::vector<double> simulations(std::vector<double> correlatedRNs);
		std::vector<double> getTableQUTL1() const { return spots_;};
		std::vector<double> getTableQUTL2() const { return qutls_;};
		std::vector<double> getTableQUTL3() const { return cdfs_;};
		double getVol(double strike) const;
		double getPremium(double strike, std::string callORput="call") const;
		void calibratorNormalVol(std::vector<double> strikes, std::vector<double> marketQuotes);
		void calibrator(std::vector<double> strikes, std::vector<double> marketQuotes, 
			std::string quoteType="premium");
		//void calibratorWithInitial(std::vector<double> strikes, std::vector<double> marketQuotes, 
		//	std::string quoteType, std::vector<double> initialParams);
		void calibratorWithInitial(std::vector<double> strikes, std::vector<double> marketQuotes, 
			std::string quoteType);
		void calibratorWithInitialATM(std::vector<double> strikes, std::vector<double> marketQuotes, std::string quoteType);

		void setParameterAlpha(double alpha) {alpha_=alpha;};
		void setParameterNu(double nu) {nu_=nu;};
		void setParameterRho(double rho) {if (rho*rho<=1) rho_=rho;else rho_=sin(rho);};
		double getParameterAlpha() const { return alpha_;};
		double getParameterNu() const { return nu_;};
		double getParameterRho() const { return rho_;};
		void setMaturity(double maturity) {maturity_=maturity;};
		void setForward(double forward) {forward_=forward;};
		void setBeta(double beta) {beta_=beta;};
		void setATMvol(double atmVol) {atmVol_ = atmVol;};

		double getMaturity() const { return maturity_;};
		double getForward() const { return forward_;};
		double getBeta() const { return beta_;};
		double getShift() const { return shift_;};

	private:
		bool calibrated_;
		double maturity_, forward_, beta_, shift_;
		double alpha_, nu_, rho_;
		double atmVol_;
		std::string type_;
		std::vector<double> strikes_; 
		std::vector<double> marketQuotes_;

		
		double ATMvolPoly(double alpha);
		double ATMvolRoots(double lBound, double uBound, double tol);

		bool notQUTL_;
		std::vector<double> spots_;
		std::vector<double> qutls_; 
		std::vector<double> cdfs_; 
		void objFcn(int m, int n, double* x, double* fvec, int* iflag, std::string quoteType ="premium");
		void objFcnATM(int m, int n, double* x, double* fvec, int* iflag, std::string quoteType);

		void parameters()
			{
			if( !((alpha_>=0.0) && (nu_>=0.0) && (rho_*rho_<=1.0)) )
				throw("SABR parameters are not valid");
			calibrated_ = true;
			};

		int amin(std::vector<double> q, int position);
		int amax(std::vector<double> q, int position);
		void asort(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c);
		//Mesh initialisation
		std::vector<double> gsMesh(int);
		std::vector<double> grMesh(int);
		//Normal volatilities:
		double Nvolb0(double K) const;		//	beta == 1
		double Nvolb(double K) const;		//	0 < beta <1
		double Nvolb1(double K) const;		//	beta == 1
		//Obj function
		void objFcnNormalVolCalibration(int m, int n, double* x, double* fvec, int* iflag);
	};

}
#endif