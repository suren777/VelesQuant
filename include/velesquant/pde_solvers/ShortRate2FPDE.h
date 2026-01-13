//		ShortRate2FPDE.h

#ifndef ShortRate2FPDE_H
#define ShortRate2FPDE_H


#include <velesquant/models/utility.h>
#include <vector>
#include <list>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class ShortRate2FPDE
{
public:
	ShortRate2FPDE(double kappa1, double kappa2, double lambda, 
		std::vector<double> timeSigma1s, std::vector<double> sigma1s, std::vector<double> timeSigma2s, std::vector<double> sigma2s, 
		std::vector<double> timeAlphas, std::vector<double> alphas)
		: kappa1_(kappa1), kappa2_(kappa2), lambda_(lambda), 
		  timeSigma1s_(timeSigma1s), sigma1s_(sigma1s), timeSigma2s_(timeSigma2s), sigma2s_(sigma2s), 
		  timeAlphas_(timeAlphas), alphas_(alphas) {};

	~ShortRate2FPDE() {};
	
	void calibrator(std::vector<double> timeDFs, std::vector<double> DFs, std::vector<defSwap> swapQuotes);

	double pricingZB(double Maturity);
	double pricingSwaption(double Expiry, double Tenor, double Strike, double PayFrequency=0.5);
	std::vector<double> calculateDFs(std::vector<double>& timePoints);
	
	double getKappa1() { return kappa1_; };
	double getKappa2() { return kappa2_; };
	double getLambda() { return lambda_; };
	std::vector<double> getTimeSigma1s() { return timeSigma1s_; };
	std::vector<double> getSigma1s() { return sigma1s_; };
	std::vector<double> getTimeSigma2s() { return timeSigma2s_; };
	std::vector<double> getSigma2s() { return sigma2s_; };
	std::vector<double> getTimeAlphas() { return timeAlphas_; };
	std::vector<double> getAlphas() { return alphas_; };

private:
	double kappa1_, kappa2_, lambda_;
	std::vector<double> timeSigma1s_, sigma1s_;
	std::vector<double> timeSigma2s_, sigma2s_;
	mutable std::vector<double> timeAlphas_, alphas_;

	std::vector<double> DFs_;
	std::vector<defSwap> quoteSwap_;
	
	void termStructureCalibrator();
	void objFcnCalibrator(int m, int n, double* x, double* fvec, int* iflag);

	mutable std::vector<double> gridT_;
	mutable std::vector<double> gridX_, gridY_;
	void buildGrid(double Time, int Nt=5000);
	
	void oneStepBackwardADIDouglas(int t, const std::vector<std::vector<double> >&inM, std::vector<std::vector<double> >&outM);	
	void oneStepBackwardExplicit(int t, const std::vector<std::vector<double> >&inM, std::vector<std::vector<double> >&outM);
	void oneStepBackwardDouglasX(int t, const std::vector<std::vector<double> >&inM, const std::vector<std::vector<double> >&midM, std::vector<std::vector<double> >&outM);
	void oneStepBackwardDouglasY(int t, const std::vector<std::vector<double> >&inM, const std::vector<std::vector<double> >&midM, std::vector<std::vector<double> >&outM);
		
	void oneStepForwardADIDouglas(int T, const std::vector<std::vector<double> >&inM, std::vector<std::vector<double> >&outM);	
	void oneStepForwardExplicit(int T, const std::vector<std::vector<double> >&inM, std::vector<std::vector<double> >&outM);
	void oneStepForwardDouglasX(int T, const std::vector<std::vector<double> >&inM, const std::vector<std::vector<double> >&midM, std::vector<std::vector<double> >&outM);
	void oneStepForwardDouglasY(int T, const std::vector<std::vector<double> >&inM, const std::vector<std::vector<double> >&midM, std::vector<std::vector<double> >&outM);
		
	double trapezoidal2D(std::vector<std::vector<double> >& inV);

	double whichValue(double t, const std::vector<double> times, const std::vector<double> values);
};

}
#endif