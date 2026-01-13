//		ShortRate1FPDE.h

#ifndef ShortRate1FPDE_H
#define ShortRate1FPDE_H

#include <velesquant/models/utility.h>
#include <vector>
#include <list>

namespace velesquant {

// using namespace std; // REMOVED by refactor
#pragma warning (disable:4018)

class ShortRate1FPDE
{
public:
	ShortRate1FPDE(double R0, double kappa, double alpha, double beta, double gamma, 
		std::vector<double> timeSigmas, std::vector<double> sigmas, 
		std::vector<double> timeThetas, std::vector<double> thetas)
		:  R0_(R0), kappa_(kappa), alpha_(alpha), beta_(beta), gamma_(gamma), 
		timeSigmas_(timeSigmas), sigmas_(sigmas), 
		timeThetas_(timeThetas), thetas_(thetas)	{in_calibration_=false; buildGrid();};

	ShortRate1FPDE(std::vector<double> timeDFs, std::vector<double> DFs, std::vector<defSwap> quoteSwap)
	{
		R0_ = -log(DFs[0])/timeDFs[0];
		kappa_ = 0.02;
		alpha_ = 1.0;
		beta_ = 0.0001;
		gamma_ = 1.0;
		buildGrid();
		in_calibration_=false;
		calibrator(timeDFs, DFs, quoteSwap);
	};

	~ShortRate1FPDE() {};

	void calibrator(std::vector<double> timeDFs, std::vector<double> DFs, std::vector<defSwap> swapQuotes);

	double pricingZB(double Maturity);
	void discountBack(double Expiry, double Maturity, std::vector<double> &f);
	double pricingZBO(double Expiry, double Maturity, double Strike, const std::string &type);
	double pricingCouponBond(double Expiry, double Tenor, double Coupon, double PayFrequency);
	double pricingCBO(double Expiry, double Tenor, double Coupon, double Strike, double PayFrequency, const std::string &type);
	double pricingSwap(double Expiry, double Tenor, double Strike, double PayFrequency);
	double pricingCallableSwap(double Expiry, double Tenor, std::vector<double> Exercises, double Coupon, double Strike, double PayFrequency, const std::string &type);
	void pricingCouponBondt(double Expiry, double Tenor, double Coupon, double PayFrequency, std::vector<double> &f);
	double pricingSwaption(double Expiry, double Tenor, double Strike, double PayFrequency=0.5);
	std::vector<double> calculateDFs(std::vector<double>& timePoints);
	double swaptionATM(double Expiry, double Tenor, double VolATM)
	{
		return (getDFinterp(Expiry)-getDFinterp(Expiry+Tenor))*(2.0*cdf_normal(0.5*VolATM*sqrt(Expiry))-1.0);
	};
	double pricingBermudan(double Expiry, double Tenor, std::vector<double> Exercises, double Strike, double PayFrequency);
	double getSwapRate(double Expiry, double Tenor, double PayFrequency);
	double getR0() { return R0_; };
	double getKappa() { return kappa_; };
	double getAlpha() { return alpha_; };
	double getBeta() { return beta_; };
	double getGamma() { return gamma_; };
	std::vector<double> getTimeSigmas() { return timeSigmas_; };
	std::vector<double> getSigmas() { return sigmas_; };
	std::vector<double> getTimeThetas() { return timeThetas_; };
	std::vector<double> getThetas() { return thetas_; };
	std::vector<double> getGrid() { return gridR_; };
	void change_grid_factor(double Rmax, double factor){buildGrid(Rmax, factor);};
	std::vector<double> SwaptionDiagnostic(double Expiry, double Tenor, double Strike, double PayFrequency);
	std::vector<double> RiskNeutralDensity(double T1, double T2);
	double getImpVolATM(double Expiry, double Tenor, double PayFrequency);
private:
	double R0_, kappa_, alpha_, beta_, gamma_;
	double vleft_, vright_;
	bool in_calibration_;
	int iR0_;
	mutable std::vector<double> timeSigmas_, sigmas_;
	mutable std::vector<double> timeThetas_, thetas_, DFinterp_;
	mutable std::vector<double> a_,b_,c_,e_,f_,g_;
	std::vector<double> DFs_, timeDFs_;
	std::vector<defSwap> quoteSwap_;

	mutable std::vector<double> gridT_;
	mutable std::vector<double> gridR_;
	void buildGrid(double Rmax = 1 ,double factor = 0.7);
	void oneStepBackward(const int t, std::vector<double> &inV);
	void oneStepForward(const int T, std::vector<double> &inV);

	void termStructureCalibrator();
	void objFcnCalibrator(int m, int n, double* x, double* fvec, int* iflag);

	double trapezoidal(std::vector<double>& inV);
	double whichSigma(double t) const;
	double whichTheta(double t) const;
	double getDFinterp(double t);
	std::vector<double> grMesh(int,double Rmax = 0.5 ,double factor = 0.7);
	double pen_fun(double*x, int n);
	void calibrateKappa();
	double objKappa(double a, const std::vector<double> &ivr, const std::vector<double> &br, const std::vector<int> &ind);
	double Bratio(double a, double Mi, double Tj, double Tk);
	double F1(int j, double p, double theta);
	double F2(int j, double sigma2);
	double F3(int j, double p);


};

}
#endif