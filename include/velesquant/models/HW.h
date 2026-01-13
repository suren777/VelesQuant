//		HW.h

#ifndef HW_H
#define HW_H

#include <vector>
#include <velesquant/models/interpolation.h>
#include <velesquant/models/utility.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class HW
{
public:
	HW(double kappa, std::vector<double> timeSigmas, std::vector<double> sigmas, std::vector<double> timeDFs, std::vector<double> DFs)
		:  kappa_(kappa), timeSigmas_(timeSigmas), sigmas_(sigmas), timeDFs_(timeDFs), DFs_(DFs) {sigmas0_=sigmas_;kappa0_=kappa_;};

	~HW() {};

	double optionBond(double Expiry, double Maturity, double Strike, std::string callORput);
	double swaption(double Expiry, double Tenor, double Strike, double PayFrequency=0.5);
	double BlackStrike(double T0, double TN, double impVol);
	double BlackStrikePlain(double T0, double TN);
	std::vector<double> simulationHW(std::vector<double> times) const;

	double ZC(double Expiry);
	void calibrator(std::vector<defSwap> swapQuotes, const std::string &type);	
	void calibratorBstrp(std::vector<defSwap> swapQuotes, const std::string &type);
	double getKappa() { return kappa_; };
	std::vector<double> getTimeSigmas() { return timeSigmas_; };
	std::vector<double> getSigmas() { return sigmas_; };
	double getSwapRate(double Expiry, double Tenor, double PayFrequency=0.5);
	double swaptionIVblackPub(double Expiry, double Tenor, double swap_price);
	double get_swaptionATM(double Expiry, double Tenor, double VolATM){	return swaptionATM(Expiry, Tenor, VolATM);};
	std::vector<double> getDFs(){return DFs_;};
	std::vector<double> getDFsTimes(){return timeDFs_;};
private:
	double kappa_, kappa0_;
	std::vector<double> timeDFs_, DFs_;  //time_[0] first point
	mutable std::vector<double> timeSigmas_, sigmas_, sigmas0_;
	int iter_;
	mutable std::vector<double> marketSwaption_;
	std::vector<defSwap> quoteSwap_;
	std::vector<defSwap> qSB_;

	double totalVariance(double T0);
	double equalityCP(double vX, double vO, double Expiry, double Tenor, double Strike, double PayFrequency);
	double criticalPoint(double Expiry, double Tenor, double Strike, double PayFrequency);

	double formulaBlack(double dfT0, double dfTN, double Strike, double impVol, double T0, std::string callORput);

	double getDF(double T) { return interpolation("CubicNaturalSpline", timeDFs_, DFs_, T); };
	double swaptionATM(double Expiry, double Tenor, double VolATM)
	{
		return (getDF(Expiry)-getDF(Expiry+Tenor))*(2.0*cdf_normal(0.5*VolATM*sqrt(Expiry))-1.0);
	};

	double whichSigma(double t) const;
	double whichFwdRate(double t) const;

	double pen_fun(double*x, double*lb, double*ub, int n);
	void objFcn(int m, int n, double* x, double* fvec, int* iflag);	
	void objFcnIV(int m, int n, double* x, double* fvec, int* iflag, double* lb, double* ub);	
	//	void objFcnPrice(int m, int n, double* x, double* fvec, int* iflag);
	void objFcnPrice(int m, int n, double* x, double* fvec, int* iflag, double* lb, double* ub);	
	double objFcnBswap(double x, int m);
	double objFcnBIV(double x, int m);
	//not used
	double swap(double Expiry, double Tenor, double Strike, double PayFrequency=0.5);
	double swaptionIVblack(double Expiry, double Tenor, double swaption_price);
	double Bratio(double a, double Mi, double Tj, double Tk);
	void calibrateKappa();
	mutable std::vector<double> iv_ratio,bond_ratio,index;
	double objKappa(double x, const std::vector<double> &ivr, const std::vector<double> &br, const std::vector<int> &ind);
};

}
#endif