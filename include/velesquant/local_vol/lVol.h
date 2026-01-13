//		lVol.h

#ifndef LVOL_H
#define LVOL_H

#include <vector>
#include <velesquant/local_vol/sabr.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class lVol
{
public:
	lVol() {};
	lVol(std::vector<sabr> sabrModels) : sabrModels_(sabrModels) {};
	lVol(std::vector<double> maturities, std::vector<double> forwards, std::vector<double> betas, 
		 std::vector<double> alphas, std::vector<double> nus, std::vector<double> rhos, double spot);
	~lVol() {};
	
	std::vector<double> simulation(std::vector<double> times) const;
	std::vector<double> simulationTEST(std::vector<double> times) const;
	std::vector<double> simulation(std::vector<double> times, std::vector<double> rands) const;
	
	void buildGrid(const std::vector<double> times);
	void buildGrid(double lastTime, int N=100);
		
	void oneStepBackward(const int t, const std::vector<double> &inV, std::vector<double> &outV);
	void oneStepBackwardDirichlet(const int t, const std::vector<double> &inV, std::vector<double> &outV);

	double callPDE(double maturity, double strike, const int N=100);
	double putPDE(double maturity, double strike, const int N=100);
	double dntPDE(double maturity, double upperBarrier, double lowerBarrier, const int Nt=100);
	
	void backwardOneStep(const int t, std::vector<double> &y, const std::vector<double> &z);
	void forwardOneStep(const int t, const std::vector<double> &y, std::vector<double> &z);
    std::vector<double> density(double maturity, const int Nt);

	std::vector<std::vector<double> > exportLV(const std::vector<double> times);

private:
	double spot_;
	std::vector<sabr> sabrModels_;
	
	mutable std::vector<double> gridT_;
	mutable std::vector<std::vector<double> > gridFwd_;
	mutable std::vector<std::vector<double> > gridLV_;
		
	double premiumCALL(double strike, double forward, double maturity, double vol) const;
    double interpolatedCall(const sabr& nextModel, double strike, double forward, double maturity) const;
    double interpolatedCall(const sabr& preModel, const sabr& nextModel, double strike, double forward, double maturity) const;
	
	double getForward(double time) const;
	double getLocVol(double preTime, double preSpot, double preFwd) const;

	void getLocVol(double time, double preSpot, double preFwd, double &locvol, double &forward) const;
	void getLocVolOLD(double time, double preSpot, double preFwd, double &locvol, double &forward) const;
};

}
#endif