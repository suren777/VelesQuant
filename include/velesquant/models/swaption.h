//		swaption.h

#ifndef SWAPTION_H
#define SWAPTION_H

#include <string>
#include <velesquant/local_vol/sabr.h>
#include <velesquant/xlw/CellMatrix.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class swaption
{
public:
	swaption(double expiry, double tenor, double forward, double annuity, 
		double beta=0.85, double alpha=0.5, double nu=0.25, double rho=-0.75);

	swaption(double expiry, double tenor, double forward, double annuity,
		double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType="premium");
	
	swaption(double expiry, double tenor, double forward, double annuity, 
		double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType, 
		CellMatrix& initialParams);

	~swaption() { delete swaptionSABR_; };
	
	double swaptionFairValue(double strike, std::string callORput="call") const;
	double swapFairValue(double strike) const;
		
	double simulation(double corrRN);

	double getImpliedVol(double strike) const;

	double getParameterAlpha() const;
	double getParameterNu() const;
	double getParameterRho() const;

private:
	double annuity_;
	sabr* swaptionSABR_;
};

}
#endif