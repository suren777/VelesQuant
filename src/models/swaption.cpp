#include <velesquant/models/swaption.h>

#include <vector>
#include <string>
#include <velesquant/local_vol/sabr.h>
#include <velesquant/xlw/CellMatrix.h>
#include <algorithm>


using namespace std;
using namespace velesquant::xlw;

namespace velesquant {


swaption::swaption(double expiry, double tenor, double forward, double annuity,
		double beta, double alpha, double nu, double rho)
	: annuity_(annuity)
{
	swaptionSABR_ = new sabr(expiry,forward,beta,alpha,nu,rho);
};

swaption::swaption(double expiry, double tenor, double forward, double annuity, 
		double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType)
	: annuity_(annuity)
{
	swaptionSABR_ = new sabr(expiry,forward,beta);
	int m = strikes.RowsInStructure();
	vector<double> theStrikes(m), theQuotes(m);
	for(int i=0; i<m; ++i) 
	{
		theStrikes[i] = strikes(i,0).NumericValue();
		theQuotes[i] = marketQuotes(i,0).NumericValue();
	}
	swaptionSABR_->calibrator(theStrikes,theQuotes,quoteType);
};
	
swaption::swaption(double expiry, double tenor, double forward, double annuity, 
		double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType, 
		CellMatrix& initialParams)
	: annuity_(annuity)
{
	swaptionSABR_ = new sabr(expiry,forward,beta);
	int m = strikes.RowsInStructure();
	vector<double> theStrikes(m), theQuotes(m);
	for(int i=0; i<m; ++i) 
	{
		theStrikes[i] = strikes(i,0).NumericValue();
		theQuotes[i] = marketQuotes(i,0).NumericValue();
	}
	vector<double> theParams(3);
	theParams[0] = initialParams(0,0).NumericValue();
	theParams[1] = initialParams(1,0).NumericValue();
	theParams[2] = initialParams(2,0).NumericValue();
	//swaptionSABR_->calibratorWithInitial(theStrikes,theQuotes,quoteType,theParams);
	swaptionSABR_->calibratorWithInitial(theStrikes,theQuotes,quoteType);
};

double swaption::swaptionFairValue(double strike, std::string callORput) const
{
	return annuity_*swaptionSABR_->premiumBlackScholes(strike,callORput);
};
double swaption::swapFairValue(double strike) const
{
	return annuity_*(swaptionSABR_->getForward()-strike);
};

double swaption::simulation(double corrRN)
{
	return swaptionSABR_->simulation(corrRN);
};

double swaption::getImpliedVol(double strike) const
{
	return swaptionSABR_->impliedVol(strike);		
};

double swaption::getParameterAlpha() const
{
	return swaptionSABR_->getParameterAlpha();		
};

double swaption::getParameterNu() const
{
	return swaptionSABR_->getParameterNu();		
};

double swaption::getParameterRho() const
{
	return swaptionSABR_->getParameterRho();		
};
}