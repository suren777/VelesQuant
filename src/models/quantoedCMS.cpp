#include <velesquant/models/quantoedCMS.h>

#include <velesquant/local_vol/sabr.h>
#include <vector>
#include <string>
#include <velesquant/xlw/CellMatrix.h>
#include <algorithm>


using namespace std;
using namespace velesquant::xlw;

namespace velesquant {


quantoedCMS::quantoedCMS(double expirySR, double tenorSR, double forwardSR, double annuitySR, 
		 double payCMS, double discountCMS, double corFX, double atmVolFX, 
		 double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType)
{
	sabr* swaptionSABR = new sabr(expirySR,forwardSR,beta);
	int m = strikes.RowsInStructure();
	vector<double> theStrikes(m), theQuotes(m);
	for(int i=0; i<m; ++i) 
	{
		theStrikes[i] = strikes(i,0).NumericValue();
		theQuotes[i] = marketQuotes(i,0).NumericValue();
	}
	swaptionSABR->calibrator(theStrikes,theQuotes,quoteType);
	double alpha = swaptionSABR->getParameterAlpha();
	double nu = swaptionSABR->getParameterNu();
	double rho = swaptionSABR->getParameterRho();
	double atmVol = swaptionSABR->impliedVol(forwardSR);
	delete swaptionSABR;

	double quantoedForwardSR = quantoAdj(expirySR,forwardSR,atmVol,corFX,atmVolFX);
	double quantoeForwardCMS = convAdj(expirySR,tenorSR,forwardSR,quantoedForwardSR,annuitySR,discountCMS,atmVol);
	quantoedCMSSABR_ = new sabr(expirySR, quantoeForwardCMS, beta,alpha,nu,rho);
};    


double quantoedCMS::fairValue(double strike, std::string callORput) const
{
	return quantoedCMSSABR_->premiumBlackScholes(strike,callORput);
};

double quantoedCMS::simulation(double corrRN)
{
	return quantoedCMSSABR_->simulation(corrRN);
};

double quantoedCMS::getForward() const 
{
	return quantoedCMSSABR_->getForward();
};

double quantoedCMS::quantoAdj(double expirySR, double forwardSR, double atmVol, double corFX, double atmVolFX) const
{
	return forwardSR*exp(corFX*atmVolFX*atmVol*expirySR);
};

double quantoedCMS::convAdj(double expirySR, double tenorSR, double forwardSR, double quantoedForwardSR,
		                    double annuitySR, double discountCMS, double atmVol) const
{
	double A = 1.0/tenorSR;
	double B = (discountCMS/annuitySR - A)/forwardSR;
	return quantoedForwardSR*(A+B*quantoedForwardSR*exp(atmVol*atmVol*expirySR))/(A+B*quantoedForwardSR);
};
}