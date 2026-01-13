#include <velesquant/models/cms.h>

#include <velesquant/local_vol/sabr.h>
#include <vector>
#include <string>
#include <velesquant/xlw/CellMatrix.h>
#include <algorithm>


using namespace std;
using namespace velesquant::xlw;

namespace velesquant {


cms::cms(double expirySR, double tenorSR, double forwardSR, double annuitySR, 
		 double payCMS, double discountCMS, 
		 double beta, double alpha, double nu, double rho)
	: discountCMS_(discountCMS)
{
	sabr* swaptionSABR = new sabr(expirySR,forwardSR,beta,alpha,nu,rho);
	double atmVol = swaptionSABR->impliedVol(forwardSR);
	delete swaptionSABR;
	double forwardCMS = convAdj(expirySR,forwardSR,discountCMS,annuitySR,tenorSR,atmVol);
	cmsSABR_ = new sabr(expirySR, forwardCMS, beta,alpha,nu,rho);
};

cms::cms(double expirySR, double tenorSR, double forwardSR, double annuitySR, 
		 double payCMS, double discountCMS, 
		 double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType)
	: discountCMS_(discountCMS)
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
	double forwardCMS = convAdj(expirySR,tenorSR,forwardSR,annuitySR,discountCMS,atmVol);
	cmsSABR_ = new sabr(expirySR, forwardCMS, beta,alpha,nu,rho);
};

// Added function by James Carey on 19-11-2015 to allow for alternative convexity adjustment calculation
cms::cms(double expirySR, double tenorSR, double forwardSR, double freqSR, double freqCMS,
		 double payCMS, double discountCMS, 
		 double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType)
	: discountCMS_(discountCMS),
	  freqSR_(freqSR),
	  freqCMS_(freqCMS)
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

	double forwardCMS = convAdjAlt(expirySR, tenorSR, forwardSR, freqSR, freqCMS, atmVol);
	cmsSABR_ = new sabr(expirySR, forwardCMS, beta,alpha,nu,rho);
};

cms::cms(double expirySR, double tenorSR, double forwardSR, double annuitySR, 
		double payCMS, double discountCMS, 
		double beta, CellMatrix& strikes, CellMatrix& marketQuotes, string quoteType,
		CellMatrix& initialParams)
	: discountCMS_(discountCMS)
{
	sabr* swaptionSABR = new sabr(expirySR,forwardSR,beta);
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
	//swaptionSABR->calibratorWithInitial(theStrikes,theQuotes,quoteType,theParams);
	swaptionSABR->calibratorWithInitial(theStrikes,theQuotes,quoteType);
	double alpha = swaptionSABR->getParameterAlpha();
	double nu = swaptionSABR->getParameterNu();
	double rho = swaptionSABR->getParameterRho();
	double atmVol = swaptionSABR->impliedVol(forwardSR);
	delete swaptionSABR;
	double forwardCMS = convAdj(expirySR,tenorSR,forwardSR,annuitySR,discountCMS,atmVol);
	cmsSABR_ = new sabr(expirySR, forwardCMS, beta,alpha,nu,rho);
};

double cms::fairValue(double strike, std::string callORput) const
{
	return discountCMS_*cmsSABR_->premiumBlackScholes(strike,callORput);
};
    
double cms::simulation(double corrRN)
{
	return cmsSABR_->simulation(corrRN);
};

double cms::getForward() const 
{
	return cmsSABR_->getForward();
};

double cms::getDiscountCMS() const 
{	
	return discountCMS_;	
};

double cms::getMaturity() const
{
	return cmsSABR_->getMaturity();		
};

double cms::getATMvol() const
{
	double forwardCMS = cmsSABR_->getForward();
	return cmsSABR_->impliedVol(forwardCMS);		
};

double cms::getImpliedVol(double strike) const
{
	return cmsSABR_->impliedVol(strike);		
};

double cms::getParameterAlpha() const
{
	return cmsSABR_->getParameterAlpha();		
};

double cms::getParameterNu() const
{
	return cmsSABR_->getParameterNu();		
};

double cms::getParameterRho() const
{
	return cmsSABR_->getParameterRho();		
};

double cms::convAdj(double expirySR, double tenorSR, double forwardSR,
		            double annuitySR, double discountCMS, double atmVol) const
{
	double A = 1.0/tenorSR;
	double B = (discountCMS/annuitySR - A)/forwardSR;
	if (atmVol*atmVol*expirySR < 1.0 * atan(1.0)) 
	return forwardSR*(A+B*forwardSR*exp(atmVol*atmVol*expirySR))/(A+B*forwardSR);
	else
	return forwardSR*(A+B*forwardSR*(1.0+atmVol*atmVol*expirySR))/(A+B*forwardSR); 
};

// Added function by James Carey on 19-11-2015 to allow for alternative convexity adjustment calculation
double cms::convAdjAlt(double expirySR, double tenorSR, double forwardSR,
					   double freqSR, double freqCMS, double atmVol) const
{
	double n = tenorSR*freqSR;
	double S0f = forwardSR/freqSR;
	double omega = 1 - S0f*(freqSR/freqCMS + n/(pow(1.0 + S0f,n) - 1.0))/(1.0 + S0f);
	double convAdj = forwardSR*omega*atmVol*atmVol*expirySR;
	return forwardSR + convAdj;
};
}