
#include <velesquant/local_vol/sabr.h>
#include <velesquant/local_vol/lm.h>
#include <velesquant/models/utility.h>
#include <assert.h>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <ql/quantlib.hpp>
#include <velesquant/models/solver.h>
#include <algorithm>


using namespace std;
using namespace QuantLib;
using namespace boost::placeholders;

namespace velesquant {


sabr::sabr(double maturity, double forward, double beta, double alpha, double nu, double rho)
	: maturity_(maturity), forward_(forward), beta_(beta), alpha_(alpha), nu_(nu), rho_(rho) 
{ 
	parameters();
	shift_=0;
	notQUTL_ = true;
	type_="BS";
};

sabr::sabr(double maturity, double forward, double beta, double alpha, double nu, double rho, double shift)
	: maturity_(maturity), forward_(forward), beta_(beta), alpha_(alpha), nu_(nu), rho_(rho) , shift_(shift)
{ 
	parameters();
	notQUTL_ = true;
	type_="Normal";
};

double sabr::normalVol(double K) const
{
	if(beta_==0.0)
		if (K==forward_)
			return (Nvolb0(K*.99999)+Nvolb0(K*1.00001))/2;
		else
			return Nvolb0(K);
	else if (beta_ == 1.0) 
		if (K==forward_)
			return (Nvolb1(K*.99999)+Nvolb1(K*1.00001))/2;
		else
			return Nvolb1(K);
	else
		if (K==forward_)
			return (Nvolb(K*.99999)+Nvolb(K*1.00001))/2;
		else
			return Nvolb(K);
}

double sabr::Nvolb0(double K) const
{
	double zeta = nu_/alpha_*(forward_-K);
	double x = log((sqrt(1-rho_*zeta+zeta*zeta)-rho_+zeta)/(1-rho_));
	return alpha_*zeta/x*(1+(.25*rho_*nu_*alpha_+1/24.0*(2-3*rho_*rho_)*nu_*nu_)*maturity_);
};

double sabr::Nvolb(double K) const
{
	double fbk = pow(forward_+shift_,beta_)-pow(K+shift_,beta_);
	double fmbk = pow(forward_+shift_,1-beta_)-pow(K+shift_,1-beta_);
	double lfk = log((forward_+shift_)/(K+shift_));
	double z = nu_/alpha_*fmbk/(1-beta_);
	double x = log((sqrt(1-rho_*z+z*z)-rho_+z)/(1-rho_));
	double p1 = alpha_*(forward_-K)*(1-beta_)/fmbk*z/x;
	double p2 = -1.0/24*(beta_*(2-beta_)*(1-beta_)*(1-beta_)*alpha_*alpha_*lfk*lfk)/(fmbk*fmbk);
	double p3 = 0.25*rho_*nu_*alpha_*fbk/(forward_-K)+(2-3*rho_*rho_)/24*nu_*nu_;
	return p1*(1+(p2+p3)*maturity_);
};

double sabr::Nvolb1(double K) const
{
	double lfbk = log(forward_+shift_) + log(K+shift_);
	double z = nu_/alpha_*lfbk;
	double x = log((sqrt(1-rho_*z+z*z)-rho_+z)/(1-rho_));
	double p1 = alpha_*(forward_-K)/lfbk*z/x;
	double p3 = 0.25*rho_*nu_*alpha_+(2-3*rho_*rho_)/24*nu_*nu_-1.0/24*alpha_*alpha_;
	return p1*(1.0+p3*maturity_);
};

double sabr::getVol(double strike) const
{
	if(type_ == "BS") return impliedVol(strike);
	else return normalVol(strike);
};
double sabr::getPremium(double strike, string callOrPut) const
{
	if(type_ == "BS") return premiumBlackScholes(strike,callOrPut);
	else return premiumBachelier(strike,callOrPut);
};
double sabr::impliedVol(double strike) const
{
	if(!calibrated_) 
		throw ("SABR model is not calibrated yet!");
	double mean = std::pow(forward_*strike,0.5*(1.0-beta_));
	double relR = std::log(forward_/strike);
	double tm = (1.0-beta_)*(1.0-beta_)*relR*relR;
	double dem = mean*(1.0+tm/24.0+tm*tm/1920.0);
	double coe = (1.0-beta_)*alpha_/mean*(1.0-beta_)*alpha_/mean/24.0
		+rho_*beta_*nu_*alpha_/mean/4.0 +(2.0-3.0*rho_*rho_)*nu_*nu_/24.0;
	double z = nu_/alpha_*mean*relR;
	double mltr;
	if(std::fabs(z*z)>1.0E-20) 
	{
		double xz = std::log( (std::sqrt(1.0-2.0*rho_*z+z*z)+z-rho_)/(1.0-rho_) );
		mltr = z/xz;
	} 
	else 
		mltr = 1.0 -0.5*rho_*z +(2.0-3.0*rho_*rho_)*z*z/12.0;
	double impVol = alpha_/dem*mltr*(1.0+maturity_*coe);
	//assert(impVol>0.0);
	if(impVol < 0.0) impVol=1.0E-09;
	return impVol;
}

double sabr::premiumBachelier(double strike, std::string callORput) const
{
	int I = (callORput=="Call")?1:-1;
	double vol = normalVol(strike)*std::sqrt(maturity_);
	double q = (forward_-strike)/vol;
	return I*(forward_-strike)*cdf_normal(I*q)+vol*pdf_normal(I*q);
};

double sabr::premiumBlackScholes(double strike, std::string callORput) const
{
	double vol = impliedVol(strike) * std::sqrt(maturity_);
	double d1 = std::log(forward_/strike)/vol +0.5*vol;
	double d2 = d1-vol;
	double premium;
	if(forward_ >= strike)
	{
		if((d1 != d1) || (d2 != d2))
			premium = forward_ - strike;
		else
			premium = forward_*cdf_normal(d1) - strike*cdf_normal(d2);
		if(callORput=="put" || callORput=="Put" || callORput=="p" || callORput=="P") 
			premium -= (forward_-strike);
	}
	else
	{
		if((d1 != d1) || (d2 != d2))
			premium = 0.0;
		else
			premium = strike*cdf_normal(-d2) - forward_*cdf_normal(-d1);
		if( !(callORput=="put" || callORput=="Put" || callORput=="p" || callORput=="P") ) 
			premium += (forward_-strike);
	}
	return premium;
};

double sabr::localVol(double spot) const
{	
	double lv = localVolCall(spot);
	//double lv = localVolIV(spot);
	//double lv = localVolzabr(spot);
	return lv;
};

double sabr::localVolCall(double spot) const   //based on the call premium
{
	double cLft = premiumBlackScholes(spot*0.999);
	double cVal = premiumBlackScholes(spot);
	double cRgt = premiumBlackScholes(spot*1.001);
	double dS2val = (cRgt+cLft-2.0*cVal)/(0.001*spot*0.001*spot);
	if(dS2val<1.0E-10) 
		dS2val = 1.0E-10;
	assert(dS2val>0.0);
	sabr *incrT = new sabr(maturity_*1.001,forward_,beta_,alpha_,nu_,rho_);
	double cBigT = incrT->premiumBlackScholes(spot);
	delete incrT;
	double dTval = (cBigT-cVal)/(0.001*maturity_);
	if(dTval<1.0E-10) 
		dTval = 1.0E-10;
	assert(dTval>0.0);
	double lv = std::sqrt(2.0*dTval/dS2val)/spot;
	//if(lv != lv) lv = 0.0;	if(lv > 111.111) lv = 111.111;
	return lv;
};

double sabr::localVolIV(double spot) const   //based on the implied vol
{
	double volL = impliedVol(spot*0.999);
	double vol = impliedVol(spot);
	double volR = impliedVol(spot*1.001);
	double dSvol = (volR-volL)/(2.0*0.001*spot);
	double dS2vol = (volR+volL-2.0*vol)/(0.001*spot*0.001*spot);
	sabr *incrT = new sabr(maturity_*1.001,forward_,beta_,alpha_,nu_,rho_);
	double volT = incrT->impliedVol(spot);
	delete incrT;
	double dTvol = (volT-vol)/(0.001*maturity_);
	double derivation = vol*sqrt(maturity_);
	double d1 = log(forward_/spot)/derivation +0.5*derivation;
	double d2 = d1-derivation;
	double fz = vol*vol +2.0*vol*maturity_*dTvol;
	double fm = 1.0 +2.0*d1*spot*sqrt(maturity_)*dSvol +spot*spot*maturity_*(d1*d2*dSvol*dSvol+vol*dS2vol);
	double lv = sqrt(fz/fm);
	return lv;
};

double sabr::localVolzabr(double spot) const
{	
	double y = nu_/alpha_ * (std::pow(forward_,(1.0-beta_))-std::pow(spot,(1.0-beta_)))/(1.0-beta_);	
	double jy = std::sqrt(1.0-2.0*rho_*y+y*y);	
	//double x = std::log((jy+y-rho_)/(1-rho_))/nu_;
	//double xi = std::abs(x)/std::sqrt(maturity_);
	//double px = 2.0 * (1.0 - xi * cdf_normal(-xi) / pdf_normal(xi));
	//double lv = jy * alpha_ * std::pow(spot,beta_)/spot * std::sqrt(px);
	double lv = jy * alpha_ * std::pow(spot,beta_)/spot;
	//if(lv != lv) lv = 0.0;	if(lv > 5.0) lv = 5.0;
	return lv;
};


void sabr::objFcn(int m, int n, double* x, double* fvec, int* iflag, std::string quoteType)
{
	setParameterAlpha(x[0]);
	setParameterNu(x[1]);
	setParameterRho(x[2]);
	double penalty = ((alpha_<=0)?1e32*alpha_*alpha_:0)+((nu_<=0)?1e32*nu_*nu_:0);
	for(int i=0; i<m; i++)
	{
		if(quoteType == "premium")
			//fvec[i] = (premiumBlackScholes(strikes_[i])-marketQuotes_[i])/(0.1+abs(1.0-strikes_[i]/forward_));
				//	fvec[i] = (premiumBlackScholes(strikes_[i])-marketQuotes_[i]);
					fvec[i] = (getPremium(strikes_[i])-marketQuotes_[i])+penalty;
		else
			//fvec[i] = (impliedVol(strikes_[i])-marketQuotes_[i])/(0.1+abs(1.0-strikes_[i]/forward_));
			fvec[i] = (getVol(strikes_[i])-marketQuotes_[i])+penalty;
	}
};

void sabr::objFcnATM(int m, int n, double* x, double* fvec, int* iflag, std::string quoteType)
{

	setParameterNu(x[0]);
	setParameterRho(x[1]);

	if (beta_ == 1)
	{
		setParameterAlpha(ATMvolRoots(0.0, 1.0, 0.0001));
	}
	else
	{
		setParameterAlpha(ATMvolRoots(0.0, 100.00, 0.0001));
	}
	double penalty = ((nu_<=0)?1e8*nu_*nu_:0);
	for(int i=0; i<m; i++)
	{
		if(quoteType == "premium")
			//fvec[i] = (premiumBlackScholes(strikes_[i])-marketQuotes_[i])/(0.1+abs(1.0-strikes_[i]/forward_));
				fvec[i] = (premiumBlackScholes(strikes_[i])-marketQuotes_[i])+penalty;
		else
			//fvec[i] = (impliedVol(strikes_[i])-marketQuotes_[i])/(0.1+abs(1.0-strikes_[i]/forward_));
			fvec[i] = (impliedVol(strikes_[i])-marketQuotes_[i])+penalty;
	}
};

void sabr::calibratorWithInitialATM(std::vector<double> strikes, std::vector<double> marketQuotes, std::string quoteType)
{
	int m = strikes.size();			//no. of observations
	strikes_.resize(m);		
	strikes_ = strikes; 	
	marketQuotes_.resize(m);		
	marketQuotes_ = marketQuotes;
	int n = 2;					//no. of SABR model paremeters 
	double* x = new double[n];	//initial estimate of parameters vector

	x[0] = nu_;	//nu
	x[1] = rho_;	//rho

	double* fvec = new double[m]; //no need to populate 
	double ftol = 1e-10; //tolerance
	double xtol = 1e-10; //tolerance
	double gtol = 1e-10; //tolerance
	int maxfev = 5000; //maximum function evaluations
	double epsfcn = 1e-10; //tolerance
	double* diag = new double[n]; //some internal thing
	int mode = 1; //some internal thing
	double factor = 1; // a default recommended value
	int nprint = 0; //don't know what it does
	int info = 0; //output variable
	int nfev = 0; //output variable will store no. of function evals
	double* fjac = new double[m*n]; //output array of jacobian
	int ldfjac = m; //recommended setting
	int* ipvt = new int[n]; //for internal use
	double* qtf = new double[n]; //for internal use
	double* wa1 = new double[n]; //for internal use
	double* wa2 = new double[n]; //for internal use
	double* wa3 = new double[n]; //for internal use
	double* wa4 = new double[m]; //for internal use

	boost::function<void (sabr*, int,int,double*,double*,int*,std::string)> obj = &sabr::objFcnATM;
	lmfcn fcn = boost::bind(obj, this, _1,_2,_3,_4,_5, quoteType);

	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);

	//*   info is an integer output variable. if the user has terminated execution, info is set to the (negative)
	//*     value of iflag. see description of fcn. otherwise, info is set as follows.
	//	*     info = 0  improper input parameters.
	//	*     info = 1  both actual and predicted relative reductions in the sum of squares are at most ftol.
	//	*     info = 2  relative error between two consecutive iterates is at most xtol.
	//	*     info = 3  conditions for info = 1 and info = 2 both hold.
	//	*     info = 4  the cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.
	//	*     info = 5  number of calls to fcn has reached or exceeded maxfev.
	//	*     info = 6  ftol is too small. no further reduction in the sum of squares is possible.
	//	*     info = 7  xtol is too small. no further improvement in the approximate solution x is possible.
	//	*     info = 8  gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.

	QL_ENSURE(info!=4, "SABR caliration fails with the info = " << info);

	//the below is output result
	setParameterNu(x[0]);
	setParameterRho(x[1]);
	if (beta_ == 1)
	{
		setParameterAlpha(ATMvolRoots(0.0, 1.0, 0.0001));
	}
	else
	{
		setParameterAlpha(ATMvolRoots(0.0, 100.00, 0.0001));
	}

	parameters();  
	delete[] x;
	delete[] fvec;
	delete[] diag;
	delete[] fjac;
	delete[] ipvt;
	delete[] qtf;
	delete[] wa1;
	delete[] wa2;
	delete[] wa3;
	delete[] wa4;
};


void sabr::calibratorWithInitial(std::vector<double> strikes, std::vector<double> marketQuotes, std::string quoteType)
{
	int m = strikes.size();			//no. of observations
	strikes_.resize(m);		
	strikes_ = strikes; 	
	marketQuotes_.resize(m);		
	marketQuotes_ = marketQuotes;
	int n = 3;					//no. of SABR model paremeters 
	double* x = new double[n];	//initial estimate of parameters vector
	//x[0] = initialParams[0];	//alpha
	//x[1] = initialParams[1];	//nu
	//x[2] = initialParams[2];	//rho

	x[0] = alpha_;	//alpha
	x[1] = nu_;	//nu
	x[2] = rho_;	//rho

	double* fvec = new double[m]; //no need to populate 
	double ftol = 1e-10; //tolerance
	double xtol = 1e-10; //tolerance
	double gtol = 1e-10; //tolerance
	int maxfev = 5000; //maximum function evaluations
	double epsfcn = 1e-10; //tolerance
	double* diag = new double[n]; //some internal thing
	int mode = 1; //some internal thing
	double factor = 1; // a default recommended value
	int nprint = 0; //don't know what it does
	int info = 0; //output variable
	int nfev = 0; //output variable will store no. of function evals
	double* fjac = new double[m*n]; //output array of jacobian
	int ldfjac = m; //recommended setting
	int* ipvt = new int[n]; //for internal use
	double* qtf = new double[n]; //for internal use
	double* wa1 = new double[n]; //for internal use
	double* wa2 = new double[n]; //for internal use
	double* wa3 = new double[n]; //for internal use
	double* wa4 = new double[m]; //for internal use

	boost::function<void (sabr*, int,int,double*,double*,int*,std::string)> obj = &sabr::objFcn;
	lmfcn fcn = boost::bind(obj, this, _1,_2,_3,_4,_5, quoteType);

	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);

	//*   info is an integer output variable. if the user has terminated execution, info is set to the (negative)
	//*     value of iflag. see description of fcn. otherwise, info is set as follows.
	//	*     info = 0  improper input parameters.
	//	*     info = 1  both actual and predicted relative reductions in the sum of squares are at most ftol.
	//	*     info = 2  relative error between two consecutive iterates is at most xtol.
	//	*     info = 3  conditions for info = 1 and info = 2 both hold.
	//	*     info = 4  the cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.
	//	*     info = 5  number of calls to fcn has reached or exceeded maxfev.
	//	*     info = 6  ftol is too small. no further reduction in the sum of squares is possible.
	//	*     info = 7  xtol is too small. no further improvement in the approximate solution x is possible.
	//	*     info = 8  gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.

	QL_ENSURE(info!=4, "SABR caliration fails with the info = " << info);

	//the below is output result
	setParameterAlpha(x[0]);
	setParameterNu(x[1]);
	setParameterRho(x[2]);
	parameters();  
	delete[] x;
	delete[] fvec;
	delete[] diag;
	delete[] fjac;
	delete[] ipvt;
	delete[] qtf;
	delete[] wa1;
	delete[] wa2;
	delete[] wa3;
	delete[] wa4;
};

void sabr::calibrator(std::vector<double> strikes, std::vector<double> marketQuotes, std::string quoteType) 
{
	int m = strikes.size();			//no. of observations
	strikes_.resize(m);		
	strikes_ = strikes; 	
	marketQuotes_.resize(m);		
	marketQuotes_ = marketQuotes;
	int n = 3;					//no. of SABR model paremeters 
	double* x = new double[n];	//initial estimate of parameters vector
	//x[0] = 0.05;  //alpha
	//x[1] = 0.65;  //nu
	//x[2] = -0.05; //rho
	x[0] = alpha_;
	x[1] = nu_;
	x[2] = rho_;
	//x[0] = 0.75;  //alpha
	//x[1] = 1.25;  //nu
	//x[2] = -0.75; //rho

	double* fvec = new double[m]; //no need to populate 
	double ftol = 1e-10; //tolerance
	double xtol = 1e-10; //tolerance
	double gtol = 1e-10; //tolerance
	int maxfev = 5000; //maximum function evaluations
	double epsfcn = 1e-10; //tolerance
	double* diag = new double[n]; //some internal thing
	int mode = 1; //some internal thing
	double factor = 1; // a default recommended value
	int nprint = 0; //don't know what it does
	int info = 0; //output variable
	int nfev = 0; //output variable will store no. of function evals
	double* fjac = new double[m*n]; //output array of jacobian
	int ldfjac = m; //recommended setting
	int* ipvt = new int[n]; //for internal use
	double* qtf = new double[n]; //for internal use
	double* wa1 = new double[n]; //for internal use
	double* wa2 = new double[n]; //for internal use
	double* wa3 = new double[n]; //for internal use
	double* wa4 = new double[m]; //for internal use

	boost::function<void (sabr*, int,int,double*,double*,int*,std::string)> obj = &sabr::objFcn;
	lmfcn fcn = boost::bind(obj, this, _1,_2,_3,_4,_5, quoteType);

	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);

	//*   info is an integer output variable. if the user has terminated execution, info is set to the (negative)
	//*     value of iflag. see description of fcn. otherwise, info is set as follows.
	//	*     info = 0  improper input parameters.
	//	*     info = 1  both actual and predicted relative reductions in the sum of squares are at most ftol.
	//	*     info = 2  relative error between two consecutive iterates is at most xtol.
	//	*     info = 3  conditions for info = 1 and info = 2 both hold.
	//	*     info = 4  the cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.
	//	*     info = 5  number of calls to fcn has reached or exceeded maxfev.
	//	*     info = 6  ftol is too small. no further reduction in the sum of squares is possible.
	//	*     info = 7  xtol is too small. no further improvement in the approximate solution x is possible.
	//	*     info = 8  gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.

	QL_ENSURE(info!=4, "SABR caliration fails with the info = " << info);

	//the below is output result
	setParameterAlpha(x[0]);
	setParameterNu(x[1]);
	setParameterRho(x[2]);
	parameters();  
	delete[] x;
	delete[] fvec;
	delete[] diag;
	delete[] fjac;
	delete[] ipvt;
	delete[] qtf;
	delete[] wa1;
	delete[] wa2;
	delete[] wa3;
	delete[] wa4;
};

//void sabr::calibratorNormalVol(std::vector<double> strikes, std::vector<double> marketQuotes)
//{
//	int n = 3;
//	int m = strikes.size();			//no. of observations
//	strikes_.resize(m);		
//	strikes_ = strikes; 	
//	marketQuotes_.resize(m);	
//	double* x = new double[n];		//initial estimate of parameters vector
//	x[0] = alpha_;
//	x[1] = rho_;
//	x[2] = nu_;
//
//	QL_ENSURE(m>=n, "too much freedom in Calibration  " << m-n);
//	double* fvec = new double[m];	//no need to populate 
//	double ftol = 1e-10; //tolerance
//	double xtol = 1e-10; //tolerance
//	double gtol = 1e-10; //tolerance
//	int maxfev = 10000;  //maximum function evaluations
//	double epsfcn = 1e-10; //tolerance
//	double* diag = new double[n]; //some internal thing
//	int mode = 1; //some internal thing
//	double factor = 1; // a default recommended value
//	int nprint = 0; //don't know what it does
//	int info = 0; //output variable
//	int nfev = 0; //output variable will store no. of function evals
//	double* fjac = new double[m*n]; //output array of jacobian
//	int ldfjac = m; //recommended setting
//	int* ipvt = new int[n]; //for internal use
//	double* qtf = new double[n]; //for internal use
//	double* wa1 = new double[n]; //for internal use
//	double* wa2 = new double[n]; //for internal use
//	double* wa3 = new double[n]; //for internal use
//	double* wa4 = new double[m]; //for internal use	
//	boost::function<void (sabr*, int,int,double*,double*,int*)> obj = &sabr::objFcnNormalVolCalibration;
//	lmfcn fcn = boost::bind(obj, this, _1,_2,_3,_4,_5);
//	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
//		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);
//
//	QL_ENSURE(info != 4, "SABR Model Calibration Fails " << info);
//	//the below is output result
//
//	alpha_ = x[0];
//	rho_ = sin(x[1]);
//	nu_ = x[2];
//
//	delete[] x;
//	delete[] fvec;
//	delete[] diag;
//	delete[] fjac;
//	delete[] ipvt;
//	delete[] qtf;
//	delete[] wa1;
//	delete[] wa2;
//	delete[] wa3;
//	delete[] wa4;
//};

//void sabr::objFcnNormalVolCalibration(int m, int n, double* x, double* fvec, int* iflag)
//{
//	alpha_ = x[0];
//	rho_ = sin(x[1]);
//	nu_ = x[2];
//	double penalty = ((alpha_<=0)?1e8*alpha_*alpha_:0)+((nu_<=0)?1e8*nu_*nu_:0);
//	double error;
//	for (int i = 0; i<m; i ++)
//	{
//		error = marketQuotes_[i]-normalVol(strikes_[i]);
//		fvec[i] = error+penalty;
//	}
//};




boost::math::normal_distribution<> stdnormal(0.0,1.0);
void sabr::qutlTable()
{
	int N = 1200;
	vector<double> spots(N),qutls,cdfss;
	if (forward_ < 1.0) spots = grMesh(N);
	else  spots = gsMesh(N);
	vector<double> isp;

	for(int i=0; i<N; i++)
	{
		double cVal = getPremium(spots[i]);
		double cRgt = getPremium(spots[i]*1.001);
		double cdf = 1.0+1000.0*(cRgt-cVal)/spots[i];
		if ((cdf > 0.0)&&(cdf<1.0))
		{
			qutls.push_back(boost::math::quantile(stdnormal,cdf));
			cdfss.push_back(cdf);
			isp.push_back(spots[i]);
		}
	}
	//sorting the cdf distribution in an increasing order
	vector<double> scdfs; 
	asort(qutls,isp,cdfss);
	double current = qutls[0];
	int k = isp.size();
	for(int i = 0; i < k; i++)
		if((qutls[i] > current))
		{	
			current = qutls[i];
			spots_.push_back(isp[i]);
			qutls_.push_back(qutls[i]);
			cdfs_.push_back(cdfss[i]);
		}
		//adjusting (rebiaseing) the cdf to meet the forward
		N = spots_.size();
		double mean = 0.0;
		for(int i=0; i<N; i++)
		{
			if(i == 0) 	mean += 0.5*spots_[i]*cdfs_[i];
			else		mean += 0.5*(spots_[i]+spots_[i-1])*(cdfs_[i]-cdfs_[i-1]);
		}
		double bias = forward_ - mean;
		for(int i=0; i<N; i++) spots_[i] += bias;

		notQUTL_ = false;
};

double sabr::simulation(double corrRN)
{
	if(notQUTL_)
		qutlTable();
	QuantLib::LinearInterpolation interp(qutls_.begin(),qutls_.end(),spots_.begin());
	double qmin = qutls_[0];
	double qmax = qutls_[qutls_.size()-1];
	bool allowExtrapolation = true;
	if (corrRN>qmax) corrRN=qmax; 
	else if(corrRN<qmin) corrRN=qmin;
	return interp((qmin>corrRN)?qmin:corrRN, allowExtrapolation);
};
std::vector<double> sabr::simulations(std::vector<double> correlatedRNs)
{
	if(notQUTL_)
		qutlTable();
	QuantLib::LogLinearInterpolation interp(qutls_.begin(),qutls_.end(),spots_.begin());
	int N = correlatedRNs.size();
	std::vector<double> spots(N);
	bool allowExtrapolation = true;
	double qmin = qutls_[0];
	double qmax = qutls_[qutls_.size()-1];
	for(int i=0; i<N; i++)
		spots[i] = interp((((qmin>correlatedRNs[i])?qmin:correlatedRNs[i])>qmax)?qmax:correlatedRNs[i], allowExtrapolation);
	return spots;
};

int sabr::amin(vector<double> q, int position=0)
{
	int N = q.size();
	double min_q = q[position];
	int mpos=position;
	for (int i = position+1; i < N; i++)
		if (min_q>q[i])
		{
			min_q = q[i]; 
			mpos = i;
		}
		return mpos;	
};

int sabr::amax(vector<double> q, int position=0)
{
	int N = q.size();
	double max_q = q[position];
	int mpos=position;
	for (int i = position+1; i < N; i++)
		if (max_q<q[i])
		{
			max_q = q[i]; 
			mpos = i;
		}
		return mpos;	
};

void sabr::asort(vector<double> &a, vector<double> &b, vector<double> &c)
{
	int N = a.size();
	int imin = amin(a);
	double dmin = a[imin];

	for (int i = 1; i < N; i++)
	{
		double auxa = a[imin];
		double auxb = b[imin];
		double auxc = c[imin];		
		a[imin] = a[i-1];
		b[imin] = b[i-1];
		c[imin] = c[i-1];
		a[i-1] = auxa;
		b[i-1] = auxb;
		c[i-1] = auxc;
		imin=amin(a,i);
	}
};

vector<double> sabr::gsMesh (int Ms)
{
	vector<double> Mesh(Ms);
	double Sleft = .7*forward_;
	double Sright = 1.4*forward_;
	double d1 = forward_/20.0;
	double Smax = forward_*14.0;
	double zmin = asinh(-Sleft/d1);
	double zint = (Sright-Sleft)/d1;
	double zmax = zint+asinh((Smax-Sright)/d1);
	double dz = (zmax-zmin)/Ms;
	for (int i = 0; i < Ms ; i++)
	{
		double z = zmin + i*dz;
		if ( z < 0 )
			Mesh[i] = Sleft + d1 * sinh(z);
		else  
			Mesh[i] = ( z <= zint )?(Sleft+d1*z):(Sright+d1*sinh(z-zint));		
	}
	return Mesh;
};



vector<double> sabr::grMesh (int Mv)
{
	vector<double> Mesh(Mv);
	double Rmax = 1.0;
	double d3 = Rmax/(Mv*.9);
	double start = asinh((-Rmax-forward_)/d3);
	double dx = 1.0/Mv*(asinh((Rmax - forward_)/d3)-start);
	for (int i = 0; i < Mv ; i++)
		Mesh[i] =forward_+d3*sinh(start + i*dx);
	return Mesh;
};

double sabr::ATMvolPoly(double alpha)
{
	double a = (pow(1-beta_,2)*maturity_)/(24*pow(forward_,2-2*beta_));
	double b = (rho_*beta_*nu_*maturity_)/(4*pow(forward_,1-beta_));
	double c = (1 + (2-3*pow(rho_,2))/(24)*nu_*nu_*maturity_);
	double d = atmVol_*pow(forward_,1-beta_);

	return a*pow(alpha,3) + b*pow(alpha,2) + c*alpha - d; //=0
};

double sabr::ATMvolRoots(double lBound, double uBound, double tol)
{
	//Bisection algorithm
	double alpha = 0.5*(lBound+uBound);
	double y = ATMvolPoly(alpha);
	int Step = 0;
	do
	{       
		if (y < 0.0) lBound = alpha;
		if (y > 0.0) uBound = alpha;
		alpha = 0.5*(lBound+uBound);
		y = ATMvolPoly(alpha);
		if ((uBound-lBound)<tol/ 10) Step++;
	} while ( (fabs(y-0.0) > tol)&&(Step < 50) );
	return alpha;
}
}