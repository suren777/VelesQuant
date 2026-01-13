#include <velesquant/models/HW.h>
#include <velesquant/models/utility.h>
#include <velesquant/local_vol/lm.h>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions.hpp>
#include <cmath>
#include <ql/quantlib.hpp>
#include <algorithm>

const double DT = 0.001;
const double SDT = sqrt(DT);
using namespace std;
using namespace QuantLib;
using namespace boost::placeholders;
#pragma warning (disable:4715)

namespace velesquant {



double HW::optionBond(double Expiry, double Maturity, double Strike, string callORput)
	{
	double v0 = totalVariance(Expiry); 
	double impVol = (1-exp(-kappa_*(Maturity-Expiry)))/kappa_*sqrt(v0/Expiry); 
	return formulaBlack(getDF(Expiry),getDF(Maturity),Strike,impVol,Expiry,callORput);
	};

double HW::ZC(double Expiry)
{
	double B = (1-exp(-kappa_*Expiry))/kappa_;
	return getDF(Expiry)*exp(-B*(whichFwdRate(0)+.5*B*totalVariance(Expiry)));
}

double HW::swaption(double Expiry, double Tenor, double Strike, double PayFrequency)
	{
	double CP = criticalPoint(Expiry,Tenor,Strike,PayFrequency);
	double vO = totalVariance(Expiry);
	double dfT0 = getDF(Expiry);
	int nC = int(Tenor/PayFrequency+0.5);
	double swaptionV = 0.0;
	for(int i=1; i<=nC; i++)
		{
		double Ti = Expiry+i*PayFrequency;
		double dfTi = getDF(Ti);
		double Gi = (1-exp(-kappa_*(Ti-Expiry)))/kappa_;
		double Ki = dfTi/dfT0*exp(-CP*Gi-0.5*vO*Gi*Gi);
		double impVol = (1-exp(-kappa_*(Ti-Expiry)))/kappa_*sqrt(vO/Expiry); 
		double Puti = formulaBlack(dfT0,dfTi,Ki,impVol,Expiry,"put");
		swaptionV += Strike*PayFrequency*Puti;
		if(i==nC) swaptionV += Puti;
		}
	return swaptionV;
	};

vector<double> HW::simulationHW(vector<double> times) const
	{
	int N=times.size();
	vector<double> path(N);
	int T = 0;
	double var = 0.0;
	double x = 0.0;
	for(int i=0; i<=int(times[N-1]/DT+0.5); i++)
		{
		if(i*DT <= times[T] && (i+1)*DT > times[T])
			{
			path[T] = x + whichFwdRate(i*DT);
			T++;
			}
		double sigma = whichSigma(i*DT);
		var += sigma*sigma*(exp(2*kappa_*i*DT)-exp(2*kappa_*(i-1)*DT))/(2*kappa_);
		x += (exp(-2*kappa_*i*DT)*var-kappa_*x)*DT + sigma*SDT*random_normal();  
		}
	return path;
	};

void HW::calibrator(std::vector<defSwap> swapQuotes, const string &type="swap")
	{
	quoteSwap_ = swapQuotes;
	double pos;
	timeSigmas_.resize(0);
	timeSigmas_.push_back( quoteSwap_[0].Expiry );
	pos=timeSigmas_[0];
	for (int i=1; i < int(quoteSwap_.size());i++)
		{
		if (quoteSwap_[i].Expiry>pos)
			{
			timeSigmas_.push_back( quoteSwap_[i].Expiry );
			pos = quoteSwap_[i].Expiry;
			}
		}
	sigmas_.resize( timeSigmas_.size());
	int n = sigmas_.size()+1;	    //no. of HW model volatility term paremeters 
	double* x = new double[n];		//initial estimate of parameters vector
	double* lb = new double[n];		//Constrains on parameters lb - lower and lb - upper
	double* ub = new double[n];
	for(int i=0; i<n-1; i++)
		{
		x[i] = sigmas_[i];		// sigmas initial value
		lb[i] = 1e-6;
		ub[i] = 1;
		}
	x[n-1] = kappa0_;			// kappa initial value
	lb[n-1] = 1e-6;
	ub[n-1] = 1;
	int m = quoteSwap_.size();		//no. of observations
	double* fvec = new double[m];	//no need to populate 
	double ftol = 1e-10; //tolerance
	double xtol = 1e-10; //tolerance
	double gtol = 1e-10; //tolerance
	int maxfev = 5000;  //maximum function evaluations
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

	lmfcn fcn;

	marketSwaption_.resize(m);
	if (type=="swap")
		{
		for(int i = 0; i < m; i++) 
			marketSwaption_[i] = swaptionATM(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].VolATM);
		fcn = boost::bind(&HW::objFcnPrice, this, _1,_2,_3,_4,_5,lb,ub);
		}
	else
		{
		for(int i = 0; i < m; i++) 
			marketSwaption_[i] = quoteSwap_[i].VolATM;
		fcn = boost::bind(&HW::objFcnIV, this, _1,_2,_3,_4,_5,lb,ub);
		}

	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);
	//the below is output result
	for(int i=0; i<n-1; i++)
		sigmas_[i] = fabs(x[i]);	// sigmas final value
	kappa_ = x[n-1];				// kappa final value
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

//void HW::objFcnPrice(int m, int n, double* x, double* fvec, int* iflag)
//	{	// int n = sigmas_.size()+1;
//	for(int i=0; i<n-1; i++)
//		sigmas_[i] = fabs(x[i]);	// sigmas
//	kappa_ = x[n-1];				// kappa
//	for(int i=0; i<m; i++)
//		{
//		double modelSwaption = swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
//		fvec[i] = modelSwaption - marketSwaption_[i];
//		}
//	};
double HW::pen_fun(double*x, double*lb, double*ub, int n)
	{
	double penalty = 0.0;
	double lambda = 1e6;
	for(int i=0; i<n; i++)penalty+=((x[i]<lb[i])?lambda*pow(x[i],2):0)+((x[i]>ub[i]))?lambda*pow(x[i],2):0;
	return penalty;
	}

void HW::objFcnPrice(int m, int n, double* x, double* fvec, int* iflag, double* lb, double* ub)
	{
	double penalty = pen_fun(x, lb, ub, n);
	for(int i=0; i<n-1; i++)
		{
		sigmas_[i] = x[i];		
		}
	//	sigmas_[i] = fabs(x[i]);	// sigmas
	kappa_ = x[n-1];				// kappa
	for(int i=0; i<m; i++)
		{
		double modelSwaption = swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
		fvec[i] = modelSwaption-marketSwaption_[i]+penalty;
		}
	};

void HW::objFcnIV(int m, int n, double* x, double* fvec, int* iflag, double* lb, double* ub)
	{	// int n = sigmas_.size()+1;
	double penalty = pen_fun(x, lb, ub, n);
	for(int i=0; i<n-1; i++)
		sigmas_[i] =x[i];	
	kappa_ = x[n-1];				
	for(int i=0; i<m; i++)
		{
		double modelSwaption = swaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
		fvec[i] = swaptionIVblack(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, modelSwaption) - marketSwaption_[i]+penalty;
		}
	};

double HW::totalVariance(double T0)
	{
	double var = 0.0;
	int nP = timeSigmas_.size();
	for(int n=0; n<nP; n++)
		{
		double ti = 0.0;
		if(n>0) ti = timeSigmas_[n-1];
		double te = timeSigmas_[n];
		if(T0>te) 
			var += sigmas_[n]*sigmas_[n]*(exp(2*kappa_*te)-exp(2*kappa_*ti))/(2*kappa_);
		else
			{	
			var += sigmas_[n]*sigmas_[n]*(exp(2*kappa_*T0)-exp(2*kappa_*ti))/(2*kappa_);
			break;
			}
		}

	if(T0>timeSigmas_[nP-1]) 
		var += sigmas_[nP-1]*sigmas_[nP-1]*(exp(2*kappa_*T0)-exp(2*kappa_*timeSigmas_[nP-1]))/(2*kappa_);

	return exp(-2*kappa_*T0)*var;
	};

double HW::equalityCP(double x, double vO, double Expiry, double Tenor, double Strike, double PayFrequency)
	{
	double dfT0 = getDF(Expiry);
	int nC = int(Tenor/PayFrequency+0.5);
	double CP = 1.0;
	for(int i=1; i<=nC; i++)
		{
		double Ti = Expiry+i*PayFrequency;
		double Gi = (1-exp(-kappa_*(Ti-Expiry)))/kappa_;
		double Ki = getDF(Ti)/dfT0*exp(-x*Gi-0.5*vO*Gi*Gi);
		CP -= Strike*PayFrequency*Ki;
		if(i==nC) CP -= Ki;
		}
	return CP;
	};

double HW::criticalPoint(double Expiry, double Tenor, double Strike, double PayFrequency)
	{
	double vO = totalVariance(Expiry);
	double xl = -0.5;
	double xu = 0.5;
	double x = 0.5*(xl+xu);
	double b = equalityCP(x, vO, Expiry, Tenor, Strike, PayFrequency);
	int Niter = 0;
	do {
		Niter++;      
		if (b < 0.0) xl = x;
		if (b > 0.0) xu = x;
		x = 0.5*(xl+xu);
		b = equalityCP(x, vO, Expiry, Tenor, Strike, PayFrequency);
		} while ( fabs(b) >= 1.0E-10 && Niter < 20 && (xu-xl)>=10e-6 );
	return x;
	};

double HW::formulaBlack(double dfT0, double dfTN, double Strike, double impVol, double T0, string callORput)
	{
	double sd = impVol*sqrt(T0);
	double d1 = log(dfTN/(dfT0*Strike))/sd + 0.5*sd;
	double d2 = d1-sd;
	double price;
	if(callORput=="call" || callORput=="Call" || callORput=="c" || callORput=="C")
		price = dfTN*cdf_normal(d1)-(dfT0*Strike)*cdf_normal(d2);
	else if(callORput=="put" || callORput=="Put" || callORput=="p" || callORput=="P")
		price = (dfT0*Strike)*cdf_normal(-d2)-dfTN*cdf_normal(-d1);
	return price;
	};

double HW::whichFwdRate(double t) const
	{
	int n = timeDFs_.size();
	if(t < timeDFs_[0])
		return -log(DFs_[0])/timeDFs_[0];
	if(t >= timeDFs_[n-1])
		return -log(DFs_[n-1]/DFs_[n-2])/(timeDFs_[n-1]-timeDFs_[n-2]);
	for(int i=1; i<n; i++)
		if((t >= timeDFs_[i-1]) && (t < timeDFs_[i]))
			return -log(DFs_[i]/DFs_[i-1])/(timeDFs_[i]-timeDFs_[i-1]);
	};
double HW::whichSigma(double t) const
	{
	int n = sigmas_.size();
	if(t < timeSigmas_[0])
		return sigmas_[0];
	if(t >= timeSigmas_[n-1])
		return sigmas_[n-1];
	for(int i=1; i<n; i++)
		if((t >= timeSigmas_[i-1]) && (t < timeSigmas_[i]))
			return sigmas_[i];
	};


// not used
double HW::swap(double Expiry, double Tenor, double Strike, double PayFrequency)
	{
	double dfT0 = getDF(Expiry);
	int nC = int(Tenor/PayFrequency+0.5);
	for(int i=1; i<=nC; i++)
		{
		double dfTi = getDF(Expiry+i*PayFrequency);
		dfT0 -= Strike*PayFrequency*dfTi;
		if(i==nC) dfT0 -= dfTi;
		}
	return dfT0;
	};
double HW::getSwapRate(double Expiry, double Tenor, double PayFrequency)
	{
	double dfT0 = getDF(Expiry);
	int nC = int(Tenor/PayFrequency+0.5);
	double Level = 0.0;
	for(int i=1; i<=nC; i++)
		{
		double dfTi = getDF(Expiry+i*PayFrequency);
		Level += PayFrequency*dfTi;
		if(i==nC) dfT0 -= dfTi;
		}
	return dfT0/Level;
	};

void HW::calibratorBstrp(std::vector<defSwap> swapQuotes, const string &type)
	{
	quoteSwap_ = swapQuotes;

	int MAX_ITER = 20;
	double X_TOL = 1e-8;
	double F_TOL  = 1e-12;
	double dF, F;
	double x=sigmas0_[0];
	int k = 0;
	boost::function<double (double)> objFcnB;
	double min = 1e-6;
	double max = 1;

	calibrateKappa();
	double pos;
	timeSigmas_.resize(0);
	timeSigmas_.push_back( quoteSwap_[0].Expiry );
	pos=timeSigmas_[0];
	for (int i=1; i < int(quoteSwap_.size());i++)
		{
		if (quoteSwap_[i].Expiry>pos)
			{
			timeSigmas_.push_back( quoteSwap_[i].Expiry );
			pos = quoteSwap_[i].Expiry;
			}
		}
	sigmas_.resize( timeSigmas_.size());
	for(int j = 0; j < int( timeSigmas_.size() );j++)
		{		
		iter_ = j;

		qSB_.resize(0);
		for (int l = k; l < int (quoteSwap_.size()) ; l++,k++)
			{			
			if(l <= int(quoteSwap_.size()) )
				{
				if (!(quoteSwap_[l].Expiry <= timeSigmas_[j])) {break;}
				else qSB_.push_back(quoteSwap_[l]);
				}
			}
		int m = qSB_.size();		//no. of observations
		marketSwaption_.resize(m);
		if (type=="swap")
			{
			for (int i = 0; i < m; i++) marketSwaption_[i] = swaptionATM(qSB_[i].Expiry, qSB_[i].Tenor, qSB_[i].VolATM);
			objFcnB = boost::bind(&HW::objFcnBswap, this,  _1, m);
			}
		else
			{
			for (int i = 0; i < m; i++) marketSwaption_[i] = qSB_[i].VolATM;
			objFcnB = boost::bind(&HW::objFcnBIV, this, _1, m);
			}
		double delta;
		for (int i = 0; i < MAX_ITER ; i++)
			{
			F = objFcnB(x);
			double Fu = objFcnB(x+X_TOL);
			double Fd = objFcnB(x-X_TOL);
			dF = 0.5*(Fu-Fd)/X_TOL;
			if ((fabs(F)<F_TOL)||(fabs(dF)<F_TOL)) 
				break;
			double x1 = x;
			delta = F/dF;
			x1 -= delta;
			if (x1 <= min )
				{
				delta = 0.5F * (x - min);
				x1 = x - delta;
				if((x1 == min) || (x1 == max))
					break;
				}
			else if (x1 >= max )
				{
				delta = 0.5F * (x - max);
				x1 = x - delta;
				if((x1 == min) || (x1 == max))
					break;
				}
			if(delta > 0)
				max = x;
			else
				min = x;
			if (abs(x1-x) < X_TOL*fabs(x1)||(x1!=x1)) break;
			x = x1;
			if ( fabs(max - min) < 1e-7 ) break;
			}
		sigmas_[j] = x;
		min = 1e-6;
		max = 1.0;
		}		
	};

double HW::objFcnBIV(double x, int m)
	{	
	double sum = 0.0;
	sigmas_[iter_] = x;	// sigmas
	for(int i=0; i<m; i++)
		{
		double modelSwaption = swaption(qSB_[i].Expiry, qSB_[i].Tenor, qSB_[i].SwapRate, qSB_[i].Frequency);
		double a = 1.0 - swaptionIVblack(qSB_[i].Expiry, qSB_[i].Tenor, modelSwaption)/marketSwaption_[i];
		sum+= fabs(a);
		}
	return sum/m*1e4;
	};
double HW::objFcnBswap(double x, int m)
	{	
	double sum = 0.0;
	sigmas_[iter_] = x;	// sigmas
	for(int i=0; i<m; i++)
		{
		double modelSwaption = swaption(qSB_[i].Expiry, qSB_[i].Tenor, qSB_[i].SwapRate, qSB_[i].Frequency);
		double a = modelSwaption/marketSwaption_[i]-1.0;
		sum+= a*a;
		}
	return sum/m;
	};

void HW::calibrateKappa()
	{
	int MAX_ITER = 100;
	double X_TOL = 1e-7;
	double F_TOL  = 1e-12;
	double a0;
	double a1;
	double dF,F;
	double min = 1e-6;
	double max = 1.0;
	int N = quoteSwap_.size();
	vector<double> iv_ratio;
	vector<double> bond_ratio;
	vector<int> index;

	for (int i = 0 ; i < N-1; i++)
		{
		if (quoteSwap_[i].Expiry==quoteSwap_[i+1].Expiry)
			{
			double ex = quoteSwap_[i].Expiry;
			double t1 = quoteSwap_[i].Tenor + ex;
			double t2 = quoteSwap_[i+1].Tenor + ex;
			index.push_back(i);
			iv_ratio.push_back(quoteSwap_[i+1].VolATM/quoteSwap_[i].VolATM);
			bond_ratio.push_back(
				(getDF(ex) - getDF(t1)) / 
				(getDF(ex) - getDF(t2))
				);
			}		
		}

	//double best, min_val=1000;
	a0 = kappa0_;
	double delta;
	for (int i = 0; i < MAX_ITER ; i++)
		{
		F = objKappa(a0,iv_ratio,bond_ratio,index);
		double Fu = objKappa(a0+X_TOL,iv_ratio,bond_ratio,index);
		double Fd = objKappa(a0-X_TOL,iv_ratio,bond_ratio,index);
		dF = 0.5*(Fu-Fd)/X_TOL;
		if ((fabs(F)<F_TOL)||(fabs(dF)<F_TOL)) 
			break;
		a1 = a0;
		delta = F/dF;
		a1 -= delta;
		if (a1 <= min)
			{
			delta = 0.5F * (a0 - min);
			a1 = a0 - delta;
			if((a1 == min) || (a1 == max))
				break;
			}
		else if (a1 >= max )
			{
			delta = 0.5F * (a0 - max);
			a1 = a0 - delta;
			if((a1 == min) || (a1 == max))
				break;
			}
		if(delta > 0)
			max = a0;
		else
			min = a0;
		if (abs(a1-a0) < X_TOL*fabs(a1)||(a1!=a1)) break;
		a0 = a1;
		if ( fabs(max - min) < 1e-5 ) break;
		}
	QL_ENSURE(a0==a0, "Failed to calibrate Kappa \n");
	kappa_ = a0;
	};

double HW::Bratio(double a, double M, double t1, double t2)
	{
	return (1.0 - exp(-a*(t2-M)))/(1.0 - exp(-a*(t1-M)));
	};



double HW::objKappa(double x, const vector<double> &ivr, const vector<double> &br, const vector<int> &ind)
	{
	int N = ivr.size();
	double sum = 0.0;
	if (x==0) return 1000;
	for (int i = 0; i< N; i++)
		{
		double ex = quoteSwap_[ind[i]].Expiry;
		double t1 = quoteSwap_[ind[i]].Tenor + ex;
		double t2 = quoteSwap_[ind[i]+1].Tenor + ex;
		double a = br[i]*Bratio(x, ex, t1, t2)-ivr[i];
		sum += a*a;
		}
	return sum/N;
	};

double HW::swaptionIVblack(double Expiry, double Tenor, double swaption_price)
	{
	double xl = 1e-6;
	double xu = 1;
	double x = 0.5*(xl+xu);

	double DF = getDF(Expiry)-getDF(Expiry+Tenor);
	double b = DF*(2.0*cdf_normal(0.5*x*sqrt(Expiry))-1.0)-swaption_price;	 

	int Niter = 0;
	do {
		Niter++;      
		if (b < 0.0) xl = x;
		if (b > 0.0) xu = x;
		x = 0.5*(xl+xu);
		b = DF*(2.0*cdf_normal(0.5*x*sqrt(Expiry))-1.0)-swaption_price;
		} while ( (fabs(b) >= 1.0E-12) && (Niter < 40) && (xu-xl > 1E-7) );
	return x;
	};

double HW::BlackStrike(double T0, double TN, double impVol)
	{
	double xl = 1e-6;
	double xu = 1;
	double x = 0.5*(xl+xu);	
	double dfT0 = getDF(T0);
	double dfTN = getDF(T0+TN);
	double DF = dfT0-dfTN;
	double atm_swap = DF*(2.0*cdf_normal(0.5*impVol*sqrt(T0))-1.0);
	double b = formulaBlack(dfT0, dfTN, x, impVol, T0, "call")-atm_swap;	 

	int Niter = 0;
	do {
		Niter++;      
		if (b < 0.0) xl = x;
		if (b > 0.0) xu = x;
		x = 0.5*(xl+xu);
		b = formulaBlack(dfT0, dfTN, x, impVol, T0, "call")-atm_swap;
		} while ( (fabs(b) >= 1.0E-12) && (Niter < 40) && (xu-xl > 1E-7) );
	return x;
	};
double HW::BlackStrikePlain(double T0, double TN)
	{
	double dfT0 = getDF(T0);
	double dfTN = getDF(T0+TN);
	return dfTN/dfT0;
	};

double HW::swaptionIVblackPub(double Expiry, double Tenor, double swap_price)
	{
	return swaptionIVblack(Expiry, Tenor, swap_price);
	};
}