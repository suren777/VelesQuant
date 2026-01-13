
#include <velesquant/pde_solvers/HWPDE.h>
#include <velesquant/pde_solvers/CyclicReduction.h>
#include <velesquant/local_vol/lm.h>
#include <boost/function.hpp>
#include <boost/math/distributions.hpp>
#include <boost/bind.hpp>
#include <ql/errors.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <omp.h>
#include <algorithm>

using namespace std;
#pragma warning (disable:4715)
#pragma warning (disable:4018)

namespace velesquant {

#define MR  512 //int(std::pow(2.0,9));
static const int DAYS_A_YEAR = 10001;
static const double delT = 10.0/(DAYS_A_YEAR-1.0);
const int nt = omp_get_max_threads();


void HWPDE::buildGrid(double Rmax ,double factor)
{
	gridR_=grMesh(MR-1, Rmax, factor);
	a_.resize(MR);
	b_.resize(MR);
	c_.resize(MR);
	e_.resize(MR);
	f_.resize(MR);
	g_.resize(MR);

#pragma omp parallel for
	for (int i=0; i < MR ; i++)
	{

		double R = gridR_[i];
		double Ru = (i	==	MR-1)? vright_	: gridR_[i+1];
		double Rl = (i	==	0	)? vleft_	: gridR_[i-1];
		a_[i] = 2.0/(R-Rl)/(Ru-Rl);
		b_[i] = -2.0/(R-Rl)/(Ru-R);
		c_[i] = 2.0/(Ru-R)/(Ru-Rl);
		e_[i] = -(Ru-R)/(R-Rl)/(Ru-Rl);
		f_[i] = (Ru-2*R+Rl)/(R-Rl)/(Ru-R);
		g_[i] = (R-Rl)/(Ru-Rl)/(Ru-R);	
	}
	for(int r=0; r<MR; r++) 
	{
		if(gridR_[r] >= R0_) {iR0_=r; break;}
	}
};

void HWPDE::discountBack(double t0, double Tn, vector<double> &f)
{
	int Nt = int(Tn/delT+1); 
	int N0 = int(t0/delT+1);
	for(int t=Nt-2; t>=N0; t--) oneStepBackward(t, f); 
};

double HWPDE::pricingZB(double Maturity)
{
	vector<double> payoff(MR,1.0);	
	discountBack(0, Maturity,payoff);
	return	payoff[iR0_];
};


void HWPDE::pricingCouponBondt(double t0, double Tn, double Coupon, double PayFrequency,vector<double> &f)
{
	int Nt = int(Tn/delT+1); //working daily time grid
	int iExpiry = int(t0/delT+1);
	int Ncoupon = int((Tn-t0)/PayFrequency+0.5);
	for(int t=Nt-2; t>=iExpiry; t--) //swap
	{
		double couponTime = t0 + Ncoupon *PayFrequency;
		if(couponTime > t*delT && couponTime <= (t+1)*delT)
		{
			for(int r=0; r<MR; r++) 
				f[r] += PayFrequency*Coupon;
			Ncoupon--;
		}
		oneStepBackward(t, f); //from t+1 to t
	}
};

double HWPDE::pricingCouponBond(double Expiry, double Tenor, double Coupon, double PayFrequency)
{
	vector<double> payoff(MR,1);
	pricingCouponBondt(Expiry,Expiry+Tenor,Coupon,PayFrequency,payoff);
	if (Expiry>0) discountBack(0,Expiry,payoff);
	return	payoff[iR0_];
};

double HWPDE::pricingCBO(double Expiry, double Tenor, double Coupon, double Strike, double PayFrequency, const std::string &type)
{
	vector<double> f(MR,1);
	int Nt =  int(Expiry/delT+1);
	pricingCouponBondt(Expiry,Expiry+Tenor,Coupon,PayFrequency,f); // Price CB
	for(int i = 0; i < MR ; i++) f[i] = (type=="Call")?max(0.0,f[i] - Strike):max(0.0,Strike - f[i]); //Apply Payoff
	if (Expiry>0) discountBack(0,Expiry,f);//Discount to time t=0
	return	f[iR0_]; //Find the correct value
	return -1; // If you ever hit here, you are doing something completely wrong, believe me man.
};

double HWPDE::pricingZBO(double Expiry, double Maturity, double Strike, const std::string &type)
{
	int Nt = int(Expiry/delT+1); //working daily time grid?
	vector<double> f(MR,1);
	//Price ZB from T to t
	discountBack(Expiry,Expiry+Maturity,f);
	//Apply Payoff to f: max(0,f-K) || max(0,K-f)
	for(int i =0; i < MR; i++) f[i] = (type=="Call")?max(0.0,f[i]-Strike):max(0.0,-f[i]+Strike);
	if (Expiry>0) discountBack(0,Expiry,f);//Discount to time t=0
	return	f[iR0_];
};


double HWPDE::pricingSwap(double Expiry, double Tenor, double Strike, double PayFrequency)
{
	vector<double> payoff(MR,-1);
	pricingCouponBondt(Expiry,Expiry+Tenor, -Strike, PayFrequency, payoff);
#pragma omp parallel for
	for(int r=0; r<MR; r++) payoff[r] += 1.0; //swap payoff
	if (Expiry>0) discountBack(0,Expiry,payoff);
	return	payoff[iR0_];
};

double HWPDE::pricingCallableSwap(double Expiry, double Tenor, vector<double> Exercises, double Coupon, double Strike, double PayFrequency, const std::string &type)
{
	vector<double> payoff(MR,-1),call_value(MR);
	Exercises.push_back(Expiry+Tenor);
	int Ne = Exercises.size();	
	int itype = (type=="Call")?1:-1; 
	for(int e=Ne-2; e>=0; e--)
	{
		pricingCouponBondt(Exercises[e],Exercises[e+1], -Coupon, PayFrequency, payoff);
		if (e==Ne-2)
		{
#pragma omp parallel for
			for(int r=0; r<MR; r++)	call_value[r] =max(0.0, ((payoff[r]+1.0)-Strike)*itype);// payoff	
		}
		else
		{
			discountBack(Exercises[e],Exercises[e+1],call_value);
#pragma omp parallel for
			for(int r=0; r<MR; r++)	call_value[r] =max(call_value[r], ((payoff[r]+1.0)-Strike)*itype);//payoff	
		}
	}
	if (Expiry>0) discountBack(0,Exercises[0],call_value);
	return	-itype*call_value[iR0_];
};


double HWPDE::pricingSwaption(double Expiry, double Tenor, double Strike, double PayFrequency)
{
	vector<double> payoff(MR,-1.0);
	pricingCouponBondt(Expiry,Expiry+Tenor, -Strike, PayFrequency, payoff);
	for(int r=0; r<MR; r++)	payoff[r] = max(0.0, payoff[r]+1.0); //swaption payoff
	if (Expiry>0) discountBack(0,Expiry,payoff);
	return	payoff[iR0_];
};


double HWPDE::pricingBermudan(double Expiry, double Tenor, vector<double> Exercises, double Strike, double PayFrequency)
{	
	std::vector<double> payoff(MR,-1.0), swaption(MR);
	Exercises.insert(Exercises.begin(),Expiry);
	Exercises.push_back(Expiry+Tenor);
	int Ne = Exercises.size();	
	for(int e=Ne-2; e>=0; e--)
	{
		pricingCouponBondt(Exercises[e],Exercises[e+1], -Strike, PayFrequency, payoff);
		if (e==Ne-2)
		{
#pragma omp parallel for
			for(int r=0; r<MR; r++)	swaption[r] = max(0.0, payoff[r]+1.0);//swaption payoff	
		}
		else
		{
			discountBack(Exercises[e],Exercises[e+1],swaption);
#pragma omp parallel for
			for(int r=0; r<MR; r++)	swaption[r] = max(swaption[r], payoff[r]+1.0);//bermudan swaption payoff	
		}
	}
	if (Expiry>0) discountBack(0,Expiry,swaption);
	return	swaption[iR0_];
};


vector<double> HWPDE::getDFs(vector<double>& timePoints)
{
	int N = timePoints.size();
	vector<int> iTs(N);
	vector<double> inV(MR,0.0), DFs(N,0.0);
	for(int i=0; i<N; i++)	iTs[i] = int(timePoints[i]/delT+1);

	inV[iR0_] = 2.0/(gridR_[iR0_+1]-gridR_[iR0_-1]);

	for(int i=0; i<N; i++)
	{
		int sT = 0;
		if(i > 0)	sT = iTs[i-1]-1;
		for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);
		DFs[i] = trapezoidal(inV);
		QL_ENSURE(DFs[i] <= 1.3, "pricingDF Functor Fails " << DFs[i]);
	}
	return DFs;
};

double HWPDE::getAnnuity(double Expiry, double Tenor, double PayFrequency) 
{
	double SwapLow = pricingSwap(Expiry,Tenor, 0.01, PayFrequency);
	double SwapHigh = pricingSwap(Expiry,Tenor, 0.02, PayFrequency);
	double Level = (SwapLow-SwapHigh)/0.01;
	QL_ENSURE(Level >= 0.0, "getAnnuity Functor Fails " << Level);
	return Level;
};
double HWPDE::getSwapRate(double Expiry, double Tenor, double PayFrequency) 
{
	double SwapLow = pricingSwap(Expiry,Tenor, 0.01, PayFrequency);
	double SwapHigh = pricingSwap(Expiry,Tenor, 0.02, PayFrequency);
	double SwapRate = 0.01 + SwapLow/(SwapLow-SwapHigh)*0.01;
//	QL_ENSURE(SwapRate >= -0.0005, "getSwapRate Functor Fails " << SwapRate); //Negative 50BP Rate
	return SwapRate;
};
double HWPDE::swaptionValueATM(double Expiry, double Tenor, double PayFrequency, double SwapRate, double VolATM)
{
	double Level = getAnnuity(Expiry,Tenor,PayFrequency);
	return Level*SwapRate*(2.0*cdf_normal(0.5*VolATM*sqrt(Expiry))-1.0);
};

double HWPDE::getImpVolATM(double Expiry, double Tenor, double PayFrequency)
{
	double SwapLow = pricingSwap(Expiry,Tenor, 0.00, PayFrequency);
	double SwapHigh = pricingSwap(Expiry,Tenor, 0.01, PayFrequency);
	double Level = (SwapLow-SwapHigh)/0.01;
	double SwapRate = SwapLow/Level;
	double swaptionATM = pricingSwaption(Expiry, Tenor, SwapRate, PayFrequency);
	double Lo = 0.001;	// 0.1 PP
	double Hi = 4.999;	// 500 PP
	double Vol = 0.5*(Lo+Hi);
	double swaptionVal = Level*SwapRate*(2.0*cdf_normal(0.5*Vol*sqrt(Expiry))-1.0);
	int Niter = 0;
	do {
		Niter++;      
		if (swaptionVal < swaptionATM) Lo = Vol;
		if (swaptionVal > swaptionATM) Hi = Vol;
		Vol = 0.5*(Lo+Hi);
		swaptionVal = Level*SwapRate*(2.0*cdf_normal(0.5*Vol*sqrt(Expiry))-1.0);
	} while ( fabs(1.0-swaptionVal/swaptionATM) >= 1.0E-6 && Niter < 40 );
	return Vol;
};


void HWPDE::calibrator(vector<double> timeDFs, vector<double> DFs, vector<defSwap> swapQuotes)
{
	int newThetas = timeDFs.size();
	timeDFs_ = timeDFs;
	timeThetas_.resize(1);
	timeThetas_[0] = .5;
	vector<double> aux;
	for (int i = 0; i < swapQuotes.size(); i++)
	{
		aux.push_back(swapQuotes[i].Expiry);
		aux.push_back(swapQuotes[i].Tenor);
		aux.push_back(swapQuotes[i].Expiry+swapQuotes[i].Tenor);
	}
	sort(aux.begin(), aux.end());

	timeDFs_ = timeDFs;
	timeThetas_.resize(1);
	timeThetas_[0] = aux[0];
	int i = 1;
	int k =1;
	do{
		if (timeThetas_[i-1]<aux[k])
		{
			timeThetas_.push_back(aux[k]);
			i++;
		}
		k++;
	}while((timeThetas_[i-1]<=swapQuotes[swapQuotes.size()-1].Expiry+swapQuotes[swapQuotes.size()-1].Tenor)&&(aux.size()>k));


	thetas_.resize(timeThetas_.size());
	for (int i = 0; i < timeThetas_.size();i++) thetas_[i] = 0.005;
	DFs_ = DFs;
	quoteSwap_ = swapQuotes;

	double pos;
	timeSigmas_.resize(0);
	timeSigmas_.push_back( quoteSwap_[0].Expiry );
	pos=timeSigmas_[0];
	for (int i=1; i < quoteSwap_.size();i++)
	{
		if (quoteSwap_[i].Expiry>pos)
		{
			timeSigmas_.push_back( quoteSwap_[i].Expiry );
			pos = quoteSwap_[i].Expiry;
		}
	}
	sigmas_=sigmas0_;
	sigmas_.resize(timeSigmas_.size());



	int n = sigmas_.size()+1;	    //no. of HW model volatility term paremeters 
	double* x = new double[n];		//initial estimate of parameters vector
	double* lb = new double[n];
	double* ub = new double[n];
	for(int i=0; i<n-1; i++)
	{
		x[i] = sigmas_[i];		// sigmas initial value
		lb[i] = 0;
		ub[i] = 1;
	}
	lb[n-1] = 0;
	ub[n-1] = 1.0;
	x[n-1] = kappa0_;			// kappa initial value

	int m = quoteSwap_.size();		//no. of observations
	int N = timeThetas_.size();
	marketSwaption_.resize(m);
	for (int i = 0; i<m;i++) marketSwaption_[i] = swaptionATM(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].VolATM);
	QL_ENSURE(m>=n, "too much freedom in Calibration  " << m-n);
	double* fvec = new double[m];	//no need to populate 
	double ftol = 1e-10; //tolerance
	double xtol = 1e-10; //tolerance
	double gtol = 1e-10; //tolerance
	int maxfev = 10000;  //maximum function evaluations
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
	boost::function<void (HWPDE*, int,int,double*,double*,int*, double *, double*)> obj = &HWPDE::objFcnCalibration;
	lmfcn fcn = boost::bind(obj, this, _1,_2,_3,_4,_5,lb,ub);
	cal_time_=omp_get_wtime();
	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcn);
	cal_time_=omp_get_wtime()-cal_time_;
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
	QL_ENSURE(info != 4, "Hull-White Model Calibration Fails " << info);
	//the below is output result
	for(int i=0; i<n-1; i++)
		sigmas_[i] = x[i];	// sigmas final value
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
	delete[] lb;
	delete[] ub;
};

double HWPDE::pen_fun(double*x, double*lb, double*ub, int n)
{
	double penalty = 0.0;
	double lambda = 1e8;
	for(int i=0; i<n; i++)penalty+=((x[i]<=lb[i])||(x[i]>=ub[i]))?lambda*pow(x[i],2):0;
	return penalty;
}

void HWPDE::objFcnCalibration(int m, int n, double* x, double* fvec, int* iflag, double *lb, double *ub)
{
	double penalty = pen_fun(x, lb, ub, n);
	for(int i=0; i<n-1; i++) sigmas_[i] = x[i];		
	kappa_ = x[n-1];				// kappa
	termStructureCalibrator();	// TEREM STRUCTURE CALIBRATION
	int i;
	for(i=0; i<m; i++)
	{
		double modelSwaption = pricingSwaption(quoteSwap_[i].Expiry, quoteSwap_[i].Tenor, quoteSwap_[i].SwapRate, quoteSwap_[i].Frequency);
		fvec[i] = modelSwaption - marketSwaption_[i]+penalty;
	}
};


void HWPDE::termStructureCalibrator()
{
	int N = timeThetas_.size();
	std::vector<int> iTs(N);
	for(int i=0; i<N; i++)
		iTs[i] = int(timeThetas_[i]/delT+1);
	std::vector<double> DFinterp(N);
	for (int i = 0 ; i < N; i++ ) DFinterp[i] = getDFinterp(timeThetas_[i]);
	//buildGrid(timeThetas_[N-1], iTs[N-1]);
	std::vector<double> inV(MR,0.0), lastV(MR,0);
	for(int r=0; r<MR; r++)
		if(gridR_[r] >= R0_)
		{ // Dirac delta as initial distribution 
			inV[r] = 2.0/(gridR_[r+1]-gridR_[r-1]);
			break;
		}
#pragma omp parallel for
		for (int i = 0; i < MR ; i++)
			lastV[i] = inV[i]; // remember the starting density
		for(int i=0; i<N; i++)
		{
			int sT = 0;
			if(i > 0) 
				sT = iTs[i-1] - 1;
			double df = 0.0;
			int niter = 0;
			do {
				niter++;
#pragma omp parallel for
				for (int i = 0; i < MR ; i++) inV[i] = lastV[i]; // reset to the starting density

				if(thetas_[i] == 0.0) thetas_[i] = 0.0005;	// case of 0 value handler
				thetas_[i] *= 1.001;  // 0.1% up for theta parameter
				for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);
				double dfUP = trapezoidal(inV);

#pragma omp parallel for
				for (int i = 0; i < MR ; i++)  inV[i] = lastV[i]; // reset to the starting density

				thetas_[i] /= 1.001;  // back to theta 
				for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);
				df = trapezoidal(inV);
				thetas_[i] *= (1 + 0.001*(DFinterp[i]-df)/(dfUP-df));	// Newton iteration
				thetas_[i] = max(-0.0050, thetas_[i]);	// Cutoff the nagative value (BETTER)
			} while( fabs(1.0-df/DFinterp[i]) >= 1.0E-6 && niter < 10 );
#pragma omp parallel for
			for (int i = 0; i < MR ; i++) inV[i] = lastV[i]; // reset to the starting density
			for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);
#pragma omp parallel for
			for (int i = 0; i < MR ; i++)
				lastV[i] = inV[i]; // remember the starting density
		};
		return;
};

double HWPDE::trapezoidal(std::vector<double>& inV)
{
	double value = 0.0;
#pragma omp parallel for reduction(+:value)
	for(int r=1; r<MR; r++) 
		value += 0.5*(gridR_[r]-gridR_[r-1])*(inV[r]+inV[r-1]);
	return value;
};
double HWPDE::whichSigma(double t) const
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
double HWPDE::whichTheta(double t) const
{
	int n = thetas_.size();
	if(t < timeThetas_[0])
		return thetas_[0];
	if(t >= timeThetas_[n-1])
		return thetas_[n-1];
	for(int i=1; i<n; i++)
		if((t >= timeThetas_[i-1]) && (t < timeThetas_[i]))
			return thetas_[i];
};

const double DT = 0.001;
const double SDT = sqrt(DT);
vector<double> HWPDE::simulationHWPDE(vector<double> times) const
{
	int N=times.size();
	vector<double> path(N);
	int T = 0;
	double r = R0_;
	for(int i=0; i<=int(times[N-1]/DT+0.5); i++)
	{
		if(i*DT <= times[T] && (i+1)*DT > times[T])
		{
			path[T] = r;
			T++;
		}
		double sigma = whichSigma(i*DT);
		double theta = whichTheta(i*DT);
		r += (theta-kappa_*r)*DT + sigma*SDT*random_normal();  
	}
	return path;
};


void HWPDE::oneStepBackward(const int t, vector<double> &inV)
{
	//std::vector<double> l(MR-2), c(MR-2), u(MR-2), d(MR-2), V(MR-2);
	array<double, MR> l,c,u,d;
	double tm = (t+.5)*delT;
	double sigma = whichSigma(tm);
	double f2 = sigma*sigma/4;
	double theta = whichTheta(tm);
	double ndt = 1.0/delT;

#pragma omp parallel for 
	for(int r=0; r<MR; r++)
	{
		double f1 =((r==0)||(r==MR-1))? 0 : 0.5*(theta-kappa_*gridR_[r]);
		l[r] = f1*e_[r]+f2*a_[r];
		c[r] = -ndt+f1*f_[r]+f2*b_[r]-gridR_[r]/2;
		u[r] = f1*g_[r]+f2*c_[r];
		if (r==0) 
		{
			u[r] +=l[r];
			d[r] = -((2*ndt+c[r])*inV[r]+u[r]*inV[r+1]);
		}
		else if (r==MR-1) 
		{
			l[r] += u[r];
			d[r] = -(l[r]*inV[r-1]+(2.0*ndt+c[r])*inV[r]);
		}
		else
			d[r] = -(l[r]*inV[r-1]+(2.0*ndt+c[r])*inV[r]+u[r]*inV[r+1]);
	}
	l[0] = 0.0;
	u[MR-1] = 0.0;
	TriDiagonalSolve(MR,l,c,u,d,inV);	
};

void HWPDE::oneStepForward(const int T, vector<double> &inV)
{
	//std::vector<double> l(MR-2), c(MR-2), u(MR-2), d(MR-2), V(MR-2);
	array<double, MR> l,c,u,d;
	double Tm = (T+.5)*delT;
	double sigma = whichSigma(Tm);
	double f2 = -sigma*sigma/4;
	double theta = whichTheta(Tm);
	double ndt = 1.0/delT;

#pragma omp parallel for 



	for(int r=0; r<MR; r++)
	{
		double f1 =((r==0)||(r==MR-1))? 0 : 0.5*(theta-kappa_*gridR_[r]);
		double f3 = (gridR_[r]-kappa_)/2;
		l[r] = f1*e_[r]+f2*a_[r];
		c[r] = ndt+f1*f_[r]+f2*b_[r]+f3;
		u[r] = f1*g_[r]+f2*c_[r];
		if (r==0) 
		{
			u[r] +=l[r];
			d[r] = -((-2*ndt+c[r])*inV[r]+u[r]*inV[r+1]);
		}
		else if (r==MR-1) 
		{
			l[r] += u[r];
			d[r] = -(l[r]*inV[r-1]+(-2.0*ndt+c[r])*inV[r]);
		}
		else
			d[r] = -(l[r]*inV[r-1]+(-2.0*ndt+c[r])*inV[r]+u[r]*inV[r+1]);
	}
	l[0] = 0.0;
	u[MR-1] = 0.0;
	TriDiagonalSolve(MR,l,c,u,d,inV);	
};

vector<double> HWPDE::grMesh (int Mv, double Rmax, double factor)
{
	int N = Mv;
	vector<double> Mesh(Mv+1);
	double d3 = Rmax/(N*factor);
	double start = asinh((-Rmax-R0_)/d3);
	double dx = 1.0/N*(asinh((Rmax - R0_)/d3)-start);

	for (int i = 0; i <= N ; i++) Mesh[i] = R0_ + d3 * sinh(start + i*dx);
	vleft_=Mesh[0]-(Mesh[1]-Mesh[0]);
	vright_=Mesh[Mv]+(Mesh[Mv]-Mesh[Mv-1]);
	for (int i = 0; i<N; i ++)
		if ((Mesh[i]<R0_)&&(R0_<Mesh[i+1])) {Mesh[i+1] = R0_;break;}
		return Mesh;
};


void HWPDE::calibrateKappa()
{
	int MAX_ITER = 1000;
	double X_TOL = 1e-5;
	double F_TOL  = 1e-5;
	double a0 = .1;
	double a1 = a0;
	double dF, F;
	double min = 1e-6;
	double max = 1;
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
				(getDFinterp(ex) - getDFinterp(t1)) / 
				(getDFinterp(ex) - getDFinterp(t2))
				);
		}		
	}

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
		if ( fabs(max - min) < X_TOL ) break;
	}
	QL_ENSURE(a0==a0, "Failed to calibrate Kappa \n");
	kappa_ = a0;
};

double HWPDE::Bratio(double a, double Mi, double Tk, double Tj)
{
	if (a == 0.0) return 100000;
	else 
	{ 
		double t1 = (Tj-Mi);
		double t2 = (Tk-Mi);
		return (1.0 - exp(-a*t1))/(1.0 - exp(-a*t2));
	}
};
double HWPDE::dBratio(double a, double Mi, double Tk, double Tj)
{
	if (a == 0.0) return 100000;
	else 
	{ 
		double t1 = (Tj-Mi);
		double t2 = (Tk-Mi);
		double dF = exp(a*(t2-t1))*(t1*(exp(t2*a)-1)-t2*exp(t1*a)+t2)/(exp(t2*a)-1.0)/(exp(t2*a)-1.0);
		return dF; 
	}
};

double HWPDE::objKappa(double a, const vector<double> &ivr, const vector<double> &br, const vector<int> &ind)
{
	int N = ivr.size();
	double sum = 0.0;
	for (int i = 1; i< N; i++)
	{
		double ex = quoteSwap_[ind[i]].Expiry;
		double t1 = quoteSwap_[ind[i]].Tenor + ex;
		double t2 = quoteSwap_[ind[i]+1].Tenor + ex;
		double vs = fabs(br[i]*Bratio(a, ex, t1, t2))-ivr[i];
		sum+=vs*vs;
	}
	return sum;
};

double HWPDE::getDFinterp(double t)
{
	QuantLib::LogLinearInterpolation interp(timeDFs_.begin(),timeDFs_.end(),DFs_.begin());
	double a = interp(t,false);
	return a;
};


void HWPDE::calibratorBootStrap(vector<double> timeDFs, vector<double> DFs, vector<defSwap> swapQuotes)
{
	int MAX_ITER = 50;
	double X_TOL = 1e-4;
	double F_TOL  = 1e-5;
	double dF;
	double F;
	double min = 1e-6;
	double max = 1;
	double x = sigmas0_[0];

	vector<double> aux;
	for (int i = 0; i < swapQuotes.size(); i++)
	{
		aux.push_back(swapQuotes[i].Expiry);
		aux.push_back(swapQuotes[i].Tenor);
		aux.push_back(swapQuotes[i].Expiry+swapQuotes[i].Tenor);
	}
	sort(aux.begin(), aux.end());

	timeDFs_ = timeDFs;
	timeThetas_.resize(1);
	timeThetas_[0] = aux[0];
	int i = 1;
	int k =1;
	do{
		if (timeThetas_[i-1]<aux[k])
		{
			timeThetas_.push_back(aux[k]);
			i++;
		}
		k++;
	}while((timeThetas_[i-1]<=swapQuotes[swapQuotes.size()-1].Expiry+swapQuotes[swapQuotes.size()-1].Tenor)&&(aux.size()>k));

	thetas_.resize(timeThetas_.size(),0.005);
	for (int i = 0; i < timeThetas_.size();i++) thetas_[i] = 0.005;
	DFs_ = DFs;

	quoteSwap_ = swapQuotes;

	calibrateKappa(); //Calibrate Kappa

	double pos;
	timeSigmas_.resize(0);
	timeSigmas_.push_back( quoteSwap_[0].Expiry );
	pos=timeSigmas_[0];
	for (int i=1; i < quoteSwap_.size();i++)
	{
		if (quoteSwap_[i].Expiry>pos)
		{
			timeSigmas_.push_back( quoteSwap_[i].Expiry );
			pos = quoteSwap_[i].Expiry;
		}
	}
	sigmas_.resize( timeSigmas_.size());
	k = 0;
	for(int j = 0; j < timeSigmas_.size();j++)
	{		
		iter_ = j;

		qSB_.resize(0);
		for (int l = k; l < quoteSwap_.size() ; l++,k++)
		{			
			if(l<=quoteSwap_.size())
			{
				if (!(quoteSwap_[l].Expiry <= timeSigmas_[j])) {break;}
				else qSB_.push_back(quoteSwap_[l]);
			}
		}

		int m = qSB_.size();		//no. of observations
		int n = 1;
		marketSwaption_.resize(m);
		for (int i = 0; i < m; i++) marketSwaption_[i] = swaptionValueATM(qSB_[i].Expiry, qSB_[i].Tenor, qSB_[i].Frequency, qSB_[i].SwapRate, qSB_[i].VolATM);


		double best, min_val=1000;

		double delta;
		cal_time_=omp_get_wtime();
		for (int i = 0; i < MAX_ITER ; i++)
		{
			F = objFcnCalibrationB(x,m);
			double Fu = objFcnCalibrationB(x+X_TOL,m);
			double Fd = objFcnCalibrationB(x-X_TOL,m);
			dF = 0.5*(Fu-Fd)/X_TOL;
			if (min_val>F){min_val = F; best = x;}
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
			if ( fabs(max - min) < X_TOL ) break;
		}
		sigmas_[j] = x;
		min = 1e-6;
		max = 1.0;
	}	
	cal_time_=omp_get_wtime()-cal_time_;
};

double HWPDE::objFcnCalibrationB(double x, int m)
{
	sigmas_[iter_] = x;	// sigmas
	termStructureCalibratorBtstrp(qSB_[m-1].Expiry+qSB_[m-1].Tenor);	// TEREM STRUCTURE CALIBRATION
	double sum = 0.0;
	for(int i=0; i<m; i++)
	{
		double modelSwaption = pricingSwaption(qSB_[i].Expiry, qSB_[i].Tenor, qSB_[i].SwapRate, qSB_[i].Frequency);
		double a = fabs(1.0-modelSwaption/marketSwaption_[i]);
		sum+=a;
	}
	return sum/m;
};

double HWPDE::swaptionIVblack(double Expiry, double Tenor, double swaption_price)
{
	double xl = 1e-6;
	double xu = 1;
	double x = 0.5*(xl+xu);

	double DF = getDFinterp(Expiry)-getDFinterp(Expiry+Tenor);
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

void HWPDE::termStructureCalibratorBtstrp(double max_maturity)
{
	int N = 0;
	if (timeThetas_[0]>max_maturity) N = 1;
	else
		for (int i = 1; i< timeThetas_.size();i++)
			if((max_maturity>timeThetas_[i-1])&&(max_maturity<=timeThetas_[i]))
			{N = i; break;}
			std::vector<int> iTs(N);
			for(int i=0; i<N; i++)
				iTs[i] = int(timeThetas_[i]/delT+1);
			//buildGrid(timeThetas_[N-1], iTs[N-1]);
			std::vector<double> inV(MR,0.0), lastV(MR);
			for(int r=0; r<MR; r++)
				if(gridR_[r] >= R0_)
				{ // Dirac delta as initial distribution 
					inV[r] = 2.0/(gridR_[r+1]-gridR_[r-1]);
					break;
				}
				lastV = inV; // remember the starting density
				for(int i=0; i<N; i++)
				{
					int sT = 0;
					if(i > 0) 
						sT = iTs[i-1] - 1;
					double df = 0.0;
					int niter = 0;
					do {
						niter++;
						inV = lastV; // reset to the starting density
						if(thetas_[i] == 0.0) thetas_[i] = 0.0005;	// case of 0 value handler
						thetas_[i] *= 1.001;  // 0.1% up for theta parameter
						for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);
						double dfUP = trapezoidal(inV);
						inV = lastV; // reset to the starting density
						thetas_[i] /= 1.001;  // back to theta 
						for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);

						df = trapezoidal(inV);
						thetas_[i] *= (1 + 0.001*(getDFinterp(timeThetas_[i])-df)/(dfUP-df));	// Newton iteration
						thetas_[i] = max(-0.0050, thetas_[i]);	// Cutoff the nagative value (BETTER)
					} while( fabs(1.0-df/getDFinterp(timeThetas_[i])) >= 1.0E-6 && niter < 7 );
					inV = lastV; // reset to the starting density
					for(int t=sT; t<iTs[i]-1; t++)	oneStepForward(t, inV);
					lastV = inV; // remember the starting density
				};
				return;
};


}