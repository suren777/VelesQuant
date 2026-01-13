#include <velesquant/local_vol/schobzhu.h>
#include <complex>
#include <velesquant/models/utility.h>
#include <velesquant/local_vol/lm.h>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <ql/quantlib.hpp>
#include <algorithm>

using namespace std;
using namespace QuantLib;
using namespace boost::placeholders;

namespace velesquant {



const double PI = 3.14159265358979323846264338327950288;
const double DT = 0.001;
const double SDT = sqrt(DT);

SchobZhu::SchobZhu(double spot, double sigma0, double kappa, double theta,  double xi,  double rho):
s0_(spot),  sigma0_(sigma0), kappa_(kappa), theta_(theta),  xi_(xi), rho_(rho){};

double SchobZhu::SchobelIntegrand(double u, double maturity, double forward, double strike, int type) const 
	{ 
	Cdoub li(0,1);
	Cdoub X = (type == 0)?li*u:(1.0+li*u);


	double xi2 = xi_*xi_;

	Cdoub a = X*(1-2*kappa_*rho_/xi_-X*(1-rho_*rho_));
	Cdoub b = X*kappa_*theta_*rho_/xi_;
	Cdoub c = X*rho_/xi_;
	Cdoub alpha = std::sqrt(xi2*a+kappa_*kappa_);
	Cdoub beta = (kappa_-xi2*c)/alpha;
	Cdoub gamma = kappa_*kappa_*theta_-xi2*b;

	Cdoub saT = std::sinh(alpha*maturity);
	Cdoub caT = std::cosh(alpha*maturity);
	Cdoub D = 1.0/xi2*(kappa_-alpha*(saT+beta*caT)/(caT+beta*saT));
	Cdoub B = 1.0/xi2/alpha*((kappa_*theta_*alpha-beta*gamma+gamma*(saT+beta*caT))/(caT+beta*saT)-kappa_*theta_*alpha);
	Cdoub C = .5*(kappa_*maturity-log(caT+beta*saT))
		+(.5*(kappa_*theta_*alpha*kappa_*theta_*alpha-gamma*gamma) * (saT/(caT+beta*saT)-alpha*maturity)
		+(kappa_*theta_*alpha-beta*gamma)*gamma*(caT-1.0)/(caT+beta*saT))/(xi2*alpha*alpha*alpha);
	Cdoub A = X*(sigma0_*sigma0_+xi2*maturity)*rho_/xi_;

	Cdoub Psi = exp(-.5*A+.5*D*sigma0_*sigma0_+B*sigma0_+C);

	return std::real(exp(li*u*log(forward/strike))*Psi/(li*u));

	}

double SchobZhu::SchobelPrice(double maturity, double forward, double strike) const
	{
	boost::function<double (double)> Phi0, Phi1;
	Phi0 = boost::bind(&SchobZhu::SchobelIntegrand, this, _1, maturity, forward, strike, 0);
	Phi1 = boost::bind(&SchobZhu::SchobelIntegrand, this, _1, maturity, forward, strike, 1);
	GaussLaguerreIntegration gLegInt(16);
	double callPrice=s0_*(0.5+(1/PI)*gLegInt(Phi1))-strike*(0.5+(1/PI)*gLegInt(Phi0));
	return callPrice;

	}


Vdoub SchobZhu::simulationSchobZhu(Vdoub times, Vdoub forwards) const
	{
	int N=times.size();
	vector<double> pathF(N);
	int T = 0;
	double V=sigma0_, S=s0_, rho2 = sqrt(1-rho_*rho_);
	for(int i=0; i<=int(times[N-1]/DT); i++)
		{
		double z1 = random_normal();
		double z2 = rho_*z1 + rho2*random_normal();
		S += S*V*SDT*z1;  // This line (S) must put BEFORE line (V)!!!
		V += kappa_*(theta_-V)*DT +xi_*SDT*z2;  //Euler scheme
		if(i*DT <= times[T] && (i+1)*DT > times[T])
			{
			pathF[T] = forwards[T]*S/s0_;
			T++;
			}
		}
	return pathF;
	
	
	}

void SchobZhu::objFcn(int m, int n, double* x, double* fvec, int* iflag)
	{
	setParameterVar0(x[0]);
	setParameterKappa(x[1]);
	setParameterTheta(x[2]);
	setParameterXi(x[3]);
	setParameterRho(x[4]);
	for(int i=0; i<m; i++)
		fvec[i] = SchobelPrice(maturitys_[i], forwards_[i], strikes_[i]) - marketQuotes_[i];
	};



void SchobZhu::calibrator(Vdoub maturitys, Vdoub forwards, Vdoub strikes, Vdoub marketQuotes) 
	{
	int m = strikes.size();		//no. of observations
	maturitys_.resize(m);		
	maturitys_ = maturitys; 
	forwards_.resize(m);			
	forwards_ = forwards;
	strikes_.resize(m);			
	strikes_ = strikes;
	marketQuotes_.resize(m);
	marketQuotes_ = marketQuotes;

	int n = 5;					//no. of SZ model paremeters 
	double* x = new double[n];	//initial estimate of parameters vector
	x[0] = getParameterVar0();
	x[1] = getParameterKappa();
	x[2] = getParameterTheta();
	x[3] = getParameterXi();
	x[4] = getParameterRho();

	double* fvec = new double[m]; //no need to populate 
	double ftol = 1e-08; //tolerance
	double xtol = 1e-08; //tolerance
	double gtol = 1e-08; //tolerance
	int maxfev = 400; //maximum function evaluations
	double epsfcn = 1e-08; //tolerance
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

	boost::function<void (SchobZhu*, int,int,double*,double*,int*)> objSZ = &SchobZhu::objFcn;
	lmfcn fcnSZ = boost::bind(objSZ, this, _1,_2,_3,_4,_5);
	lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
		nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, fcnSZ);
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
	QL_ENSURE(info != 4, "Schob Zhu Model Calibration Fails " << info);
	//the below is output result
	setParameterVar0(x[0]);
	setParameterKappa(x[1]);
	setParameterTheta(x[2]);
	setParameterXi(x[3]);
	setParameterRho(x[4]);

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
}