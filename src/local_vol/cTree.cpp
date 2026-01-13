#include <velesquant/local_vol/cTree.h>
using namespace std;
#include <algorithm>

namespace velesquant {



cTree::cTree(double S, Vdoub &T, Vdoub &F, Vdoub &IV)
	{
	spot_=S;
	T_=T;
	F_=F;
	IV_=IV;
	};

double cTree::calculateBinomial(double strike, double Maturity, int Nnodes, exStyle style, pType pay, tType tree)
	{
	double	dt	= Maturity / double(Nnodes);
	if (r_.size()==0) r_.push_back(.01); //To fix
	if (q_.size()==0) q_.push_back(0.002); //To fix
	Mdoub	myTree;
	Vdoub p(Nnodes+1, 0);
	Vdoub	up(Nnodes), pu(Nnodes);

	for (int i = 0 ; i < Nnodes ; i++)
		{
		up[i] =  exp(interp_vol(T_,IV_,i*dt, (i+1)*dt)*sqrt(dt));
		double d = 1.0/up[i];
		pu[i] = (exp((r_[0]-q_[0])*dt)-d)/(up[i]-d);
		}

	double type = (pay==Call)?-1.0:1.0;
	for (int i = 0 ; i <= Nnodes; i++)
		{
		p[i] = max((strike - spot_ * pow(up[Nnodes-1], 2 * i - Nnodes))*type,0.0);
		}
	for (int j = Nnodes - 1; j>=0 ; j--)
		{
		double rdt = exp(-r_[0]*dt);
		for (int i = 0; i <= j ; i++)
			{
			if (style==European)
				p[i] =	rdt	*	(	(1.0 - pu[j]) * p[i] + pu[j] * p[i+1]) ;
			else
				{
				double aux = max( (strike - spot_ * pow(up[j], double (2 * i - j) ) )*type, 0.0);
				p[i] =	max( aux,	rdt	*	(	(1.0 - pu[j]) * p[i] + pu[j] * p[i+1]) );

				}

			}
		}
	return p[0];
	};

double cTree::calculateBinomial(double strike, double Maturity, int Nnodes, exStyle style, pType pay, tType tree, Vdoub schedule)
	{
	double	dt	= Maturity / double(Nnodes);
	if (r_.size()==0) r_.push_back(.01); //To fix
	if (q_.size()==0) q_.push_back(0.002); //To fix
	Mdoub	myTree;
	Vdoub p(Nnodes+1, 0);
	Vdoub	up(Nnodes), pu(Nnodes);

	for (int i = 0 ; i < Nnodes ; i++)
		{
		up[i] =  exp(interp_vol(T_,IV_,i*dt, (i+1)*dt)*sqrt(dt));
		double d = 1.0/up[i];
		pu[i] = (exp((r_[0]-q_[0])*dt)-d)/(up[i]-d);
		}

	double type = (pay==Call)?-1.0:1.0;
	for (int i = 0 ; i <= Nnodes; i++)
		{
		p[i] = max((strike - spot_ * pow(up[Nnodes-1], 2 * i - Nnodes))*type,0.0);
		}
	for (int j = Nnodes - 1; j>=0 ; j--)
		{
		double rdt = exp(-r_[0]*dt);
		for (int i = 0; i <= j ; i++)
			{
			if (style==European)
				p[i] =	rdt	*	(	(1.0 - pu[j]) * p[i] + pu[j] * p[i+1]) ;
			else
				{
				double aux = max( (strike - spot_ * pow(up[j], double (2 * i - j) ) )*type, 0.0);
				p[i] =	max( aux,	rdt	*	(	(1.0 - pu[j]) * p[i] + pu[j] * p[i+1]) );

				}

			}
		}
	return p[0];
	};

double cTree::interp_vol(vector<double> &T, vector<double> &F, double t1, double t2) const
	{
	double sig0=interpolate(T,F,t1);
	double sigT=interpolate(T,F,t2);
	double dt = t2-t1;
	if (dt<=0.0) throw("Incorrect dt");
	if (sig0==sigT)
		return sig0;
	else 
		return sqrt((t2*sigT*sigT-t1*sig0*sig0)/dt);
	};

double cTree::interpolate(vector<double> &T, vector<double> &F, double t) const
	{
	int i,N = T.size();
	for (i = 0; i<N ; i++) if (t<T[i]) break;
	if (i==0) 
		return F[i]/T[i]*t;
	if ((i>0)&&(i<N))
		return F[i-1]+(F[i]-F[i-1])/(T[i]-T[i-1])*(t-T[i-1]);
	else 
		return F[N-2]+(t-T[N-2])/(T[N-1]-T[N-2])*(F[N-1]-F[N-2]);

	};

double cTree::calculateTrinomial(double strike, double Maturity, int Nnodes, exStyle style, pType pay, tType tree)
	{
	double	dt	= Maturity / double(Nnodes);
	if (r_.size()==0) r_.push_back(.01); //To fix
	if (q_.size()==0) q_.push_back(0.02); //To fix
	int p_size = 2*Nnodes+1;
	Vdoub	p(p_size, 0);
	Vdoub	up(Nnodes), pu(Nnodes), pd(Nnodes);
	double pm;
	double C;
	double dt2 = sqrt(0.5*dt);
	for (int i = 0 ; i < Nnodes ; i++)
		{
		double  sig = interp_vol(T_,IV_,i*dt, (i+1)*dt);
		double  s2 =  sig * dt2;
		double  rt =  (r_[0]-q_[0])*dt*0.5 ;
		double dmt =  exp(s2) -  exp(-s2) ;
		up[i] = exp( sig * sqrt(2.0*dt));
		pu[i] = ( exp( rt ) - exp( -s2	) ) / dmt;
		pu[i]*= pu[i];
		pd[i] = ( exp( s2 ) - exp( rt	) ) / dmt;
		pd[i]*= pd[i];
		}

	double type = (pay==Call)?-1.0:1.0;
	for (int i = 0 ; i < p_size; i++)
		{
		double payoff = (strike - spot_ * pow(up[Nnodes-1], i - Nnodes))*type;
		p[i] = max(payoff,0.0);
		}

	for (int j = Nnodes - 1; j>=0 ; j--)
		{
		double rdt = exp(-r_[0]*dt);
		for (int i = 0; i <= 2*j ; i++)
			{
			pm = (1.0 - pd[j]-pu[j]);
			C = rdt*(pd[j]*p[i]+pm*p[i+1]+pu[j]*p[i+2]);
			if (style == American)
				{
				double aux = max((strike - spot_ * pow(up[j], double (i - j) ))*type,0.0);
				p[i] = max(aux, C);
				}
			else 
				p[i] = C;
			}
		}
	return p[0];
	};

}