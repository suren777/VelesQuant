#include <velesquant/models/logBasket.h>
#include <velesquant/models/utility.h>
#include <algorithm>
#include <cmath>
#include <ql/quantlib.hpp>
using namespace std;

namespace velesquant {



lBasket::lBasket(Vdoub Spot, 
	Vdoub Strike, 
	Vdoub Maturities, 
	Mdoub Forwards, 
	Mdoub IV, 
	Mdoub correlation)
	{
	maturity_ = Maturities;
	spot_=Spot; 
	strike_=Strike;
	forward_=Forwards; 
	iv_=IV;
	correlation_=cholesky(correlation);
	Nassets_=Spot.size();
	};

lBasket::~lBasket()
	{
	maturity_.clear();
	spot_.clear(); 
	strike_.clear();
	forward_.clear(); 
	iv_.clear();
	correlation_.clear();
	};
void lBasket::update_seed(int seed)
	{
	set_seed(seed);
	}

double lBasket::generate_spot(double sigma, double dT, double Noise)
	{
	double sigT = sigma*sqrt(dT);
	return (Noise-0.5*sigT)*sigT;
	}

Mdoub lBasket::simulate_basket(Vdoub schedule)
	{
	int N = schedule.size();
	double dT,t1 = 0.0;
	Vdoub Noise(Nassets_);
	Vdoub S(Nassets_);
	Vdoub cNoise(Nassets_);
	Mdoub pathF(Nassets_,vector<double>(N));
	S = spot_;
	for (int k = 0; k<N;k++)
		{
		for (int i = 0; i<Nassets_ ; i++)
			Noise[i] = random_normal();
		lAxy(correlation_, Noise, cNoise);
		if (k==0) dT = schedule[k];
		else dT=(schedule[k]-schedule[k-1]);
		for (int i = 0; i < Nassets_; i++)
			{	
			double sig = interp_vol(maturity_, iv_[i],t1,t1+dT);
			S[i]=S[i]*exp(generate_spot(sig, dT, cNoise[i]));				
			pathF[i][k]=S[i]*interpolate(maturity_,forward_[i],t1+dT,i)/spot_[i];
			}
		t1+=dT;
		}

	return pathF;
	};

double lBasket::interp_vol(vector<double> &T, vector<double> &F, double t1, double t2) const
	{
	double sig0=interpolate(T,F,t1,-1);
	double sigT=interpolate(T,F,t2,-1);
	return sqrt((t2*sigT*sigT-t1*sig0*sig0)/(t2-t1));
	};

double lBasket::interpolate(vector<double> &T, vector<double> &F, double t, int asset) const
	{
	int i,N = T.size();
	double spot =0.0;
	if (asset == -1) spot = F[0];
	else spot = spot_[asset];
	if (t==0.0) return F[0];
	for (i = 0; i<N ; i++) if (t<T[i]) break;
	if (i==0) 
		return spot+(F[i]-spot)/(T[i])*(t);
	if ((i>0)&&(i<N))
		return F[i-1]+(F[i]-F[i-1])/(T[i]-T[i-1])*(t-T[i-1]);
	else 
		return F[N-2]+(t-T[N-2])/(T[N-1]-T[N-2])*(F[N-1]-F[N-2]);

	};

Vdoub lBasket::simulate_basketWR(Vdoub schedule)
	{
	int N = schedule.size();
	double dT,t1 = 0.0;
	Vdoub Noise(Nassets_);
	Vdoub S(Nassets_);
	Vdoub cNoise(Nassets_);
	Vdoub pathF(N);
	S = spot_;
	for (int k = 0; k<N;k++)
		{
		for (int i = 0; i<Nassets_ ; i++)
			Noise[i] = random_normal();
		lAxy(correlation_,Noise,cNoise);
		if (k==0) dT = schedule[k];
		else dT=(schedule[k]-schedule[k-1]);
		for (int i = 0; i < Nassets_; i++)
			{	
			double sig = interp_vol(maturity_, iv_[i],t1,t1+dT);
			S[i]=S[i]*exp(generate_spot(sig, dT, cNoise[i]));
			double aux = S[i]*interpolate(maturity_,forward_[i],t1+dT, i)/spot_[i];
			aux/=spot_[i];
			if (i==0) 
				pathF[k]=aux;
			else if (pathF[k]>aux) pathF[k] = aux;
			}
		t1+=dT;
		}		
	return pathF;
	};

double lBasket::sim_basket_with_removal(double tT, double size, Vint &ishares, const Vdoub &barriers, const Vdoub &coupons,const Vdoub &initial_price, Vdoub &S, int Nsteps, int called)
	{
	double dT=size/Nsteps;
	double t1 = tT;
	double coupon = 0.0;
	Vdoub Noise(Nassets_);
	Vdoub cNoise(Nassets_);
	Vdoub step(Nassets_);
	int takeout;
	double min,min2;
	min = min2 = 1000.0;
	for (int k = 0; k<Nsteps;k++)
		{
		for (int i = 0; i<Nassets_ ; i++)
			Noise[i] = random_normal();
		lAxy(correlation_, Noise, cNoise);
		for (int i = 0; i < Nassets_; i++)
			{			
			if (ishares[i] == 1)
				{
				double sig = interp_vol(maturity_, iv_[i],t1,t1+dT);
				S[i]*=exp(generate_spot(sig, dT, cNoise[i]));				
				step[i]=S[i]/initial_price[i]*interpolate(maturity_,forward_[i],t1+dT,i)/spot_[i];
				min = (min>step[i])?step[i]:min;
				}
			}
		t1+=dT;
		}
	for (int i=0; i < Nassets_ ; i++)
		{
		if ((min2>step[i])&&(ishares[i]==1))
			{
			min2 = step[i] ; 
			takeout = i;
			}
		}
	if (vsum(ishares) == 0) return 0.0;
	if (min<barriers[1]) coupon=coupons[2];
	else if (min<barriers[0]) coupon = coupons[1];
	else if (min>=barriers[0]) 
		{
		coupon = coupons[0];
		if ( called == 1 ) for (int i=0; i<int(ishares.size());ishares[i++]=0);
		}
	ishares[takeout]=0;	
	return coupon;
	};

double lBasket::vsum(Vdoub &vec)
	{
	double sum=0.0;
	for (int i = 0; i < int(vec.size());i++)	sum+=vec[i];
	return sum;
	};

int lBasket::vsum(Vint &vec)
	{
	int sum=0;
	for (int i = 0; i < int(vec.size());i++)	sum+=vec[i];
	return sum;
	};
void lBasket::get_spots(Vdoub &S)
	{
	for (int i = 0; i < int (spot_.size()) ;i++) S[i]=spot_[i] ;
	};
}