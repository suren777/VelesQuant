//	lognormal_basket.h

#ifndef LOGBASKET_H
#define LOGBASKET_H

#include <vector>
#include <complex>

namespace velesquant {

// using namespace std; // REMOVED by refactor

typedef std::vector<double>			Vdoub;
typedef std::vector<int>				Vint;
typedef std::vector<std::vector<double>>	Mdoub;
class lBasket
{
public:
	lBasket(Vdoub Spot, Vdoub Strike, Vdoub Maturities, Mdoub Forwards, Mdoub IV, Mdoub correlation);
	~lBasket();
	
	Mdoub simulate_basket(Vdoub schedule);
	Vdoub simulate_basketWR(Vdoub schedule);
	int get_nassets(){return Nassets_;};
	double generate_spot(double sigma, double dT, double Noise);
	void update_seed(int seed);
	void get_spots(Vdoub &S);
	double sim_basket_with_removal(double tT, double size, Vint &ishares, const Vdoub &barriers, const Vdoub &coupons,const Vdoub &initial_price, Vdoub &spot, int Nsteps, int called);
private:
	Vdoub maturity_, spot_, strike_;
	Mdoub forward_, iv_, correlation_;
	int Nassets_;
	double interp_vol(std::vector<double> &T, std::vector<double> &F, double t1, double t2) const;
	double interpolate(std::vector<double> &T, std::vector<double> &F, double t, int asset) const;
	double vsum(Vdoub &vec);
	int vsum(Vint &vec);
};

}
#endif