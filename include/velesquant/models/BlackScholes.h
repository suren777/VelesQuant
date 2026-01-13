
namespace velesquant {
// BlackScholes.h

#ifndef BLACKSCHOLES_H
#define BLACKSCHOLES_H

double BlackScholesCall(double Spot, double Strike, double r, double d, double Vol, double Expiry);
double BlackScholesCallVega(double Spot, double Strike, double r, double d, double Vol, double Expiry);

class BSCall
{
public:
	BSCall(double r_, double d_, double T_, double Spot_, double Strike_) 
		: r(r_), d(d_), T(T_), Spot(Spot_), Strike(Strike_) {};
	double operator()(double Vol) const { return BlackScholesCall(Spot,Strike,r,d,Vol,T); };
    double Value(double Vol) const { return operator()(Vol); };
    double Vega(double Vol) const { return BlackScholesCallVega(Spot,Strike,r,d,Vol,T); };
private:
    double r, d, T, Spot, Strike;
};

}
#endif