#ifndef PDESABR_H
#define PDESABR_H

#include <velesquant/models/utility.h>
#include <boost/function.hpp>
#include <vector>

namespace velesquant {

// using namespace std; // REMOVED by refactor
typedef boost::function<double(double,double)> pde_func;
typedef boost::function<double(double)> vol_func;
class pdeSABR
{
public:
	pdeSABR(double alpha, double beta, double nu, double rho, double eps, double shift, double maturity, double F,double Fmin, double Fmax,int sizeX, int sizeT): 
		alpha_(alpha),beta_(beta), rho_(rho), eps_(eps), nu_(nu),shift_(shift),T_(maturity),f_(F), Fmin_(Fmin),Fmax_(Fmax),sizeX_(sizeX), sizeT_(sizeT)
	{
		if (beta_==0)
			type_="Normal";
		else 
		{
			type_="Shifted";
		}
	};
	~pdeSABR();
	std::vector<double> calculateDensity();
	std::vector<double> getDensity()
	{
		setDensity(); 
		return Density_;
	};
	std::vector<double> getFgrid(){return Fgrid_;}
	void setDensity();
	double sabrVol(double K);
	double getAlpha(){return alpha_;};
	double getBeta(){return beta_;};
	double getNu(){return nu_;};
	double getRho(){return rho_;};
	double getShift(){return shift_;};
	
	double sabr_option(double strike, const std::string &type="Call");
	void calibrator(std::vector<volQuote> &quotes);

private:
	double eps_, alpha_, beta_, rho_, nu_,f_,dT_,h_,shift_,T_;
	double Fmax_,Fmin_, zmax_, zmin_;
	int sizeX_,sizeT_;
	double Mshifted(double F, double t);
	double Mnormal(double F,double t = 0);
	double volNormal(double K);
	double volShifted(double K);
	double volShifted1(double K);
	std::vector<volQuote>quotes_;
	double QL_, QR_;
	mutable std::vector<double> Fgrid_,Density_,Mgrid_;
	mutable std::vector<double> Fm_,Cm_,Em_;
	std::string type_;
	bool if_calibrated;
	void oneStepForward(double t, std::vector<double>&inV, pde_func &Mf);
//	void oneStepForward(double t, vector<double>&inV);
	void initialisation();	
	void objFcnCalibration(int m, int n, double* x, double* fvec, int* iflag);
	void setFbounds();
	void setZbounds();
	double Y(double z);
	std::vector<double> Y(std::vector<double> z);
	
	double F(double z);
	std::vector<double> F(std::vector<double> z);
	
	double C(double y, double f);
	std::vector<double> C(std::vector<double> y, std::vector<double> f);
	void oneStepForwardN(std::vector<double>&inV);
	double G(double z);
	std::vector<double> G(std::vector<double> z);
	double L(double F);
	double Lm(double F);
	std::vector<double> solveLawsonSwayne();
};


}
#endif