#ifndef HHW_H
#define HHW_H

#include <vector>
#include <complex>

namespace velesquant {


typedef std::vector<double> Vdoub;
typedef std::complex<double> Cdoub;
class HHW 
	{
	public:
		HHW(double s0, double v0, double r0, double kappa, double eta, double rho, double sigma1,
	double sigma2, double a);


		~HHW(){};


		double HHWPrice(double maturiy, double strike) const;
		Vdoub HHWPrice(Vdoub maturiy, Vdoub spot, Vdoub strike) const;
		void calibrator() const;

		Vdoub simulateHHWPrice() const;

		

	private:
		double s0_,r0_,v0_;
		double c1_,c2_,c3_;
		double kappa_;
		double eta_;
		double rho_;
		double sigma1_;
		double sigma2_;
		double a_;
	
		double HHWIntegrand(double y, double T, double K, int type) const;
		const double b(double t) const {return c1_-c2_*exp(-t*c3_);};
		/*double b(double t, double T)
			{
			double a = .5*T;
			return (c1_-c2_*exp(-a*(t+1.0)*c3_))*a;
			};
			*/
		double c(double T) const;
		const double ib(const double &t, const double &T) const;
		const double intb(const double &T) const;
	
	};





}
#endif