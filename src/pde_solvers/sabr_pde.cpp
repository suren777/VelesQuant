#include <velesquant/pde_solvers/sabr_pde.h>
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/CyclicReduction.h>
using namespace std;
#include <algorithm>


namespace velesquant {

void sabr_pde::oneStepForward(vector<double>&inV, double dt, double &PL, double &PR)
{
	vector<double> l(sizeX_),c(sizeX_),u(sizeX_);
	vector<double>aux(inV.begin(),inV.end());
	double ndt =dt/(2*h_);
#pragma omp parallel for 



	for(int i=0; i<sizeX_; i++)
	{		
		u[i] = (i==sizeX_-1)?0:-ndt*Cm_[i+1]/(Fm_[i+1]-Fm_[i])*Em_[i+1];
		l[i] = (i==0)?0:-ndt*Cm_[i-1]/(Fm_[i]-Fm_[i-1])*Em_[i-1];
		if (i == 0) 
			c[i] = Cm_[0]/(Fm_[1]-Fm_[0])*Em_[0];
		else if (i == sizeX_-1) 
			c[i] = Cm_[i]/(Fm_[i]-Fm_[i-1])*Em_[i];
		else
			c[i] =1+ndt*Em_[i]*Cm_[i]*(1.0/(Fm_[i+1]-Fm_[i])+1.0/(Fm_[i]-Fm_[i-1]));
	}
	u[0] = Cm_[1]/(Fm_[1]-Fm_[0])*Em_[1];
	l[sizeX_-1] = Cm_[sizeX_-2]/(Fm_[sizeX_-1]-Fm_[sizeX_-2])*Em_[sizeX_-1];
	TriDiagonalSolve(sizeX_,l,c,u,inV,inV);	
	inV[0] = inV[sizeX_-1]=0;
	PL+=dt*Cm_[1]/(Fm_[1]-Fm_[0])*Em_[1]*inV[1];
	PR+=dt*Cm_[sizeX_-2]/(Fm_[sizeX_-2]-Fm_[sizeX_-1])*Em_[sizeX_-2]*inV[sizeX_-2];
};
vector<double> sabr_pde::Y(vector<double> z)
{
	vector<double> y;
	for(int i = 0; i < sizeX_; i ++) y.push_back(Y(z[i]));
	return y;
}
double	sabr_pde::Y(double z)
{
	return alpha_/nu_*(sinh(nu_*z)+rho_*(cosh(nu_*z)-1));
};
vector<double> sabr_pde::F(vector<double> z)
{
	vector<double> y;
	for(int i = 0; i < sizeX_ ; i ++) y.push_back(F(z[i]));
	return y;
}
vector<double> sabr_pde::C(vector<double> y, vector<double> f)
{
	vector<double> aux;
	for(int i = 0; i < sizeX_; i ++) aux.push_back(C(y[i],f[i]));
	return aux;
}
vector<double> sabr_pde::G(vector<double> z)
{
	vector<double> y;
	for(int i = 0; i < sizeX_; i ++) y.push_back(G(z[i]));
	return y;
}
void sabr_pde::calculateDensity()
{
	setZbounds();
	int J = sizeX_-2;
	double h0 = (zmax_-zmin_)/J;
	j0_=int(-zmin_/h0);
	h_=-zmin_/(j0_-.5);
	vector<double> z(sizeX_,0);
	zmax_=(J+1)*h_+zmin_; 	
	vector<double> ym(sizeX_,0);
	for (int i = 0; i<=J+1 ; i ++) ym[i] =Y(i*h_+zmin_- 0.5*h_);
	Fm_ = F(ym);
	Fm_[0] = 2* F(Y(zmin_))-Fm_[1]; Fm_[J+1]=2* F( Y(zmax_) )-Fm_[J];
	Cm_ = C(ym,Fm_);
	Cm_[0] = Cm_[1]; Cm_[J+1]=Cm_[J];
	vector<double> Gamma=G(Fm_);
	dT_=T_/sizeT_;
	double b = 1.0-sqrt(2)/2;
	double dt1 = dT_*b, dt2 = dT_*(1-2*b);
	Em_.resize(sizeX_);
	for (int i = 0; i < sizeX_; i++) Em_[i] = 1.0;
	vector<double> Emdt1 = Emd(dt1,Gamma); Emdt1[0] = Emdt1[1]; Emdt1[sizeX_-1]=Emdt1[sizeX_-2];
	vector<double> Emdt2 = Emd(dt2,Gamma); Emdt2[0] = Emdt2[1]; Emdt2[sizeX_-1]=Emdt2[sizeX_-2];
	QL_=QR_=0;
	double PL1,PL2,PR1,PR2;
	PL1=PL2=PR1=PR2=0;
	vector<double> inV(sizeX_,0.0);
	vector<double> inV2;
	inV[j0_]=1.0/h_;
	for (int i  = 1; i <= sizeT_; i++)
	{
		upE(Em_,Emdt1);	oneStepForward(inV,dt1, PL1, PR1);
		inV2=inV;		
		PL2=PL1,PR2=PR1;
		upE(Em_,Emdt1); oneStepForward(inV2,dt1, PL2, PR2);
		for(int i = 0; i < sizeX_ ; i++) inV[i] = (sqrt(2)+1)*inV2[i]-sqrt(2)*inV[i];
		PL1 = (sqrt(2)+1)*PL2-sqrt(2)*PL1;
		PR1 = (sqrt(2)+1)*PR2-sqrt(2)*PR1;
		upE(Em_,Emdt2);
	}
	QL_=PL1,QR_=PR1;
	Q_ = inV;
};
void sabr_pde::upE(vector<double> &inE, vector<double> &upE)
{
	for (int i=0; i < sizeX_; i++) inE[i] *=upE[i];
};
double sabr_pde::Emd(double dt, double Gamma)
{
	return exp(rho_*nu_*alpha_*Gamma*dt);
};
vector<double> sabr_pde::Emd(double dt, vector<double> Gamma)
{
	vector<double> aux(sizeX_,0);
	for (int i = 0; i < sizeX_;i++) aux[i] = Emd(dt,Gamma[i]);
	return aux;
}

double sgn(double x)
{
	return (x>0)-(x<0);
}

double	afSABR::F(double z)
{
	return pow(Lm(f_)+(1-beta_)*z,1.0/(1-beta_));
};
double	afSABR::C(double y, double f)
{
	return sqrt(alpha_*alpha_+2*rho_*alpha_*nu_*y+nu_*nu_*y*y)*L(f);
};
double	afSABR::G(double F)
{
	return (F==f_)?beta_/Lm(f_):(L(F)-L(f_))/(F-f_);
};
double	afSABR::L(double F)
{
	return pow(F+shift_,beta_);
}
double	afSABR::Lm(double F)
{
	return pow(F+shift_,1-beta_);
}
void	afSABR::setZbounds()
{
	zmax_=zmin_=-nd_*sqrt(T_);zmax_*=-1;
	if (beta_<1)
	{
		double aux = -Lm(f_)/(1-beta_);
		double zb = -1/nu_*log((sqrt(1-rho_*rho_+pow(rho_+nu_*aux/alpha_,2))-rho_-nu_*aux/alpha_)/(1-rho_));
		zmin_ = (zb>zmin_)?zb:zmin_;	
	}
}
double	afSABR::yStrike(double z)
{
	return (sgn(z)*Lm(z)-sgn(f_)*Lm(f_))/(1-beta_);
};



void antonovSABR::setZbounds()
{
	zmax_=zmin_=-nd_*sqrt(T_); 
	zmax_*=-1;
};
double antonovSABR::yStrike(double z)
{
	return (sgn(z)*Lm(z)-sgn(f_)*Lm(f_))/(1-beta_);
};
double antonovSABR::F(double z)
{
	double u = sgn(f_)*Lm(f_)+(1-beta_)*z;
	return sgn(u)*pow(abs(u),1/(1-beta_));
};
double antonovSABR::C(double y, double f)
{
	return sqrt(alpha_*alpha_+2*rho_*alpha_*nu_*y+nu_*nu_*y*y)*L(f);
};
double antonovSABR::G(double z)
{
	return (z==f_)?sgn(f_)*beta_/Lm(f_):(L(z)-L(f_))/(z-f_);
};
double antonovSABR::L(double F)
{
	return pow(abs(F),beta_);
};
double antonovSABR::Lm(double F)
{
	return pow(abs(F),1-beta_);
};
}