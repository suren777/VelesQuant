#ifndef HHWPDE_H
#define HHWPDE_H

#include <vector>
#include <complex>

namespace velesquant {

typedef std::vector<double> Vdoub;
typedef std::complex<double> Cdoub;



class HHWPDE
	{
	public:
		HHWPDE(int Mr, int Mv, int Ms);
		~HHWPDE();








	private:
		int Mr_, Mv_, Ms_;
		int Nt_;
		int N_;
		int Spot_, K_, T_;
		double Smax_, Vmax_, Rmax_;
		double r0_;

		void setMr(int Mr){Mr_=Mr;};
		void setMv(int Mv){Mv_=Mv;};
		void setMs(int Ms){Ms_=Ms;};
		void setNt(int Nt){Nt_=Nt;};

		//Model Parameters
		
		double kappa_, eta_;
		double rho12_,rho23_,rho13_;
		double a_;
		double c1_,c2_,c3_;
		double sigma1_,sigma2_;

		double b(double t){return c1_-c2_*exp(-c3_*t);};

		//Meshes

		Vdoub sM_,rM_,vM_; // Coordinate meshes
		Vdoub bSm1, bRm1, bVm1; // beta meshes - 1
		Vdoub bS01, bR01, bV01; // beta meshes - 0
		Vdoub bSp1, bRp1, bVp1; // beta meshes + 1

		Vdoub dSm1, dRm1, dVm1; // delta meshes - 1
		Vdoub dS01, dR01, dV01; // delta meshes - 0
		Vdoub dSp1, dRp1, dVp1; // delta meshes + 1


		void aAbB(double a, Vdoub &Ain, double b, Vdoub &Bin, Vdoub &Cout);


		Vdoub dS, dR, dV; //delta meshes

		//Coefficients for finite difference schemes.
		const double alpha	(	Vdoub aVec,	int i, int type	);
		const double beta	(	Vdoub aVec,	int i, int type	);
		const double gamma	(	Vdoub aVec, int i, int type	);
		const double delta	(	Vdoub aVec, int i, int type	);
	
		const void generateMeshes();

		const double deltax	(	Vdoub aVec, int i			);


		//Mesh initialisation
		const Vdoub gsMesh();
		const Vdoub grMesh();
		const Vdoub gvMesh();
		
		double fU(Vdoub &U, int i, int j, int k){
			if ( i > Mr_ ) return 0;
			if ( j > Mv_ ) return 0;
			if ( k > Ms_ ) return 0;
			return U[i+Mv_*j+Mv_*Ms_*k];
			};

		const double extrapolate(double y1, double y2, double x1, double x2, double x)
			{
				return y1 + (x-x1)/(x2-x1)*(y2-y1);
			};


		//Matrix Vector operations A_j*U+g_j
		void A0Ug(Vdoub &Uin, Vdoub &Uout);
		void A1Ug(Vdoub &Uin, Vdoub &Uout);
		void A2Ug(Vdoub &Uin, Vdoub &Uout);
		void A3Ug(double t, Vdoub &Uin, Vdoub &Uout);

		const void fillA1(Vdoub &A1l, Vdoub &A1d, Vdoub &A1u);
		const void fillA2(Vdoub A2l1, Vdoub &A2l2, Vdoub &A2d, Vdoub &A2u1, Vdoub &A2u2);
		const void fillA3(Vdoub &A3l, Vdoub &A3d, Vdoub &A3u, double t);


		//LU vectors
		
		//Main 
		void solvePDE(Vdoub &Uin, Vdoub &Uout);

		
		const void swap(Vdoub &Uin, Vdoub &Uout){Uin.swap(Uout);};

		const void LU3dec(Vdoub &a, Vdoub &b, Vdoub &c, Vdoub &l, Vdoub &d);
		const void LU5dec(Vdoub &l, Vdoub &m, Vdoub &n, Vdoub &o, Vdoub &p,
			/*Outputs:*/Vdoub &b, Vdoub &c, Vdoub &d, Vdoub &e);	
		const void LU3decSparce(Vdoub &a, Vdoub &b, Vdoub &c, Vdoub &l, Vdoub &d, int NS);

		const void pentaSolve(Vdoub &E,Vdoub &A, Vdoub &D,Vdoub &C,Vdoub &F,Vdoub &B,Vdoub &X, int t);

		//Calculate elements
		double dudxdy(Vdoub &a, Vdoub &b,Vdoub &U, int i, int j,int k, int type);

	};




}
#endif