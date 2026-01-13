//		cTree.h

#ifndef CTREE_H
#define CTREE_H

#include <vector>

namespace velesquant {

// using namespace std; // REMOVED by refactor
typedef std::vector<double>			Vdoub;
typedef std::vector<std::vector<double>>	Mdoub;

enum pType{Call=0, Put};
enum tType{recomb=0, nonrecomb};
enum exStyle{American = 0, European, Bermudan};


class cTree
	{
	public:
		cTree();
		cTree(double S, Vdoub &T, Vdoub &F, Vdoub &IV);
		~cTree();
		double calculateBinomial(double strike, double Maturity, int Nnodes, exStyle style, pType pay, tType tree);
		double calculateBinomial(double strike, double Maturity, int Nnodes, exStyle style, pType pay, tType tree, Vdoub schedule);
	
		double calculateTrinomial(double strike, double Maturity, int Nnodes, exStyle style, pType pay, tType tree);


	private:
		double	spot_;
		double	price_;
		bool	updated_;
		Vdoub	T_, F_, IV_, r_, q_;
		//Mdoub	IV_;
		struct tTree{pType pt; tType tt;} lastCalc_;
		double interp_vol(std::vector<double> &T, std::vector<double> &F, double t1, double t2) const;
		double interpolate(std::vector<double> &T, std::vector<double> &F, double t) const;

	};
	

}
#endif