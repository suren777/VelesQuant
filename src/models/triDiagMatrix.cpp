//Solve tridiagonal system

#include <vector>
using namespace std;
#include <algorithm>

namespace velesquant {



//From Wikipedia: Tridiagonal matrix algorithm
//	a_{i} x_{i-1} + b_{i} x_{i} + c_{i} x_{i+1} = d_{i}  with b_{0}\=0, a_{0}=0, c_{n}=0
void  TriDiagonalSolve(const int n,
					   const std::vector<double> &a, 
					   const std::vector<double> &b, 
					   const std::vector<double> &c, 
					   const std::vector<double> &d, 
					   std::vector<double> &X)
{
	std::vector<double> TempC(c);	
	std::vector<double> TempD(d);
	/* Modify the coefficients. */
	TempC[0] /= b[0];	
	TempD[0] /= b[0];	
	for(int i=1; i<n; i++){
		double id = 1/(b[i]-TempC[i-1]*a[i]);  
		TempC[i] *= id;	                      
		TempD[i] = (TempD[i]-TempD[i-1]*a[i])*id;
	}
	/* Now back substitute. */
	X[n-1] = TempD[n-1];
	for(int i=n-2; i>=0; i--)
		X[i] = TempD[i]-TempC[i]*X[i+1];
}
}