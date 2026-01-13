
#include <velesquant/local_vol/skewMC.h>
#include <velesquant/local_vol/sabr.h>
#include <velesquant/models/utility.h>
#include <ql/quantlib.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>
#include <algorithm>

using namespace std;
using namespace QuantLib;

namespace velesquant {

skewMC::skewMC(std::vector<double> maturities, std::vector<double> forwards, std::vector<double> betas, 
			   std::vector<double> alphas, std::vector<double> nus, std::vector<double> rhos)
{
	int N = maturities.size();
	sabrModels_.resize(N);
	for(int i=0; i<N; i++)
	{
		sabr sabrModel(maturities[i],forwards[i],betas[i],alphas[i],nus[i],rhos[i]);
		sabrModels_[i] = sabrModel;
	}
};

skewMC::skewMC(std::vector<double> maturities, std::vector<double> forwards, std::vector<double> betas, 
			   std::vector<double> alphas, std::vector<double> nus, std::vector<double> rhos, std::vector<double> shifts)
{
	int N = maturities.size();
	sabrModels_.resize(N);
	for(int i=0; i<N; i++)
	{
		sabrModels_[i] = sabr(maturities[i],forwards[i],betas[i],alphas[i],nus[i],rhos[i],shifts[i]);
	}
};

std::vector<double> skewMC::simulation(std::vector<double> timesPath, double spot, double kappa)
{
	vector<double> corrRandoms = autoCorrRandoms(timesPath,kappa);
	int N = timesPath.size();
	std::vector<double> spotsPath(N);
	for(int i=0; i<N; i++)
		spotsPath[i] = simulatedSpot(timesPath[i],spot,corrRandoms[i]);
	return spotsPath;
};

#pragma warning (disable:4715)




double skewMC::simulatedSpot(double time, double spot, double corrRN)
{
	double firstMaturity = sabrModels_[0].getMaturity();
	if(time < firstMaturity)
	{
		double nextSpot = sabrModels_[0].simulation(corrRN);
		return spot+time*(nextSpot-spot)/firstMaturity;
		double nextForward = sabrModels_[0].getForward();
		double theForward = spot+time*(nextForward-spot)/firstMaturity;
		double nextLogRet = nextSpot/nextForward-1.0;
		double theLogRet = time*nextLogRet/firstMaturity;
		double theSpot = theForward * (1.0+theLogRet);
		return theSpot;
		//double bBridge = sqrt(time)*random_normal() - time/firstMaturity * sqrt(firstMaturity)*random_normal();
		//return bBridge + (spot+time*(nextSpot-spot)/firstMaturity);
	}
	int N = sabrModels_.size();
	double lastMaturity = sabrModels_[N-1].getMaturity();
	if(time >= lastMaturity)
		return sabrModels_[N-1].simulation(corrRN);
	for(int i=1; i<N; i++)
	{
		double preMaturity = sabrModels_[i-1].getMaturity();
		double nextMaturity = sabrModels_[i].getMaturity();
		if( time >= preMaturity && time < nextMaturity)  
		{
			double preSpot = sabrModels_[i-1].simulation(corrRN);
			double nextSpot = sabrModels_[i].simulation(corrRN);
			return preSpot+(time-preMaturity)*(nextSpot-preSpot)/(nextMaturity-preMaturity);
			double preForward = sabrModels_[i-1].getForward();
			double nextForward = sabrModels_[i].getForward();
			double theForward = preForward+(time-preMaturity)*(nextForward-preForward)/(nextMaturity-preMaturity);
			double preLogRet = preSpot/preForward-1.0;
			double nextLogRet = nextSpot/nextForward-1.0;
			double theLogRet = preLogRet+(time-preMaturity)*(nextLogRet-preLogRet)/(nextMaturity-preMaturity);
			double theSpot = theForward * (1.0+theLogRet);
			return theSpot;
			//double bBridge = sqrt(time-preMaturity)*random_normal()
			//	- (time-preMaturity)/(nextMaturity-preMaturity)*sqrt(nextMaturity-preMaturity)*random_normal();
			//return bBridge + (preSpot+(time-preMaturity)*(nextSpot-preSpot)/(nextMaturity-preMaturity));
		}
	}
};

vector<double> skewMC::autoCorrRandoms(vector<double> times, double kappa) const
{
	vector<vector<double> > corrMatrix = autoCorrMatrix(times,kappa);
	vector<vector<double> > matrixCholesky = cholesky(corrMatrix);
	int N = times.size();
	vector<double> indRandoms(N),corRandoms(N);
	for(int i=0; i<N; i++)
		indRandoms[i] = random_normal();
	for(int i=0; i<N; i++)
	{
		corRandoms[i]=0.0;
		for(int j=0; j<N; j++)
			corRandoms[i] += matrixCholesky[i][j] * indRandoms[j];
	}
	return corRandoms;
};

vector<vector<double> > skewMC::autoCorrMatrix(vector<double> times, double kappa) const
{
	int N = times.size();
	vector<vector<double> > corrMatrix;
	corrMatrix.resize(N);
	for(int i=0; i<N; i++)
	{
		corrMatrix[i].resize(N);
		for(int j=0; j<i; j++)
		{
			if(kappa < 1.0E-10)
				corrMatrix[i][j] = std::sqrt(times[j]/times[i]);
			else
				corrMatrix[i][j] = std::sqrt((std::exp(2*kappa*times[j])-1)/(std::exp(2*kappa*times[i])-1));
			corrMatrix[j][i] = corrMatrix[i][j];
		}
		corrMatrix[i][i] = 1.0; 
	}
	return corrMatrix;		
};

vector<vector<double> > skewMC::cholesky(vector<vector<double> > corrMatrix) const
{
	int N = corrMatrix.size();
	Matrix AAA(N,N);
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			AAA[i][j] = corrMatrix[i][j];
	Matrix BBB(N,N);
	BBB = CholeskyDecomposition(AAA);
	vector<vector<double> > matrixCholesky;
	matrixCholesky.resize(N);
	for(int i=0; i<N; i++)
	{
		matrixCholesky[i].resize(N);
		for(int j=0; j<N; j++)
			matrixCholesky[i][j] = BBB[i][j];
	}
	return matrixCholesky;
};

/*
vector<vector<double> > skewMC::cholesky(vector<vector<double> > corrMatrix) const
{
    int N=corrMatrix.size();
	int M=corrMatrix[0].size();
	if(M != N) throw("input matrix is not a square matrix");
    for(int i=0; i<N; i++)
        for(int j=0; j<i; j++)
            if(corrMatrix[i][j] != corrMatrix[j][i]) throw("input matrix is not symmetric");
	vector<vector<double> > matrixCholesky;
	matrixCholesky.resize(N);
	for(int i=0; i<N; i++)
	{
		matrixCholesky[i].resize(N);
		for(int j=0; j<N; j++)
			matrixCholesky[i][j] = 0.0;
	}
	double sum;
	for(int i=0; i<N; i++)
		for(int j=0; j<i; j++)
		{
            sum = corrMatrix[i][j];
            for(int k=0; k<=i-1; k++) 
                sum -= matrixCholesky[i][k] * matrixCholesky[j][k];
            if(i==j) 
			{
                if(sum <= 0.0) throw("input matrix is not positive definite");
                // To handle positive semi-definite matrices take the square root of sum if positive, else zero.
                matrixCholesky[i][i] = std::sqrt(std::max(sum, 0.0));
            } 
			else 
                // With positive semi-definite matrices is possible to have matrixCholesky[i][i]==0.0 In this case sum happens to be zero as well
                matrixCholesky[j][i] = (sum==0.0 ? 0.0 : sum/matrixCholesky[i][i]);
        }
	return matrixCholesky;
};
*/

} // namespace velesquant
