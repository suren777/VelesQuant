
#include <velesquant/local_vol/sabr.h>
#include <velesquant/local_vol/lVol.h>
#include <velesquant/local_vol/sVol.h>
#include <velesquant/local_vol/svolt.h>
#include <velesquant/local_vol/skewMC.h>
#include <velesquant/local_vol/schobzhu.h>
#include <velesquant/local_vol/cTree.h>
#include <velesquant/models/HW.h>
#include <velesquant/pde_solvers/HWPDE.h>
#include <velesquant/pde_solvers/ShortRate1FPDE.h>
#include <velesquant/pde_solvers/ShortRate2FPDE.h>
#include <velesquant/pde_solvers/pdeSABR.h>
#include <velesquant/pde_solvers/sabr_pde.h>
#ifndef UTILITY_H
#define UTILITY_H
#include <velesquant/models/utility.h>
#endif
#include <velesquant/models/ObjectCache.h>
#ifndef FUNCTION_H
#define FUNCTION_H
#include <boost/function.hpp>
#endif
#include <ql/quantlib.hpp>

#include <velesquant/local_vol/termstructure.h>
#include <algorithm>

using namespace std;
using namespace QuantLib;

#pragma warning (disable:4996)
#pragma warning (disable:4018)

namespace velesquant {




CellMatrix theDensity(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
					  CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, double spot, double maurity)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	std::vector<double> aDensity(modelLV->density(maurity, 100));
	delete modelLV;
	int n = aDensity.size();
	CellMatrix theDensity(n,1);
	for(int i=0; i<n; i++)
		theDensity(i,0) = aDensity[i];
	return theDensity;
}



CellMatrix lvExport(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
					CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, double spot, CellMatrix theTimes)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	int n = theTimes.RowsInStructure();
	std::vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	vector<vector<double> > LV = modelLV->exportLV(times);
	delete modelLV;
	n = LV.size();
	m = LV[0].size();
	CellMatrix theLV(n,m);
	for(int i=0; i<n; i++)
		for(int j=0; j<m; j++)
			theLV(i,j) = LV[i][j];
	return theLV;
};


double dntPDElv(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
				CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, double spot, 
				double maurity, double upperBarrier, double lowerBarrier)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	double theDNT = modelLV->dntPDE(maurity, upperBarrier, lowerBarrier);
	delete modelLV;
	return theDNT;
};

double putPDElv(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
				CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, double spot, 
				double maurity, double strike)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	double thePut = modelLV->putPDE(maurity, strike);
	delete modelLV;
	return thePut;
};

CellMatrix callPDElv(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
					 CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, double spot, 
					 double maurity, double strike)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	double theCall = modelLV->callPDE(maurity, strike);
	delete modelLV;
	return theCall;
};


CellMatrix simulationMLV(CellMatrix u1Sabr, double u1Spot,
						 CellMatrix u2Sabr, double u2Spot, CellMatrix u3Sabr, double u3Spot,
						 CellMatrix corrMatrix, CellMatrix theTimes, int NP)
{
	int m = u1Sabr.RowsInStructure();
	std::vector<double> maturities1(m);
	std::vector<double> forwards1(m);
	std::vector<double> betas1(m);
	std::vector<double> alphas1(m);
	std::vector<double> nus1(m);
	std::vector<double> rhos1(m);
	for(int i=0; i<m; i++)
	{
		maturities1[i] = u1Sabr(i,0).NumericValue();
		forwards1[i] = u1Sabr(i,1).NumericValue();
		betas1[i] = u1Sabr(i,2).NumericValue();
		alphas1[i] = u1Sabr(i,3).NumericValue();
		nus1[i] = u1Sabr(i,4).NumericValue();
		rhos1[i] = u1Sabr(i,5).NumericValue();
	}
	lVol * model1LV = new lVol(maturities1,forwards1,betas1,alphas1,nus1,rhos1,u1Spot);

	m = u2Sabr.RowsInStructure();
	std::vector<double> maturities2(m);
	std::vector<double> forwards2(m);
	std::vector<double> betas2(m);
	std::vector<double> alphas2(m);
	std::vector<double> nus2(m);
	std::vector<double> rhos2(m);
	for(int i=0; i<m; i++)
	{
		maturities2[i] = u2Sabr(i,0).NumericValue();
		forwards2[i] = u2Sabr(i,1).NumericValue();
		betas2[i] = u2Sabr(i,2).NumericValue();
		alphas2[i] = u2Sabr(i,3).NumericValue();
		nus2[i] = u2Sabr(i,4).NumericValue();
		rhos2[i] = u2Sabr(i,5).NumericValue();
	}
	lVol * model2LV = new lVol(maturities2,forwards2,betas2,alphas2,nus2,rhos2,u2Spot);

	m = u3Sabr.RowsInStructure();
	std::vector<double> maturities3(m);
	std::vector<double> forwards3(m);
	std::vector<double> betas3(m);
	std::vector<double> alphas3(m);
	std::vector<double> nus3(m);
	std::vector<double> rhos3(m);
	for(int i=0; i<m; i++)
	{
		maturities3[i] = u3Sabr(i,0).NumericValue();
		forwards3[i] = u3Sabr(i,1).NumericValue();
		betas3[i] = u3Sabr(i,2).NumericValue();
		alphas3[i] = u3Sabr(i,3).NumericValue();
		nus3[i] = u3Sabr(i,4).NumericValue();
		rhos3[i] = u3Sabr(i,5).NumericValue();
	}
	lVol * model3LV = new lVol(maturities3,forwards3,betas3,alphas3,nus3,rhos3,u3Spot);

	m = corrMatrix.RowsInStructure();
	vector<vector<double> > matrixCorr;
	matrixCorr.resize(m);
	for(int i=0; i<m; i++)
	{
		matrixCorr[i].resize(m);
		for(int j=0; j<m; j++)
			matrixCorr[i][j] = corrMatrix(i,j).NumericValue();
	}
	vector<vector<double> > matrixCholesky = cholesky(matrixCorr);

	int n = theTimes.RowsInStructure();
	std::vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();

	CellMatrix theSpots(3*n,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> rands1(n), rands2(n), rands3(n);
		for(int i=0; i<n; i++)
		{
			double r1 = random_normal();
			double r2 = random_normal();
			double r3 = random_normal();
			rands1[i] = r1;
			rands2[i] = matrixCholesky[1][0]*r1 +matrixCholesky[1][1]*r2;
			rands3[i] = matrixCholesky[2][0]*r1 +matrixCholesky[2][1]*r2 +matrixCholesky[2][2]*r3;
		}
		std::vector<double> spot1s(model1LV->simulation(times,rands1));
		for(int i=0; i<n; i++)
			theSpots(i,j) = spot1s[i];
		std::vector<double> spot2s(model2LV->simulation(times,rands2));
		for(int i=0; i<n; i++)
			theSpots(i+n,j) = spot2s[i];
		std::vector<double> spot3s(model3LV->simulation(times,rands3));
		for(int i=0; i<n; i++)
			theSpots(i+n+n,j) = spot3s[i];
	}
	delete model1LV;
	delete model2LV;
	delete model3LV;
	return theSpots;
}

CellMatrix simulation2LV(CellMatrix the1Maturities, CellMatrix the1Forwards, CellMatrix the1Betas, 
						 CellMatrix the1Alphas, CellMatrix the1Nus, CellMatrix the1Rhos, double the1Spot,
						 CellMatrix the2Maturities, CellMatrix the2Forwards, CellMatrix the2Betas, 
						 CellMatrix the2Alphas, CellMatrix the2Nus, CellMatrix the2Rhos, double the2Spot,
						 double theCorr, CellMatrix theTimes, int NP)
{
	int m = the1Maturities.RowsInStructure();
	std::vector<double> maturities1(m);
	std::vector<double> forwards1(m);
	std::vector<double> betas1(m);
	std::vector<double> alphas1(m);
	std::vector<double> nus1(m);
	std::vector<double> rhos1(m);
	for(int i=0; i<m; i++)
	{
		maturities1[i] = the1Maturities(i,0).NumericValue();
		forwards1[i] = the1Forwards(i,0).NumericValue();
		betas1[i] = the1Betas(i,0).NumericValue();
		alphas1[i] = the1Alphas(i,0).NumericValue();
		nus1[i] = the1Nus(i,0).NumericValue();
		rhos1[i] = the1Rhos(i,0).NumericValue();
	}
	lVol * model1LV = new lVol(maturities1,forwards1,betas1,alphas1,nus1,rhos1,the1Spot);

	m = the2Maturities.RowsInStructure();
	std::vector<double> maturities2(m);
	std::vector<double> forwards2(m);
	std::vector<double> betas2(m);
	std::vector<double> alphas2(m);
	std::vector<double> nus2(m);
	std::vector<double> rhos2(m);
	for(int i=0; i<m; i++)
	{
		maturities2[i] = the2Maturities(i,0).NumericValue();
		forwards2[i] = the2Forwards(i,0).NumericValue();
		betas2[i] = the2Betas(i,0).NumericValue();
		alphas2[i] = the2Alphas(i,0).NumericValue();
		nus2[i] = the2Nus(i,0).NumericValue();
		rhos2[i] = the2Rhos(i,0).NumericValue();
	}
	lVol * model2LV = new lVol(maturities2,forwards2,betas2,alphas2,nus2,rhos2,the2Spot);

	int n = theTimes.RowsInStructure();
	std::vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	CellMatrix theSpots(2*n,NP);
	for(int j=0; j<NP; j++)
	{
		std::vector<double> rands1(n), rands2(n);
		double cor2 = sqrt(1.0-theCorr*theCorr);
		for(int i=0; i<n; i++)
		{
			rands1[i] = random_normal();
			rands2[i] = theCorr*rands1[i]+cor2*random_normal();
		}
		std::vector<double> spot1s(model1LV->simulation(times,rands1));
		for(int i=0; i<n; i++)
			theSpots(i,j) = spot1s[i];
		std::vector<double> spot2s(model2LV->simulation(times,rands2));
		for(int i=0; i<n; i++)
			theSpots(i+n,j) = spot2s[i];
	}
	delete model1LV;
	delete model2LV;
	return theSpots;
}

CellMatrix simoreLV(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
					CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, CellMatrix theTimes, double spot, int NP)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	int n = theTimes.RowsInStructure();
	std::vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	CellMatrix theSpots(n,NP);
	for(int j=0; j<NP; j++)
	{
		std::vector<double> spots(modelLV->simulation(times));
		for(int i=0; i<n; i++)
			theSpots(i,j) = spots[i];
	}
	delete modelLV;
	return theSpots;
}

CellMatrix simulationLV(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
						CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, CellMatrix theTimes, double spot)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	lVol * modelLV = new lVol(maturities,forwards,betas,alphas,nus,rhos,spot);
	int n = theTimes.RowsInStructure();
	std::vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	CellMatrix theSpots(n,1);
	std::vector<double> spots(modelLV->simulation(times));
	for(int i=0; i<n; i++)
		theSpots(i,0) = spots[i];
	delete modelLV;
	return theSpots;
}


std::string hestonSetSeed(const std::string &theName, int seed)
{
	if (seed<0) throw("Invalid seed");
	getObject<sVol>(theName)->update_seed(seed);
	std::string status = "OK\n";
	return status;
}


CellMatrix hestonSimulation(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, int NP)
{
	sVol *h = getObject<sVol>(theName);

	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHeston(times,forwards);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

double hestonSimulationDO(const std::string &theName, 
						  CellMatrix theTimes, 
						  CellMatrix theForwards, 
						  double barrier, 
						  double strike, 
						  int NP, 
						  int seed, 
						  const std::string &type)
{
	sVol *h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times, forwards;
	hestonSetSeed(theName,seed);
	for(int i=0; i<m; i++)
	{	
		if (!theTimes(i,0).IsANumber()) break;
		if (theTimes(i,0).NumericValue()>0.0)
		{
			times.push_back(theTimes(i,0).NumericValue());
			forwards.push_back(theForwards(i,0).NumericValue());
		}
		else
			break;
	}
	if (times.size()<1) throw("Incorrect input");
	double phi, payoff = 0.0;
	if (type == "Put") phi = -1.0; else phi = 1.0;

	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHestonDO(times,forwards,barrier);
		payoff += max(0.0,phi*(pathF[0]-strike))*pathF[1];
	}
	return payoff/double(NP);
}


CellMatrix hestonCalibrator(const std::string &theName, 
							CellMatrix theMaturitys, CellMatrix theForwards, CellMatrix theStrikes, 
							CellMatrix theQuotes, string quoteType)
{
	sVol *h = getObject<sVol>(theName);
	int m = theStrikes.RowsInStructure();
	int k = theQuotes.RowsInStructure();
	int l = theQuotes.ColumnsInStructure();
	std::vector<double> maturitys,forwards,strikes,quotes;
	if ((m == k)&&(l==1))
	{
		for(int i=0; i<m; i++)
		{
			maturitys.push_back(theMaturitys(i,0).NumericValue());
			forwards.push_back(theForwards(i,0).NumericValue());
			strikes.push_back(theStrikes(i,0).NumericValue());
			quotes.push_back(theQuotes(i,0).NumericValue());
		}
	}
	else 
	{
		m = theStrikes.ColumnsInStructure();
		if (m!=l) throw ("Error! Strikes dont match vols");

		for (int i = 0 ; i < k ; i++)
			for(int j = 0; j < l ; j++)
				if (theQuotes(i,j).NumericValue()>1e-10)
				{
					maturitys.push_back(theMaturitys(i,0).NumericValue());
					forwards.push_back(theForwards(i,0).NumericValue());
					strikes.push_back(theStrikes(0,j).NumericValue());
					quotes.push_back(theQuotes(i,j).NumericValue());
				}
	}
	h->calibrator(maturitys,forwards,strikes,quotes,quoteType);
	CellMatrix hestonPara(5,1);
	hestonPara(0,0) = h->getParameterVar0();
	hestonPara(1,0) = h->getParameterKappa();
	hestonPara(2,0) = h->getParameterTheta();
	hestonPara(3,0) = h->getParameterXi();
	hestonPara(4,0) = h->getParameterRho();
	return hestonPara;
}

CellMatrix hestonCalibratorIV(const std::string &theName, 
							  CellMatrix theMaturitys, CellMatrix theForwards, CellMatrix theStrikes, 
							  CellMatrix theQuotes)
{
	sVol *h = getObject<sVol>(theName);
	int m = theStrikes.RowsInStructure();
	std::vector<double> maturitys(m);
	std::vector<double> forwards(m);
	std::vector<double> strikes(m);
	std::vector<double> quotes(m);
	for(int i=0; i<m; i++)
	{
		maturitys[i] = theMaturitys(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		strikes[i] = theStrikes(i,0).NumericValue();
		quotes[i] = theQuotes(i,0).NumericValue();;
	}
	h->IVcalibrator(maturitys,forwards,strikes,quotes,"0");
	CellMatrix hestonPara(5,1);
	hestonPara(0,0) = h->getParameterVar0();
	hestonPara(1,0) = h->getParameterKappa();
	hestonPara(2,0) = h->getParameterTheta();
	hestonPara(3,0) = h->getParameterXi();
	hestonPara(4,0) = h->getParameterRho();
	return hestonPara;
}

CellMatrix hestonCalibratorFX(const std::string &theName, 
							  CellMatrix theMaturitys, CellMatrix theForwards, CellMatrix theStrikes, 
							  CellMatrix theQuotes)
{
	sVol *h = getObject<sVol>(theName);
	int m = theStrikes.RowsInStructure();
	std::vector<double> maturitys(m);
	std::vector<double> forwards(m);
	std::vector<double> strikes(m);
	std::vector<double> quotes(m);
	for(int i=0; i<m; i++)
	{
		maturitys[i] = theMaturitys(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		strikes[i] = theStrikes(i,0).NumericValue();
		quotes[i] = theQuotes(i,0).NumericValue();;
	}
	h->FXcalibrator(maturitys,forwards,strikes,quotes,"0");
	CellMatrix hestonPara(5,1);
	hestonPara(0,0) = h->getParameterVar0();
	hestonPara(1,0) = h->getParameterKappa();
	hestonPara(2,0) = h->getParameterTheta();
	hestonPara(3,0) = h->getParameterXi();
	hestonPara(4,0) = h->getParameterRho();
	return hestonPara;
}


std::string CreateHestonObj(const std::string &theObjName, 
							double spot, double var0, double kappa, double theta, double xi, double rho, int seed)
{
	//sVol* theObj = new sVol(spot, var0, kappa, theta, xi, rho, seed);

	return placeObject(theObjName,new sVol(spot, var0, kappa, theta, xi, rho, seed));
};

double heston(double maturity, double forward, double strike,   
			  double spot, double var0, double kappa, double theta, double xi, double rho)
{
	sVol* theObj = new sVol(spot,var0,kappa,theta,xi,rho,-1);
	double result = theObj->hestonPrice(maturity,forward,strike,"call");
	delete theObj;
	return result;
}
double hestonPrice(std::string theName, double maturity, double forward, double strike)
{
	return getObject<sVol>(theName)->hestonPriceCF(maturity,forward,strike,"call");
}
double hestonPriceCF(std::string theName, double maturity, double forward, double strike)
{
	return getObject<sVol>(theName)->hestonPriceCF(maturity,forward,strike,"call");
}





CellMatrix simoreSMC(CellMatrix theMaturities, CellMatrix theForwards, CellMatrix theBetas, 
					 CellMatrix theAlphas, CellMatrix theNus, CellMatrix theRhos, CellMatrix theTimes, 
					 double spot, double kappa, int NP)
{
	int m = theMaturities.RowsInStructure();
	std::vector<double> maturities(m);
	std::vector<double> forwards(m);
	std::vector<double> betas(m);
	std::vector<double> alphas(m);
	std::vector<double> nus(m);
	std::vector<double> rhos(m);
	for(int i=0; i<m; i++)
	{
		maturities[i] = theMaturities(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		betas[i] = theBetas(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
		nus[i] = theNus(i,0).NumericValue();
		rhos[i] = theRhos(i,0).NumericValue();
	}
	skewMC * modelSMC = new skewMC(maturities,forwards,betas,alphas,nus,rhos);
	int n = theTimes.RowsInStructure();
	std::vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	CellMatrix theSpots(n,NP);
	for(int j=0; j<NP; j++)
	{
		std::vector<double> spots(modelSMC->simulation(times,spot,kappa));
		for(int i=0; i<n; i++)
			theSpots(i,j) = spots[i];
	}
	delete modelSMC;
	return theSpots;
}


CellMatrix skewSimulation(const std::string &theName, CellMatrix theCRNs)
{
	int N = theCRNs.RowsInStructure();
	std::vector<double> correlatedRNs(N);
	for(int i=0; i<N; i++)
		correlatedRNs[i] = theCRNs(i,0).NumericValue();
	std::vector<double> spots = getObject<sabr>(theName)->simulations(correlatedRNs);
	CellMatrix theSpots(N,1);
	for(int i=0; i<N; i++)
		theSpots(i,0) = spots[i];
	return theSpots;
}

CellMatrix generateQUTL(const std::string &theName)
{
	sabr *s = getObject<sabr>(theName);
	s->qutlTable();
	std::vector<double> theSpots = s->getTableQUTL1();
	std::vector<double> theQUTLs = s->getTableQUTL2();
	int N = theQUTLs.size();
	CellMatrix tableQUTL(N,2);
	for(int i=0;i<N;i++)
	{
		tableQUTL(i,0) = theSpots[i];
		tableQUTL(i,1) = theQUTLs[i];
	}
	return tableQUTL;
}

//CellMatrix sabrCalibratorWithInitial(const std::string &theName, CellMatrix theStrikes, CellMatrix theQuotes, 
//									 std::string quoteType, CellMatrix theParams)
//{
//	std::string theKey(stripTrailingHash(theName));
//	if(cachedObject<sabr>::instance().find(theKey)==cachedObject<sabr>::instance().end()) 
//		throw("Object sabr not found in Cache");
//	int m = theStrikes.RowsInStructure();
//	std::vector<double> strikes(m);
//	std::vector<double> quotes(m);
//	for(int i=0; i<m; i++)
//	{
//		strikes[i] = theStrikes(i,0).NumericValue();
//		quotes[i] = theQuotes(i,0).NumericValue();
//	}
//	std::vector<double> params(3);
//	params[0] = theParams(0,0).NumericValue();
//	params[1] = theParams(1,0).NumericValue();
//	params[2] = theParams(2,0).NumericValue();
//	cachedObject<sabr>::instance()[theKey]->calibratorWithInitial(strikes,quotes,quoteType,params);
//	CellMatrix sabrPara(3,1);
//	sabrPara(0,0) = cachedObject<sabr>::instance()[theKey]->getParameterAlpha();
//	sabrPara(1,0) = cachedObject<sabr>::instance()[theKey]->getParameterNu();
//	sabrPara(2,0) = cachedObject<sabr>::instance()[theKey]->getParameterRho();
//	return sabrPara;
//}


CellMatrix sabrCalibratorWithInitial(const std::string &theName, CellMatrix theStrikes, CellMatrix theQuotes, 
									 std::string quoteType,std::string atmFlag)
{
	sabr *s = getObject<sabr>(theName);
	int m = theStrikes.RowsInStructure();

	std::vector<double> strikes(m);
	std::vector<double> quotes(m);
	double forward = s-> getForward();
	double T = s-> getMaturity();

	for(int i=0; i<m; i++)
	{
		strikes[i] = theStrikes(i,0).NumericValue();
		quotes[i] = theQuotes(i,0).NumericValue();

		//get ATM vol
		if (strikes[i]/ forward == 1.0)
		{
			if (quoteType == "premium" && atmFlag =="Yes" )
			{
				double atmVol = implied_vol(T, forward, strikes[i], quotes[i]);
				s-> setATMvol(atmVol);
			}
			else
			{
				s-> setATMvol(quotes[i]);
			}
		}
	}

	if (atmFlag =="Yes")
	{
		s->calibratorWithInitialATM(strikes,quotes,quoteType);
	}
	else
	{
		s->calibratorWithInitial(strikes,quotes,quoteType);
	}

	CellMatrix sabrPara(3,1);
	sabrPara(0,0) = s->getParameterAlpha();
	sabrPara(1,0) = s->getParameterNu();
	sabrPara(2,0) = s->getParameterRho();
	return sabrPara;
}

CellMatrix sabrCalibrator(const std::string &theName, CellMatrix theStrikes, CellMatrix theQuotes, 
						  std::string quoteType)
{
	sabr *s = getObject<sabr>(theName);
	int m = theStrikes.RowsInStructure();
	int k=0;
	std::vector<double> strikes(m);
	std::vector<double> quotes(m);
	for(int i=0; i<m; i++)
	{
		if (theStrikes(i,0).IsANumber()&&theQuotes(i,0).IsANumber())
		{
			k++;
			strikes[k-1] = theStrikes(i,0).NumericValue();
			quotes[k-1] = theQuotes(i,0).NumericValue();
		};
	}
	if(k==0) throw("Not enough data");
	strikes.resize(k);
	quotes.resize(k);
	s->calibrator(strikes,quotes,quoteType);
	CellMatrix sabrPara(3,1);
	sabrPara(0,0) = s->getParameterAlpha();
	sabrPara(1,0) = s->getParameterNu();
	sabrPara(2,0) = s->getParameterRho();
	return sabrPara;
}

std::string CreateSABRObj(const std::string &theObjName, 
						  double maturity, double forward, double beta, 
						  double alpha, double nu, double rho)
{
	//sabr* theObj = new sabr(maturity, forward, beta, alpha, nu, rho);
	return placeObject(theObjName,new sabr(maturity, forward, beta, alpha, nu, rho));
}
std::string dvlCreateSABRObjShifted(const std::string &theObjName, 
						  double maturity, double forward, double beta, 
						  double alpha, double nu, double rho, double shift)
{
	//sabr* theObj = new sabr(maturity, forward, beta, alpha, nu, rho);
	return placeObject(theObjName,new sabr(maturity, forward, beta, alpha, nu, rho, shift));
}

double impliedVolatility(const std::string &theName, double theStrike)
{
	return getObject<sabr>(theName)->getVol(theStrike);
}

double forwardPremium(const std::string &theName, double theStrike, std::string callORput)
{
	return getObject<sabr>(theName)->getPremium(theStrike,callORput);
}

double localVolatility(const std::string &theName, double theSpot)
{
	return getObject<sabr>(theName)->localVol(theSpot);
}

double localVolatilityzabr(const std::string &theName, double theSpot)
{
	return getObject<sabr>(theName)->localVolzabr(theSpot);
}

double sabrObj(double forward, double maturity, 
			   double alpha, double beta, double nu, double rho, double strike)
{
	sabr example(maturity, forward, beta, alpha, nu, rho);
	return example.impliedVol(strike);
}

void sabrParameters(double alpha, double beta, double nu, double rho) 
{
	if( !((alpha>0.0) && (beta>=0.0 && beta<=1.0) && (nu>=0.0) && (rho*rho<=1.0)) )
		throw("SABR parameters are not valid");
}
double sabrVolatility(double strike, double forward, double maturity, 
					  double alpha, double beta, double nu, double rho) 
{
	sabrParameters(alpha, beta, nu, rho);
	const double mean = std::pow(forward*strike, 0.5*(1.0-beta));
	const double logR = std::log(forward/strike);
	const double tem = (1.0-beta)*(1.0-beta)*logR*logR;
	const double dem = (1.0+tem/24.0+tem*tem/1920.0)*mean;
	const double coe = 1.0 + maturity * ( (1.0-beta)*(1.0-beta)*alpha*alpha/(mean*mean)/24.0
		+ rho*beta*nu*alpha/mean/4.0 + (2.0-3.0*rho*rho)*nu*nu/24.0 );
	const double z = (nu/alpha)*mean*logR;
	double multiplier;
	if (std::fabs(z*z)>1.0E-20) {
		const double xz = std::log((std::sqrt(1.0-2.0*rho*z+z*z)+z-rho)/(1.0-rho));
		multiplier = z/xz;
	}
	else {
		multiplier = 1.0 - 0.5*rho*z - (3.0*rho*rho-2.0)*z*z/12.0;
	}
	return (alpha/dem)*multiplier*coe;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


CellMatrix dvlHestonMC_TD(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, int NP)
{
	sVolT *theObj = getObject<sVolT>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = theObj->simulationHestonTd(times,forwards);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}





CellMatrix dvlHestonCalibratorTD(const std::string &theName, 
								 CellMatrix theMaturitys, CellMatrix theForwards, CellMatrix theStrikes, 
								 CellMatrix theQuotes, string quoteType)
{
	sVolT *theObj = getObject<sVolT>(theName);
	int m = theStrikes.RowsInStructure();
	int k = theQuotes.RowsInStructure();
	int l = theQuotes.ColumnsInStructure();
	std::vector<double> maturitys,forwards,strikes,quotes;
	if ((m == k)&&(l==1))
	{
		for(int i=0; i<m; i++)
		{
			maturitys.push_back(theMaturitys(i,0).NumericValue());
			forwards.push_back(theForwards(i,0).NumericValue());
			strikes.push_back(theStrikes(i,0).NumericValue());
			quotes.push_back(theQuotes(i,0).NumericValue());
		}
	}
	else 
	{
		m = theStrikes.ColumnsInStructure();
		if (m!=l) throw ("Error! Strikes dont match vols");

		for (int i = 0 ; i < k ; i++)
			for(int j = 0; j < l ; j++)
				if (theQuotes(i,j).NumericValue()>1e-10)
				{
					maturitys.push_back(theMaturitys(i,0).NumericValue());
					forwards.push_back(theForwards(i,0).NumericValue());
					strikes.push_back(theStrikes(0,j).NumericValue());
					quotes.push_back(theQuotes(i,j).NumericValue());
				}
	}
	theObj->calibratorLM(maturitys,forwards,strikes,quotes,quoteType);
	vector<double> t, xi, theta;
	t = theObj->getT();
	xi = theObj->getXi();
	theta = theObj->getTheta();
	CellMatrix hestonPara(t.size(),6);
	hestonPara(0,0) = theObj->getParameterVar0();
	hestonPara(0,1) = theObj->getParameterKappa();
	hestonPara(0,2) = theObj->getParameterRho();
	for( int i = 0; i < t.size(); i++) 		
	{
		hestonPara(i,3)=t[i];	hestonPara(i,5)=xi[i];	hestonPara(i,4)=theta[i];	
	}


	return hestonPara;
}

std::string dvlCreateHestonTDObj(const std::string &theObjName, 
								 double spot, double var0, double kappa, double rho, MyArray time, MyArray theta, MyArray xi, int seed)
{
	if (time.size()<=1) throw("Wrong Model Used");
	sVolT* theObj;
	if ((theta.size() == 1)&&(xi.size()>=1))
		theObj = new sVolT(spot, var0, kappa, rho, time,theta[1], xi, seed);
	else if ((theta.size()>= 1)&&(xi.size()==1))
		theObj = new sVolT(spot, var0, kappa, rho, time,theta, xi[1], seed);
	else if ((theta.size()>= 1)&&(xi.size()>=1))
		theObj = new sVolT(spot, var0, kappa, rho, time,theta, xi, seed);
	else throw("Something went wrong");
	return placeObject(theObjName,theObj);
};

double dvlPriceHestonCallTD(double maturity, double forward, double strike,   
							double spot, double var0, double kappa, double rho, MyArray time, MyArray theta, MyArray xi)
{
	if (time.size()<=1) throw("Wrong Model Used");
	sVolT* theObj;
	if ((theta.size() == 1)&&(xi.size()>=1))
		theObj = new sVolT(spot, var0, kappa, rho, time,theta[1], xi,	-1);
	else if ((theta.size()>= 1)&&(xi.size()==1))
		theObj = new sVolT(spot, var0, kappa, rho, time,theta,	xi[1],	-1);
	else if ((theta.size()>= 1)&&(xi.size()>=1))
		theObj = new sVolT(spot, var0, kappa, rho, time,theta,	xi,		-1);
	else throw("Something went wrong");
	double result = theObj->hestonPriceTd(maturity,forward,strike,"call");
	delete theObj;
	return result;
}

double dvlPriceHestonCallTDo(const std::string &theObjName,double maturity, double forward, double strike)
{
	return getObject<sVolT>(theObjName)->hestonPriceTdCF(maturity,forward,strike,"call");
}

///Cliquet stuff

CellMatrix hestonSimulationTimeDependentMax(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, int NP)
{
	sVolT* myObj = getObject<sVolT>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = myObj->simulationHestonTdMax(times,forwards);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

CellMatrix hestonSimulationMax(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHestonMax(times,forwards);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

CellMatrix dvlHestonMC_CliquetTD(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, 
								 double gcap, double gfloor, double lcap, double lfloor, double alpha, int NP)
{
	sVolT* myObj = getObject<sVolT>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = myObj->simulationHestonTdCliq(times,forwards,gcap,gfloor,lcap,lfloor,alpha);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

CellMatrix hestonSimulationCliq(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double gcap, double gfloor, double lcap, double lfloor, double alpha, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHestonCliq(times,forwards,gcap,gfloor,lcap,lfloor,alpha);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

CellMatrix hestonSimulationDNT(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double UP, double DOWN, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHestonDNT(times,forwards,UP,DOWN);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

CellMatrix hestonSimulationDNTdt(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double UP, double DOWN, double dt, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(1,NP);
	for(int j=0; j<NP; j++)
	{		
		thePaths(0,j) = h->simulationHestonDNTdt(times,forwards,UP,DOWN,dt);
	}
	return thePaths;
}

double hestonSimulationDNTdtS(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double maturity, double UP, double DOWN, double dt, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	long int sum=0;
	for(int j=0; j<NP; j++)
	{		
		sum +=h->simulationHestonDNTdtS(times,forwards,maturity,UP,DOWN,dt);
	}
	return double(sum)/NP;
}
double hestonSimulationDNTdtE(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double maturity, double UP, double DOWN, double dt, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	long int sum=0;
	for(int j=0; j<NP; j++)
	{		
		sum += h->simulationHestonDNTdtE(times,forwards,maturity,UP,DOWN,dt);
	}
	return double(sum)/NP;
}


CellMatrix hestonSimulationDownnOut(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double BARRIER, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHestonDownnOut(times,forwards,BARRIER);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

CellMatrix hestonSimulationUpnOut(const std::string &theName, CellMatrix theTimes, CellMatrix theForwards, double BARRIER, int NP)
{
	sVol* h = getObject<sVol>(theName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = h->simulationHestonUpnOut(times,forwards,BARRIER);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}



std::string CreateTermsStructureObj(const std::string &theObjName, 
									MyArray days, MyArray rate, MyArray qdays, MyArray divident, int calendar, int daycount) 
{
	return placeObject(theObjName, new Termstructure(days, rate, qdays, divident, calendar, daycount));

};

double forwardRate(const std::string &theName, int day, int tenor)
{
	return getObject<Termstructure>(theName)->rate(day,tenor);
}


double divident(const std::string &theName, int day, int tenor)
{
	return getObject<Termstructure>(theName)->divident(day,tenor);
}



CellMatrix hestonSimTDepCliqTerm(const std::string &theName, const std::string &termStruct, int valDay, CellMatrix theTimes, CellMatrix theForwards, 
								 double gcap, double gfloor, double lcap, double lfloor, double alpha, int NP)
{


	Termstructure *terme = getObject<Termstructure>(termStruct);
	sVolT * s = getObject<sVolT>(theName);

	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = s->simulationHestonTdCliq(times,forwards,gcap,gfloor,lcap,lfloor,alpha,valDay, terme);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}
////Schobe Zhu
CellMatrix szSimulation(const std::string &theSZObjectName, CellMatrix theTimes, CellMatrix theForwards, int NP)
{
	SchobZhu *s=getObject<SchobZhu>(theSZObjectName);
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m), forwards(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
	}
	CellMatrix thePaths(m,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> pathF = s->simulationSchobZhu(times,forwards);
		for(int i=0; i<m; i++)
			thePaths(i,j) = pathF[i];
	}
	return thePaths;
}

double szModel(double maturity, double forward, double strike,   
			   double spot, double var0, double kappa, double theta, double xi, double rho)
{
	SchobZhu* theObj = new SchobZhu(spot,var0,kappa,theta,xi,rho);
	double result = theObj->SchobelPrice(maturity,forward,strike);
	delete theObj;
	return result;
}

std::string CreateSZObj(const std::string &theSZObjName, 
						double spot, double var0, double kappa, double theta, double xi, double rho)
{
	SchobZhu* theObj = new SchobZhu(spot, var0, kappa, theta, xi, rho);	
	return placeObject(theSZObjName,theObj);
};

CellMatrix szCalibrator(const std::string &theSZObjectName, 
						CellMatrix theMaturitys, CellMatrix theForwards, CellMatrix theStrikes, 
						CellMatrix theQuotes)
{
	SchobZhu *s=getObject<SchobZhu>(theSZObjectName);
	int m = theStrikes.RowsInStructure();
	std::vector<double> maturitys(m);
	std::vector<double> forwards(m);
	std::vector<double> strikes(m);
	std::vector<double> quotes(m);
	for(int i=0; i<m; i++)
	{
		maturitys[i] = theMaturitys(i,0).NumericValue();
		forwards[i] = theForwards(i,0).NumericValue();
		strikes[i] = theStrikes(i,0).NumericValue();
		quotes[i] = theQuotes(i,0).NumericValue();;
	}
	s->calibrator(maturitys,forwards,strikes,quotes);
	CellMatrix hestonPara(5,1);
	hestonPara(0,0) = s->getParameterVar0();
	hestonPara(1,0) = s->getParameterKappa();
	hestonPara(2,0) = s->getParameterTheta();
	hestonPara(3,0) = s->getParameterXi();
	hestonPara(4,0) = s->getParameterRho();
	return hestonPara;
}

double dvlLin_interp(CellMatrix rngX, CellMatrix rngY, double Xi) 
{

	vector<double> X, Y;
	//Error checking
	int i,k=0;
	int m = rngX.RowsInStructure();
	if (m!= rngY.RowsInStructure()) throw("Error. Input data is not equal in size");
	int n = rngX.ColumnsInStructure();
	if (n!= rngY.ColumnsInStructure()) throw("Error. Input data is not equal in size");

	if (n>m)
		for(i=0; i<n; i++)
		{
			if (rngX(0,i).IsANumber()&&rngY(0,i).IsANumber())
			{
				k++;
				X.push_back(rngX(0,i).NumericValue());
				Y.push_back(rngY(0,i).NumericValue());
			};
		}
	else
		for(int i=0; i<m; i++)
		{
			if (rngX(i,0).IsANumber()&&rngY(i,0).IsANumber())
			{
				k++;
				X.push_back(rngX(i,0).NumericValue());
				Y.push_back(rngY(i,0).NumericValue());
			};
		}
		if(k==0) throw("Not enough data");

		LinearInterpolation interp(X.begin(),X.end(),Y.begin());
		return interp(Xi, true);
}



double t_func(double x)
{
	return exp(2.0*x)*sin(3.0*x);
}




double my_integral(double a, double b, int N)
{
	double f1,fN,fsum;
	double h = (b-a)/double(N);
	fsum=0.0;
	f1 = t_func(a);
	fN = t_func(b);
	for(int i = 1; i<N+1;i++)
		fsum+=t_func(a+i*h);
	return h/2.0*(f1+fN+2*fsum);
}


double my_exact(double a, double b)
{
	double Fb = exp(2.0*b)*(2.0*sin(3.0*b)-3.0*cos(3.0*b))/13.0;
	double Fa = exp(2.0*a)*(2.0*sin(3.0*a)-3.0*cos(3.0*a))/13.0;
	return Fb-Fa;
}

double black76(double forward, double strike, double maturity, double IV)
{
	double vol = IV* sqrt(maturity);
	double d1 = std::log(forward/strike)/vol +0.5*vol;
	double d2 = d1-vol;
	return forward*cdf_normal(d1) - strike*cdf_normal(d2);
}


std::string CreateTreeObj(const std::string &theTreeObjName, double spot,
						  CellMatrix maturities, CellMatrix forwards, CellMatrix vols)
{
	int N = maturities.RowsInStructure();
	Vdoub T(N), F(N), IV(N);
	for (int i = 0; i<N; i++)
	{
		T[i] = maturities(i,0).NumericValue();
		F[i] = forwards(i,0).NumericValue();
		IV[i] = vols(i,0).NumericValue();
	}
	return placeObject(theTreeObjName,new cTree(spot,  T,  F,  IV));
};

double TreePriceBinomialAmerican(std::string theName, double strike, double Maturity, int Nnodes,std::string payoff, std::string treeType)
{
	pType pay = (payoff=="Call")?Call:Put;
	tType tree = (treeType == "Recombining")?recomb:nonrecomb;
	return getObject<cTree>(theName)->calculateBinomial(strike, Maturity, Nnodes, American, pay, tree);
}

double TreePriceBinomialEuropean(std::string theName, double strike, double Maturity, int Nnodes,std::string payoff, std::string treeType)
{

	pType pay = (payoff=="Call")?Call:Put;
	tType tree = (treeType == "Recombining")?recomb:nonrecomb;
	return getObject<cTree>(theName)->calculateBinomial(strike, Maturity, Nnodes, European, pay, tree);
}
double TreePriceTrinomialEuropean(std::string theName, double strike, double Maturity, int Nnodes,std::string payoff, std::string treeType)
{
	pType pay = (payoff=="Call")?Call:Put;
	tType tree = (treeType == "Recombining")?recomb:nonrecomb;
	return getObject<cTree>(theName)->calculateTrinomial(strike, Maturity, Nnodes, European, pay, tree);
}

double TreePriceTrinomialAmerican(std::string theName, double strike, double Maturity, int Nnodes,std::string payoff, std::string treeType)
{
	pType pay = (payoff=="Call")?Call:Put;
	tType tree = (treeType == "Recombining")?recomb:nonrecomb;
	return getObject<cTree>(theName)->calculateTrinomial(strike, Maturity, Nnodes, American, pay, tree);
}



//Converts DF curve to spot rates
CellMatrix DFtoR(CellMatrix theDF, CellMatrix theT, int n)
{
	int m = theDF.RowsInStructure();
	CellMatrix R0s(m,1);

	std::vector<double> DF(m);
	std::vector<double> T(m);

	for(int i=0; i<m; i++)
	{
		DF[i] = theDF(i,0).NumericValue();
		T[i] = theT(i,0).NumericValue();
	}

	R0s = DFtoR(DF,T, n);

	return R0s;
};

//Calculates Div yield based on forward prices, spot and DF curve
CellMatrix DFtoDiv(CellMatrix theDF, CellMatrix thedfT, CellMatrix theFwd, CellMatrix thefwdT, double Spot, int n)
{
	int m1 = theDF.RowsInStructure();
	int m2 = theFwd.RowsInStructure();
	CellMatrix Divs(m2,1);

	std::vector<double> DF(m1), dfT(m1), Fwd(m2), fwdT(m2);

	for(int i=0; i<m1; i++)
	{
		DF[i] = theDF(i,0).NumericValue();
		dfT[i] = thedfT(i,0).NumericValue();
	}

	for(int i=0; i<m2; i++)
	{
		Fwd[i] = theFwd(i,0).NumericValue();
		fwdT[i] = thefwdT(i,0).NumericValue();
	}

	Divs = DFtoDiv(DF, dfT, Fwd, fwdT, Spot, n);

	return Divs;
};

//Calculates forward rates based on the DF curve
CellMatrix DFtoFwd(CellMatrix theDF,CellMatrix theT, double delta, int n)
{
	int m = theDF.RowsInStructure();
	CellMatrix Fwds(m,1);

	std::vector<double> DF(m), T(m);

	for(int i=0; i<m; i++)
	{
		DF[i] = theDF(i,0).NumericValue();
		T[i] = theT(i,0).NumericValue();
	}

	Fwds = DFtoFwd(DF, T, delta, n);

	return Fwds;
};

// Calculates the annuity for swaption
double annuity(double Expiry, double Tenor, double Freq, CellMatrix theDF)
{
	int m = theDF.RowsInStructure();
	std::vector<double> DF(m), T(m);

	for(int i=0; i<m; i++)
	{
		T[i] = theDF(i,0).NumericValue();
		DF[i] = theDF(i,1).NumericValue();
	}

	return annuity(Expiry,Tenor,Freq,DF,T);

};

CellMatrix MClogVol(double spot, double freq,  CellMatrix theTimes, CellMatrix theSims)
{
	int n = theSims.RowsInStructure();
	int m = theSims.ColumnsInStructure();
	vector<vector<double>> logReturns(m,vector<double>(n));

	vector<double> impVol(m), times(m);

	//Initial time step
	times[0] = theTimes(0,0).NumericValue(); 

	for (int j=0; j<n;j++)
	{
		logReturns[0][j] = log(theSims(j,0).NumericValue()/spot);
	}
	//Remaining time steps
	for(int i=1;i<m;i++)
	{
		times[i] = theTimes(0,i).NumericValue();

		for (int j=0; j<n;j++)
		{
			logReturns[i][j] = log(theSims(j,i).NumericValue()/theSims(j,i-1).NumericValue());
		}
	}
	impVol = MCimpVol(freq,times,logReturns);
	return impVol;
};

// Calculates the forward swap rate for swaption
double fwdSwapRate(double Expiry, double Tenor, double Freq, CellMatrix theDF)
{
	int m = theDF.RowsInStructure();
	std::vector<double> DF(m), T(m);

	for(int i=0; i<m; i++)
	{
		T[i] = theDF(i,0).NumericValue();
		DF[i] = theDF(i,1).NumericValue();
	}

	return fwdSR(Expiry,Tenor,Freq,DF,T);

};


string dvlCreateHWObj(const string &theObjName, double kappa, CellMatrix& theTimesSigmas, CellMatrix& theTimesDFs)
{
	int m = theTimesSigmas.RowsInStructure();
	vector<double> timeSigmas, sigmas;
	for(int i=0; i<m; i++)
		if (theTimesSigmas(i,0).IsANumber()&&theTimesSigmas(i,0).NumericValue()!=0)
		{
			timeSigmas.push_back( theTimesSigmas(i,0).NumericValue());
			sigmas.push_back(theTimesSigmas(i,1).NumericValue());
		}
		m = theTimesDFs.RowsInStructure();
		std::vector<double> timeDFs(m), DFs(m);
		for(int i=0; i<m; i++)
		{
			timeDFs[i] = theTimesDFs(i,0).NumericValue();
			DFs[i] = theTimesDFs(i,1).NumericValue();
		}
		return placeObject(theObjName,new HW(kappa, timeSigmas, sigmas, timeDFs, DFs));
};
double dvlSwaptionHW(const string &theName, double Expiry, double Tenor, double Strike, double PayFrequency)
{
	return getObject<HW>(theName)->swaption(Expiry, Tenor, Strike, PayFrequency);
};


CellMatrix dvlCalibrateHW(const string &theName, CellMatrix theStrikes, CellMatrix theFreq, CellMatrix theSurface, string &type)
{
	HW* myObj = getObject<HW>(theName);
	int n = theSurface.RowsInStructure();
	int m = theSurface.ColumnsInStructure();

	defSwap theQuote;
	vector<defSwap> quoteSwap;

	//throw if matrices not the same

	//loop through rows
	for(int i=1;i<n;i++)
	{
		//loop through columns
		for(int j=1;j<m;j++)
		{

			//check if instrument included
			if (theStrikes(i,j).NumericValue() != 0)
			{
				theQuote.Expiry = theSurface(i,0).NumericValue();
				theQuote.Tenor = theSurface(0,j).NumericValue();
				theQuote.Frequency = theFreq(i,j).NumericValue();
				theQuote.SwapRate = theStrikes(i,j).NumericValue();
				theQuote.VolATM = theSurface(i,j).NumericValue();

				quoteSwap.push_back(theQuote);
			}
		}
	}
	myObj->calibrator(quoteSwap,type);       
	vector<double> theTimeSigmas = myObj->getTimeSigmas();
	vector<double> theSigmas = myObj->getSigmas();
	int k = theSigmas.size();
	CellMatrix allParameters(k,3);
	allParameters(0,0) = myObj->getKappa();
	for(int i=0; i<k; i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	return allParameters;
}

CellMatrix dvlCalibrateHWwBootstrap(const string &theName, CellMatrix theStrikes, CellMatrix theFreq, CellMatrix theSurface, std::string &type)
{
	HW* myObj = getObject<HW>(theName);

	int n = theSurface.RowsInStructure();
	int m = theSurface.ColumnsInStructure();

	defSwap theQuote;
	vector<defSwap> quoteSwap;

	//throw if matrices not the same

	//loop through rows
	for(int i=1;i<n;i++)
	{
		//loop through columns
		for(int j=1;j<m;j++)
		{
			//check if instrument included
			if (theStrikes(i,j).NumericValue() != 0)
			{
				theQuote.Expiry = theSurface(i,0).NumericValue();
				theQuote.Tenor = theSurface(0,j).NumericValue();
				theQuote.Frequency = theFreq(i,j).NumericValue();
				theQuote.SwapRate = theStrikes(i,j).NumericValue();
				theQuote.VolATM = theSurface(i,j).NumericValue();
				quoteSwap.push_back(theQuote);
			}
		}
	}
	myObj->calibratorBstrp(quoteSwap,type);       
	vector<double> theTimeSigmas = myObj->getTimeSigmas();
	vector<double> theSigmas = myObj->getSigmas();
	int k = theSigmas.size();
	CellMatrix allParameters(k,3);
	allParameters(0,0) = myObj->getKappa();
	for(int i=0; i<k; i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	return allParameters;
}

CellMatrix dvlSimulatingHW(const string &theName, CellMatrix theTimes, int NP)
{
	HW* myObj = getObject<HW>(theName);
	int n = theTimes.RowsInStructure();
	vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	CellMatrix thePaths(n,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> paths(myObj->simulationHW(times));
		for(int i=0; i<n; i++)
			thePaths(i,j) = paths[i];
	}
	return thePaths;
}

string dvlCreateHWPDEObj(const string &theObjName, double R0, double kappa,       CellMatrix& theTimesSigmas, CellMatrix& theTimesThetas)
{
	int m = theTimesSigmas.RowsInStructure();
	vector<double> timeSigmas, sigmas;
	for(int i=0; i<m; i++)
		if (theTimesSigmas(i,0).IsANumber()&&theTimesSigmas(i,0).NumericValue()!=0)
		{
			timeSigmas.push_back( theTimesSigmas(i,0).NumericValue());
			sigmas.push_back(theTimesSigmas(i,1).NumericValue());
		}
		m = theTimesThetas.RowsInStructure();
		vector<double> timeThetas, thetas;
		for(int i=0; i<m; i++)
			if (theTimesThetas(i,0).IsANumber()&&theTimesThetas(i,0).NumericValue()!=0)
			{
				timeThetas.push_back( theTimesThetas(i,0).NumericValue());
				thetas.push_back( theTimesThetas(i,1).NumericValue());
			}
			return placeObject(theObjName,new HWPDE(R0, kappa, timeSigmas, sigmas, timeThetas, thetas));
};

string dvlCreateHWtoHWPDEObj(const string &theHWobjName,const string &theHWPDEobjName)
{
	HW* myHW = getObject<HW>(theHWobjName);
	return placeObject(theHWPDEobjName, 
		new HWPDE(myHW->getKappa(), myHW->getTimeSigmas(), myHW->getSigmas(), myHW->getDFsTimes(), myHW->getDFs()));
};

CellMatrix dvlHWPDEObjectContent(const string &theName)
{
	HWPDE* theOBJ = getObject<HWPDE>(theName);

	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
}

CellMatrix dvlCalibrateHWPDESurf(const string &theName, CellMatrix theTimesDFs,CellMatrix theStrikes, CellMatrix theFreq, CellMatrix theSurface, string &type)
{

	HWPDE* theOBJ = getObject<HWPDE>(theName);
	int n = theSurface.RowsInStructure();
	int m = theSurface.ColumnsInStructure();
	int l = theTimesDFs.RowsInStructure();
	vector<double> times, dfs;


	defSwap theQuote;
	vector<defSwap> quoteSwap;

	for(int i=1;i<n;i++)
	{
		//loop through columns
		for(int j=1;j<m;j++)
		{

			//check if instrument included
			if (theStrikes(i,j).NumericValue() != 0)
			{
				theQuote.Expiry = theSurface(i,0).NumericValue();
				theQuote.Tenor = theSurface(0,j).NumericValue();
				theQuote.Frequency = theFreq(i,j).NumericValue();
				theQuote.SwapRate = theStrikes(i,j).NumericValue();
				theQuote.VolATM = theSurface(i,j).NumericValue();

				quoteSwap.push_back(theQuote);
			}
		}
	}
	if (quoteSwap.size() == 0) throw("Error. Wrong inputs");
	double Max_Maturity = quoteSwap[quoteSwap.size()-1].Expiry+quoteSwap[quoteSwap.size()-1].Tenor;
	for(int i=0; i<l; i++)
	{
		times.push_back(theTimesDFs(i,0).NumericValue());
		dfs.push_back(theTimesDFs(i,1).NumericValue());
		if (dfs[i]>=Max_Maturity) break;
	}

	if (type == "Bootstrap") theOBJ->calibratorBootStrap(times,dfs,quoteSwap); 
	else theOBJ->calibrator(times,dfs,quoteSwap);       
	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
}

CellMatrix dvlCalibrateHWPDE(const string &theName, CellMatrix theTimesDFs, CellMatrix theSwaptions)
{
	HWPDE* theOBJ = getObject<HWPDE>(theName);
	int n = theTimesDFs.RowsInStructure();
	vector<double> times(n), dfs(n);
	for(int i=0; i<n; i++)
	{
		times[i] = theTimesDFs(i,0).NumericValue();
		dfs[i] = theTimesDFs(i,1).NumericValue();
	}
	int m = theSwaptions.RowsInStructure();
	vector<defSwap> swaptions(m);
	for(int i=0; i<m; i++)
	{
		swaptions[i].Expiry = theSwaptions(i,0).NumericValue();
		swaptions[i].Tenor = theSwaptions(i,1).NumericValue();
		swaptions[i].Frequency = theSwaptions(i,2).NumericValue();
		swaptions[i].SwapRate = theSwaptions(i,3).NumericValue();
		swaptions[i].VolATM = theSwaptions(i,4).NumericValue();
	}

	theOBJ->calibrator(times, dfs, swaptions);
	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
};

CellMatrix dvlCalibrateHWPDEwBtstrp(const string &theName, CellMatrix theTimesDFs, CellMatrix theSwaptions)
{
	HWPDE* theOBJ = getObject<HWPDE>(theName);
	int n = theTimesDFs.RowsInStructure();
	vector<double> times(n), dfs(n);
	for(int i=0; i<n; i++)
	{
		times[i] = theTimesDFs(i,0).NumericValue();
		dfs[i] = theTimesDFs(i,1).NumericValue();
	}
	int m = theSwaptions.RowsInStructure();
	vector<defSwap> swaptions(m);
	for(int i=0; i<m; i++)
	{
		swaptions[i].Expiry = theSwaptions(i,0).NumericValue();
		swaptions[i].Tenor = theSwaptions(i,1).NumericValue();
		swaptions[i].Frequency = theSwaptions(i,2).NumericValue();
		swaptions[i].SwapRate = theSwaptions(i,3).NumericValue();
		swaptions[i].VolATM = theSwaptions(i,4).NumericValue();
	}
	theOBJ->calibratorBootStrap(times, dfs, swaptions);
	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
};

//!!
string dvlCalibrateHWPDEObj(const string &theObjName, CellMatrix theTimesDFs, CellMatrix theSwaptions)
{
	int n = theTimesDFs.RowsInStructure();
	vector<double> times(n), dfs(n);
	for(int i=0; i<n; i++)
	{
		times[i] = theTimesDFs(i,0).NumericValue();
		dfs[i] = theTimesDFs(i,1).NumericValue();
	}
	int m = theSwaptions.RowsInStructure();
	vector<defSwap> swaptions(m);
	for(int i=0; i<m; i++)
	{
		swaptions[i].Expiry = theSwaptions(i,0).NumericValue();
		swaptions[i].Tenor = theSwaptions(i,1).NumericValue();
		swaptions[i].Frequency = theSwaptions(i,2).NumericValue();
		swaptions[i].SwapRate = theSwaptions(i,3).NumericValue();
		swaptions[i].VolATM = theSwaptions(i,4).NumericValue();
	}
	return placeObject(theObjName,new HWPDE(times, dfs, swaptions));
};
double dvlPricingZBHWPDE(const string &theName, double maturity)
{
	return  getObject<HWPDE>(theName)->pricingZB(maturity);
};
double dvlPricingSwapHWPDE(const string &theName, double Expiry, double Tenor, double Strike, double PayFrequency)
{
	return  getObject<HWPDE>(theName)->pricingSwap(Expiry, Tenor, Strike, PayFrequency);     
};
double dvlPricingSwaptionHWPDE(const string &theName, double Expiry, double Tenor, double Strike, double PayFrequency)
{
	return  getObject<HWPDE>(theName)->pricingSwaption(Expiry, Tenor, Strike, PayFrequency);
};
CellMatrix dvlPricingDFsHWPDE(const string &theName, CellMatrix& theTimes)
{
	int m = theTimes.RowsInStructure();
	vector<double> Times(m);
	for(int i=0; i<m; i++)
		Times[i] = theTimes(i,0).NumericValue();
	vector<double>  DFs = getObject<HWPDE>(theName)->getDFs(Times);
	CellMatrix theDFs(m,1);
	for(int i=0; i<m; i++)
		theDFs(i,0) = DFs[i];
	return theDFs;
};
double dvlPricingBermudanHWPDE(const string &theName, double Expiry, double Tenor, CellMatrix& theExercises, double Strike, double PayFrequency)
{
	int m = theExercises.RowsInStructure();
	vector<double> Exercises(m);
	for(int i=0; i<m; i++)	Exercises[i] = theExercises(i,0).NumericValue();
	return getObject<HWPDE>(theName)->pricingBermudan(Expiry, Tenor, Exercises, Strike, PayFrequency);
};
double dvlGetImpVolATMHWPDE(const string &theName, double Expiry, double Tenor, double PayFrequency)
{
	return getObject<HWPDE>(theName)->getImpVolATM(Expiry, Tenor, PayFrequency);      
};
double dvlGetSwapRateHWPDE(const string &theName, double Expiry, double Tenor, double PayFrequency)
{
	return getObject<HWPDE>(theName)->getSwapRate(Expiry, Tenor, PayFrequency);      
};
CellMatrix dvlSimulatingHWPDE(const string &theName, CellMatrix theTimes, int NP)
{
	HWPDE* theOBJ = getObject<HWPDE>(theName);
	int n = theTimes.RowsInStructure();
	vector<double> times(n);
	for(int i=0; i<n; i++)
		times[i] = theTimes(i,0).NumericValue();
	CellMatrix thePaths(n,NP);
	for(int j=0; j<NP; j++)
	{
		vector<double> paths(theOBJ->simulationHWPDE(times));
		for(int i=0; i<n; i++)
			thePaths(i,j) = paths[i];
	}
	return thePaths;
};

string dvlCreateShortRate1FPDEObj(const string &theObjName, double R0, double kappa, double alpha, double beta, double gamma,
								  CellMatrix& theTimesSigmas, CellMatrix& theTimesThetas)
{
	int m = theTimesSigmas.RowsInStructure();
	vector<double> timeSigmas(m), sigmas(m);
	for(int i=0; i<m; i++)
	{
		timeSigmas[i] = theTimesSigmas(i,0).NumericValue();
		sigmas[i] = theTimesSigmas(i,1).NumericValue();
	}
	m = theTimesThetas.RowsInStructure();
	vector<double> timeThetas(m), thetas(m);
	for(int i=0; i<m; i++)
	{
		timeThetas[i] = theTimesThetas(i,0).NumericValue();
		thetas[i] = theTimesThetas(i,1).NumericValue();
	}
	return placeObject(theObjName,new ShortRate1FPDE(R0,kappa,alpha,beta,gamma,timeSigmas,sigmas,timeThetas,thetas));
}



CellMatrix dvlCalibrateShortRate1FPDE(const string &theObjName, CellMatrix theTimesDFs, CellMatrix theSwaps)
{
	ShortRate1FPDE* theOBJ = getObject<ShortRate1FPDE>(theObjName);
	int n = theTimesDFs.RowsInStructure();
	vector<double> times(n), dfs(n);
	for(int i=0; i<n; i++)
	{
		times[i] = theTimesDFs(i,0).NumericValue();
		dfs[i] = theTimesDFs(i,1).NumericValue();
	}
	int m = theSwaps.RowsInStructure();
	vector<defSwap> quoteSwap(m);
	for(int i=0; i<m; i++)
	{
		quoteSwap[i].Expiry = theSwaps(i,0).NumericValue();
		quoteSwap[i].Tenor = theSwaps(i,1).NumericValue();
		quoteSwap[i].Frequency = theSwaps(i,2).NumericValue();
		quoteSwap[i].SwapRate = theSwaps(i,3).NumericValue();
		quoteSwap[i].VolATM = theSwaps(i,4).NumericValue();
		quoteSwap[i].Value = theSwaps(i,5).NumericValue();
	}

	theOBJ->calibrator(times, dfs, quoteSwap);
	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	allParameters(2,0) = theOBJ->getAlpha();
	allParameters(3,0) = theOBJ->getBeta();
	allParameters(4,0) = theOBJ->getGamma();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
};

CellMatrix dvlCalibrateShortRate1FPDESurf(const string &theName, 
										  CellMatrix theTimesDFs,
										  CellMatrix theStrikes, 
										  CellMatrix theFreq, 
										  CellMatrix theSurface)
{
	ShortRate1FPDE* theOBJ = getObject<ShortRate1FPDE>(theName);
	int n = theSurface.RowsInStructure();
	int m = theSurface.ColumnsInStructure();
	int l = theTimesDFs.RowsInStructure();
	vector<double> times, dfs;
	defSwap theQuote;
	vector<defSwap> quoteSwap;

	for(int i=1;i<n;i++)
	{
		//loop through columns
		for(int j=1;j<m;j++)
		{

			//check if instrument included
			if (theStrikes(i,j).NumericValue() != 0)
			{
				theQuote.Expiry = theSurface(i,0).NumericValue();
				theQuote.Tenor = theSurface(0,j).NumericValue();
				theQuote.Frequency = theFreq(i,j).NumericValue();
				theQuote.SwapRate = theStrikes(i,j).NumericValue();
				theQuote.VolATM = theSurface(i,j).NumericValue();

				quoteSwap.push_back(theQuote);
			}
		}
	}
	if (quoteSwap.size() == 0) throw("Error. Wrong inputs");
	double Max_Maturity = quoteSwap[quoteSwap.size()-1].Expiry+quoteSwap[quoteSwap.size()-1].Tenor;
	for(int i=0; i<l; i++)
	{
		times.push_back(theTimesDFs(i,0).NumericValue());
		dfs.push_back(theTimesDFs(i,1).NumericValue());
		if (dfs[i]>=Max_Maturity) break;
	}

	theOBJ->calibrator(times, dfs, quoteSwap);
	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	allParameters(2,0) = theOBJ->getAlpha();
	allParameters(3,0) = theOBJ->getBeta();
	allParameters(4,0) = theOBJ->getGamma();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
};


double dvlPricingZBShortRate1FPDE(const string &theName, double maturity)
{
	return getObject<ShortRate1FPDE>(theName)->pricingZB(maturity);
};
double dvlPricingSwapShortRate1FPDE(const string &theName, double Expiry, double Tenor, double Coupon, double PayFrequency)
{
	return getObject<ShortRate1FPDE>(theName)->pricingSwap(Expiry, Tenor, Coupon, PayFrequency);
};
double dvlGetSwapRate1FPDE(const string &theName, double Expiry, double Tenor, double PayFrequency)
{
	return getObject<ShortRate1FPDE>(theName)->getSwapRate(Expiry, Tenor, PayFrequency);
};


double dvlPricingSwaptionShortRate1FPDE(const string &theName, double Expiry, double Tenor, double Strike, double PayFrequency)
{
	return getObject<ShortRate1FPDE>(theName)->pricingSwaption(Expiry, Tenor, Strike, PayFrequency);
};
CellMatrix dvlPricingDFShortRate1FPDE(const string &theName, CellMatrix& theTimes)
{

	ShortRate1FPDE *myObj = getObject<ShortRate1FPDE>(theName);
	int m = theTimes.RowsInStructure();
	vector<double> Times(m);
	for(int i=0; i<m; i++)
		Times[i] = theTimes(i,0).NumericValue();
	vector<double> DFs = myObj->calculateDFs(Times);
	CellMatrix theDFs(m,1);
	for(int i=0; i<m; i++)
		theDFs(i,0) = DFs[i];
	return theDFs;
};
string dvlCalibrateShortRate1FPDEObj(const string &theObjName, CellMatrix theTimesDFs, CellMatrix theSwaptions)
{
	int n = theTimesDFs.RowsInStructure();
	vector<double> times(n), dfs(n);
	for(int i=0; i<n; i++)
	{
		times[i] = theTimesDFs(i,0).NumericValue();
		dfs[i] = theTimesDFs(i,1).NumericValue();
	}
	int m = theSwaptions.RowsInStructure();
	vector<defSwap> swaptions(m);
	for(int i=0; i<m; i++)
	{
		swaptions[i].Expiry = theSwaptions(i,0).NumericValue();
		swaptions[i].Tenor = theSwaptions(i,1).NumericValue();
		swaptions[i].Frequency = theSwaptions(i,2).NumericValue();
		swaptions[i].SwapRate = theSwaptions(i,3).NumericValue();
		swaptions[i].VolATM = theSwaptions(i,4).NumericValue();
		swaptions[i].Value = theSwaptions(i,5).NumericValue();
	}
	return placeObject(theObjName,new ShortRate1FPDE(times, dfs, swaptions));
};
CellMatrix dvlGetParametersShortRate1FPDE(const string &theName)
{
	ShortRate1FPDE* theOBJ = getObject<ShortRate1FPDE>(theName);
	vector<double> theTimeSigmas = theOBJ->getTimeSigmas();
	vector<double> theSigmas = theOBJ->getSigmas();
	vector<double> theTimeThetas = theOBJ->getTimeThetas();
	vector<double> theThetas = theOBJ->getThetas();
	CellMatrix allParameters(max(theSigmas.size(), theThetas.size()), 5);
	allParameters(0,0) = theOBJ->getR0();
	allParameters(1,0) = theOBJ->getKappa();
	allParameters(2,0) = theOBJ->getAlpha();
	allParameters(3,0) = theOBJ->getBeta();
	allParameters(4,0) = theOBJ->getGamma();
	for(int i=0; i<theSigmas.size(); i++)
	{
		allParameters(i,1) = theTimeSigmas[i];
		allParameters(i,2) = theSigmas[i];
	}
	for(int i=0; i<theThetas.size(); i++)
	{
		allParameters(i,3) = theTimeThetas[i];
		allParameters(i,4) = theThetas[i];
	}
	return allParameters;
};

string dvlCreateShortRate2FPDEObj(const string &theObjName, double kappa1, double kappa2, double lambda,
								  CellMatrix theTimesSigma1s, CellMatrix theTimesSigma2s,       CellMatrix theTimesAlphas)
{
	int m = theTimesSigma1s.RowsInStructure();
	vector<double> timeSigma1s(m), sigma1s(m);
	for(int i=0; i<m; i++)
	{
		timeSigma1s[i] = theTimesSigma1s(i,0).NumericValue();
		sigma1s[i] = theTimesSigma1s(i,1).NumericValue();
	}
	m = theTimesSigma2s.RowsInStructure();
	vector<double> timeSigma2s(m), sigma2s(m);
	for(int i=0; i<m; i++)
	{
		timeSigma2s[i] = theTimesSigma2s(i,0).NumericValue();
		sigma2s[i] = theTimesSigma2s(i,1).NumericValue();
	}
	m = theTimesAlphas.RowsInStructure();
	vector<double> timeAlphas(m), alphas(m);
	for(int i=0; i<m; i++)
	{
		timeAlphas[i] = theTimesAlphas(i,0).NumericValue();
		alphas[i] = theTimesAlphas(i,1).NumericValue();
	}
	return placeObject(theObjName,new ShortRate2FPDE(kappa1,kappa2,lambda,timeSigma1s,sigma1s,timeSigma2s,sigma2s,timeAlphas,alphas)); 
};
CellMatrix dvlCalibrateShortRate2FPDE(const string &theName, CellMatrix theTimesDFs, CellMatrix theSwaps)
{
	int n = theTimesDFs.RowsInStructure();
	vector<double> times(n), dfs(n);
	for(int i=0; i<n; i++)
	{
		times[i] = theTimesDFs(i,0).NumericValue();
		dfs[i] = theTimesDFs(i,1).NumericValue();
	}
	int m = theSwaps.RowsInStructure();
	vector<defSwap> quoteSwap(m);
	for(int i=0; i<m; i++)
	{
		quoteSwap[i].Expiry = theSwaps(i,0).NumericValue();
		quoteSwap[i].Tenor = theSwaps(i,1).NumericValue();
		quoteSwap[i].Frequency = theSwaps(i,2).NumericValue();
		quoteSwap[i].SwapRate = theSwaps(i,3).NumericValue();
		quoteSwap[i].VolATM = theSwaps(i,4).NumericValue();
		quoteSwap[i].Value = theSwaps(i,5).NumericValue();
	}
	ShortRate2FPDE* theOBJ =getObject<ShortRate2FPDE>(theName);
	theOBJ->calibrator(times, dfs, quoteSwap);
	vector<double> theTimeSigma1s = theOBJ->getTimeSigma1s();
	vector<double> theSigma1s = theOBJ->getSigma1s();
	vector<double> theTimeSigma2s = theOBJ->getTimeSigma2s();
	vector<double> theSigma2s = theOBJ->getSigma2s();
	vector<double> theTimeAlphas = theOBJ->getTimeAlphas();
	vector<double> theAlphas = theOBJ->getAlphas();
	n = max(theSigma1s.size(), theSigma2s.size());
	m = (theAlphas.size()>3) ? theAlphas.size() : 3;
	int np = max(n, m);
	CellMatrix allParameters(np, 7);
	allParameters(0,0) = theOBJ->getKappa1();
	allParameters(1,0) = theOBJ->getKappa2();
	allParameters(2,0) = theOBJ->getLambda();
	for(int i=0; i<theSigma1s.size(); i++)
	{
		allParameters(i,1) = theTimeSigma1s[i];
		allParameters(i,2) = theSigma1s[i];
	}
	for(int i=0; i<theSigma2s.size(); i++)
	{
		allParameters(i,3) = theTimeSigma2s[i];
		allParameters(i,4) = theSigma2s[i];
	}
	for(int i=0; i<theAlphas.size(); i++)
	{
		allParameters(i,5) = theTimeAlphas[i];
		allParameters(i,6) = theAlphas[i];
	}
	return allParameters;
};
double dvlPricingZBShortRate2FPDE(const string &theName, double maturity)
{
	return getObject<ShortRate2FPDE>(theName)->pricingZB(maturity);
};
double dvlPricingSwaptionShortRate2FPDE(const string &theName, double Expiry, double Tenor, double Strike, double PayFrequency)
{
	return getObject<ShortRate2FPDE>(theName)->pricingSwaption(Expiry, Tenor, Strike, PayFrequency);
};
CellMatrix dvlPricingDFShortRate2FPDE(const string &theName, CellMatrix& theTimes)
{
	ShortRate2FPDE* theOBJ =getObject<ShortRate2FPDE>(theName);
	int m = theTimes.RowsInStructure();
	vector<double> Times(m);
	for(int i=0; i<m; i++)
		Times[i] = theTimes(i,0).NumericValue();
	vector<double> DFs = theOBJ->calculateDFs(Times);
	CellMatrix theDFs(m,1);
	for(int i=0; i<m; i++)
		theDFs(i,0) = DFs[i];
	return theDFs;
};
double dvlpricingZB2F(double kappa1, double kappa2, double lambda,
					  CellMatrix theTimes, CellMatrix theSigma1s, CellMatrix theSigma2s, CellMatrix theAlphas, 
					  double maturity)
{
	int m = theTimes.RowsInStructure();
	std::vector<double> times(m);
	std::vector<double> sigma1s(m);
	std::vector<double> sigma2s(m);
	std::vector<double> alphas(m);
	for(int i=0; i<m; i++)
	{
		times[i] = theTimes(i,0).NumericValue();
		sigma1s[i] = theSigma1s(i,0).NumericValue();
		sigma2s[i] = theSigma2s(i,0).NumericValue();
		alphas[i] = theAlphas(i,0).NumericValue();
	}
	ShortRate2FPDE * modelF2 = new ShortRate2FPDE(kappa1,kappa2,lambda,times,sigma1s,times,sigma2s,times,alphas);
	double pvZB = modelF2->pricingZB(maturity);
	delete modelF2;
	return pvZB;
};


double dvlBlackATMStrikePlain(const string &theName, double T0, double TN, double PayFrequency)
{
	return getObject<HW>(theName)->getSwapRate(T0, TN, PayFrequency);	
};

double dvlSwaptionPriceToBlackVol(const string &theName, double T0, double TN, double SwapPrice)
{
	return getObject<HW>(theName)->swaptionIVblackPub(T0, TN, SwapPrice);	
};

double dvlSwaptionPriceBlackATM(const string &theName, double T0, double TN, double VolATM)
{
	return getObject<HW>(theName)->get_swaptionATM(T0, TN, VolATM);	
};

double dvlHWPDEcalTime(const string &theName)
{
	return getObject<HWPDE>(theName)->get_cal_time();	
};


CellMatrix dvlSwaptionShortRate1FPDEDiagnostics(const string &theName, double Expiry, double Tenor, double Strike, double PayFrequency)
{
	ShortRate1FPDE* theOBJ =getObject<ShortRate1FPDE>(theName);
	vector<double> surf = theOBJ->SwaptionDiagnostic(Expiry, Tenor, Strike, PayFrequency);
	vector<double> grid = theOBJ->getGrid();
	CellMatrix ret(surf.size(),2);
	for (int i = 0; i< surf.size(); i++)
	{
		ret(i,0) = grid[i];
		ret(i,1) = surf[i];
	}
	return ret;
};
double dvlSwaptionShortRate1FPDEChangeGridFactor(const string &theName,double Rmax, double factor)
{
	getObject<ShortRate1FPDE>(theName)->change_grid_factor(Rmax,factor);
	return factor;
};

CellMatrix dvlSwaptionShortRate1FPDERNDensity(const string &theName, double T1, double T2)
{
	ShortRate1FPDE* theOBJ =getObject<ShortRate1FPDE>(theName);
	vector<double> surf = theOBJ->RiskNeutralDensity(T1,T2);
	vector<double> grid = theOBJ->getGrid();
	CellMatrix ret(surf.size(),2);
	for (int i = 0; i< surf.size(); i++)
	{
		ret(i,0) = grid[i];
		ret(i,1) = surf[i];
	}
	return ret;
};

CellMatrix dvlSABRgetRND(const string &theName)
{
	sabr *myObj = getObject<sabr>(theName);
	vector<double> surf = myObj->getTableQUTL2();
	vector<double> grid = myObj->getTableQUTL1();
	vector<double> cdfs = myObj->getTableQUTL3();
	CellMatrix ret(surf.size(),3);
	for (int i = 0; i< surf.size(); i++)
	{
		ret(i,0) = grid[i];
		ret(i,1) = surf[i];
		ret(i,2) = cdfs[i];
	}
	return ret;
};

double dvlPricingZBHW(const string &theName, double maturity)
{
	return getObject<HW>(theName)->ZC(maturity);
};

double dvlPricingZBOHWPDE(const string &theName, double expiry, double maturity, double strike, const string &type)
{
	return getObject<HWPDE>(theName)->pricingZBO(expiry,maturity,strike, type);
};

double dvlPricingCouponBondHWPDE(const string &theName, double expiry, double maturity, double coupon, double frequency)
{
	return getObject<HWPDE>(theName)->pricingCouponBond(expiry,maturity,coupon,frequency);
};

double dvlPricingCBOHWPDE(const string &theName, double expiry, double maturity, double coupon, double strike, double frequency, const string &type)
{
	return getObject<HWPDE>(theName)->pricingCBO(expiry,maturity,coupon,strike,frequency,type);
};


double dvlPricingCallableSwapHWPDE(const string &theName, double expiry, double maturity, CellMatrix& theExercises, double coupon, double strike, double frequency, const string &type)
{
	int m = theExercises.RowsInStructure();
	vector<double> Exercises(m);
	for(int i=0; i<m; i++)
		Exercises[i] = theExercises(i,0).NumericValue();
	return getObject<HWPDE>(theName)->pricingCallableSwap(expiry,maturity,Exercises,coupon,strike,frequency,type);
};

double dvlPricingBermudanShortRate1FPDE(const string &theName, double Expiry, double Tenor, CellMatrix& theExercises, double Strike, double PayFrequency)
{
	int m = theExercises.RowsInStructure();
	vector<double> Exercises(m);
	for(int i=0; i<m; i++)
		Exercises[i] = theExercises(i,0).NumericValue();
	return getObject<ShortRate1FPDE>(theName)->pricingBermudan(Expiry, Tenor, Exercises, Strike, PayFrequency);
};


double dvlPricingZBOShortRate1FPDE(const string &theName, double expiry, double maturity, double strike, const string &type)
{
	return getObject<ShortRate1FPDE>(theName)->pricingZBO(expiry,maturity,strike, type);
};

double dvlPricingCouponBondShortRate1FPDE(const string &theName, double expiry, double maturity, double coupon, double frequency)
{
	return getObject<ShortRate1FPDE>(theName)->pricingCouponBond(expiry,maturity,coupon,frequency);
};

double dvlPricingCBOShortRate1FPDE(const string &theName, double expiry, double maturity, double coupon, double strike, double frequency, const string &type)
{
	return getObject<ShortRate1FPDE>(theName)->pricingCBO(expiry,maturity,coupon,strike,frequency,type);
};


double dvlPricingCallableSwapShortRate1FPDE(const string &theName, double expiry, double maturity, CellMatrix& theExercises, double coupon, double strike, double frequency, const string &type)
{
	int m = theExercises.RowsInStructure();
	vector<double> Exercises(m);
	for(int i=0; i<m; i++)
		Exercises[i] = theExercises(i,0).NumericValue();
	return getObject<ShortRate1FPDE>(theName)->pricingCallableSwap(expiry,maturity,Exercises,coupon,strike,frequency,type);
};

double dvlGetImpVolATMShortRate1FPDE(const string &theName, double Expiry, double Tenor, double PayFrequency)
{
	return getObject<ShortRate1FPDE>(theName)->getImpVolATM(Expiry, Tenor, PayFrequency);      
};

string dvlCreateSABRPDEObjNormal(const string &theObjName, 
								 double alpha, double nu, double rho, 
								 double maturity, double forward, double Fmin, double Fmax)
{
	double eps,shift,beta;
	shift=beta=0.0;
	eps = 0.01;
	int sizeX = 512;
	int sizeT = 5000;
	return placeObject(theObjName, new pdeSABR(alpha,beta,nu,rho,eps,shift,maturity, forward,Fmin,Fmax,sizeX,sizeT));
};
string dvlCreateSABRPDEObjShifted(const string &theObjName, 
								 double alpha, double beta, double nu, double rho, double shift, 
								 double maturity, double forward)
{
	double eps = 1;//1e-3;
	int sizeX = 512;
	int sizeT = 200;
	if ((beta <0)||(beta>1.0)) throw("Beta is outside the bounds [0,1)");
//	if (beta==0.0) throw("For Beta = 0 use dvlCreateSABRPDEObjNormal");
	return placeObject(theObjName, new pdeSABR(alpha,beta,nu,rho,eps,shift,maturity, forward,0,0,sizeX,sizeT));
};




CellMatrix dvlCalibrateSABRPDEObj(const string &theName, CellMatrix strikes, CellMatrix quotes)
{
	
	int m = strikes.ColumnsInStructure();
	vector<volQuote> quotes_(m);
	for (int i = 0; i < m ;i++)
	{
		quotes_[i].Strike = strikes(0,i).NumericValue();
		quotes_[i].IV = quotes(0,i).NumericValue();
	}
	pdeSABR *myObj = getObject<pdeSABR>(theName);
	myObj->calibrator(quotes_);
	CellMatrix out(1,5);
	out(0,0) = myObj->getAlpha();
	out(0,1) = myObj->getBeta();
	out(0,2) = myObj->getNu();
	out(0,3) = myObj->getRho();
	out(0,4) = myObj->getShift();
	return out;
};
CellMatrix dvlGetSABRPDEDensity(const string &theName)
{
	pdeSABR *myObj = getObject<pdeSABR>(theName);	
	vector<double> Q = myObj->getDensity();
	vector<double> F = myObj->getFgrid();
	int N = F.size();
	CellMatrix out(N,2);
	for(int i=0; i<N;i++)
	{
		out(i,0)=F[i];
		out(i,1)=Q[i];
	}
	return out;
};

double dvlPriceSABRPDE(const string &theName, double strike, const string &optionType)
{
	return getObject<pdeSABR>(theName)->sabr_option(strike, optionType);
};

double dvlGetVolSABRPDE(const string &theName, double strike)
{
	return getObject<pdeSABR>(theName)->sabrVol(strike);
};


string dvlCreateNoArbSABR(const string &theObjName, int sizeX, int sizeT, double nd)
{
	sabr * mo = getObject<sabr>(theObjName);
	string name(stripTrailingDot(theObjName));
	if ((mo->getBeta() <0)||(mo->getBeta()>1.0)) throw("Beta is outside the bounds [0,1)");
	return placeObject("NoArb"+name, new afSABR(	mo->getParameterAlpha(),	mo->getBeta(),
												mo-> getParameterNu(),	mo->getParameterRho(),
												mo->getShift(),	mo->getMaturity(),  mo->getForward(),	
												sizeX,	sizeT, nd)	);
};

string dvlCreateAntonovSABR(const string &theObjName, int sizeX, int sizeT, double nd)
{
	string name(stripTrailingDot(theObjName));
	sabr * mo = getObject<sabr>(theObjName);
	if ((mo->getBeta() <0)||(mo->getBeta()>1.0)) throw("Beta is outside the bounds [0,1)");
	return placeObject("Antonov"+name, new antonovSABR(	mo->getParameterAlpha(),	mo->getBeta(),
												mo-> getParameterNu(),	mo->getParameterRho(),
												mo->getMaturity(),  mo->getForward(),	
												sizeX,	sizeT, nd)	);
};

CellMatrix dvlSABR_RND(const string &theName)
{
	sabr_pde *myObj;
	try
	{
		myObj = getObject<afSABR>(theName);	
	}
	catch(string &e)
	{
		try
		{
			myObj = getObject<antonovSABR>(theName);
		}
		catch(string &e)
		{
			throw(e);
		}
	}
	vector<double> Q = myObj->getDensity();
	vector<double> F = myObj->getFgrid();
	int N = F.size();
	CellMatrix out(N,2);
	for(int i=0; i<N;i++)
	{
		out(i,0)=F[i];
		out(i,1)=Q[i];
	}
	return out;
};
}