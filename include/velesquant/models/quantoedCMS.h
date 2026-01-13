//		quantoedCMS.h

#ifndef QUANTOEDCMS_H
#define QUANTOEDCMS_H

#include <string>
#include <velesquant/local_vol/sabr.h>
#include <velesquant/xlw/CellMatrix.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
using namespace xlw;

class quantoedCMS
{
public:
	quantoedCMS(double expirySR, double tenorSR, double forwardSR, double annuitySR, 
		double payCMS, double discountCMS, double corFX, double atmVolFX, 
		double beta, CellMatrix& strikes, CellMatrix& marketQuotes, std::string quoteType="premium");

	~quantoedCMS() { delete quantoedCMSSABR_; };
	
	double fairValue(double strike, std::string callORput="call") const;
    double simulation(double corrRN);	
	double getForward() const;

private:
	sabr* quantoedCMSSABR_;

	double quantoAdj(double expirySR, double forwardSR, double atmVol, double corFX, double atmVolFX) const;
	double convAdj(double expirySR, double tenorSR, double forwardSR, double quantoedForwardSR,
		           double annuitySR, double discountCMS, double atmVol) const;
};

}
#endif