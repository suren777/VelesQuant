#ifndef LOCAL_VOL_FACADE_H
#define LOCAL_VOL_FACADE_H

#include <velesquant/types.h>

namespace velesquant {

Matrix theDensity(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                  Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                  double maurity);

Matrix lvExport(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                Matrix theTimes);

double dntPDElv(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                double maurity, double upperBarrier, double lowerBarrier);

double putPDElv(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                double maurity, double strike);

Matrix callPDElv(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                 Matrix theAlphas, Matrix theNus, Matrix theRhos, double spot,
                 double maurity, double strike);

} // namespace velesquant

#endif
