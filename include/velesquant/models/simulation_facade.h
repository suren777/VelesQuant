#ifndef SIMULATION_FACADE_H
#define SIMULATION_FACADE_H

#include <velesquant/types.h>

namespace velesquant {

Matrix simulationMLV(Matrix u1Sabr, double u1Spot, Matrix u2Sabr, double u2Spot,
                     Matrix u3Sabr, double u3Spot, Matrix corrMatrix,
                     Matrix theTimes, int NP);

Matrix simulation2LV(Matrix the1Maturities, Matrix the1Forwards,
                     Matrix the1Betas, Matrix the1Alphas, Matrix the1Nus,
                     Matrix the1Rhos, double the1Spot, Matrix the2Maturities,
                     Matrix the2Forwards, Matrix the2Betas, Matrix the2Alphas,
                     Matrix the2Nus, Matrix the2Rhos, double the2Spot,
                     double theCorr, Matrix theTimes, int NP);

Matrix simoreLV(Matrix theMaturities, Matrix theForwards, Matrix theBetas,
                Matrix theAlphas, Matrix theNus, Matrix theRhos,
                Matrix theTimes, double spot, int NP);
} // namespace velesquant

#endif
