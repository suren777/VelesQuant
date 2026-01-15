//      lm.h

#ifndef LM_H
#define LM_H

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

/*! \file lm.hpp
    \brief wrapper for MINPACK minimization routine
*/
#include <boost/function.hpp>
#include <ql/qldefines.hpp>

namespace velesquant {

typedef boost::function<void(int, int, double *, double *, int *)> lmfcn;

void lmdif(int m, int n, double *x, double *fvec, double ftol, double xtol,
           double gtol, int maxfev, double epsfcn, double *diag, int mode,
           double factor, int nprint, int *info, int *nfev, double *fjac,
           int ldfjac, int *ipvt, double *qtf, double *wa1, double *wa2,
           double *wa3, double *wa4, const lmfcn &);

void qrsolv(int n, double *r, int ldr, int *ipvt, double *diag, double *qtb,
            double *x, double *sdiag, double *wa);

void qrfac(int m, int n, double *a, int, int pivot, int *ipvt, int,
           double *rdiag, double *acnorm, double *wa);

} // namespace velesquant
#endif