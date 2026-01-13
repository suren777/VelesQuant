#ifndef HESTONPDE_H
#define HESTONPDE_H

#include <vector>
#include <complex>
#include "lapacke.h"

namespace velesquant {

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>


typedef std::vector<double> Vdoub;
typedef std::complex<double> Cdoub;


class HestonPDE
	{
	public:
		struct 
			{
			public CCS_Matrix derS;
			public double[,] derSS;
			public double[,] derV1;
			public double[,] derV2;
			public double[,] derVV;
			public double[,] derSV;
			public double[,] R;
			}


	}


}
#endif 