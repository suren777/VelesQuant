#ifndef PDE_UTILITY_H
#define PDE_UTILITY_H

#include <boost/function.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions.hpp>
#include <boost/assert.hpp>
#include <boost/current_function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <exception>
#include <sstream>
#include <string>

namespace velesquant {

#define QL_ENSURE(condition,message) \
if (!(condition)) { \
    std::ostringstream _ql_msg_stream; \
    _ql_msg_stream << message; \
    throw(_ql_msg_stream.str()); \
 }

struct defSwap
{
       double Expiry;
       double Tenor; 
       double Frequency; 
       double SwapRate; 
       double VolATM;             
       double Value;       
}; 

#pragma warning (disable:4996)
boost::mt19937 erng; 
boost::normal_distribution<> ndist(0.0, 1.0);
double random_normal()
{
    return ndist(erng);
}

boost::math::normal_distribution<> enormal(0,1);
double cdf_normal(double p)
{
	return cdf(enormal, p);
}

double pdf_normal(double p)
{
	return pdf(enormal, p);
}

double asinh(double a){return log(a+sqrt(a*a+1.0));};



}
#endif