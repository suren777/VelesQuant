#include <velesquant/local_vol/termstructure.h>
using namespace std;
#include <algorithm>

namespace velesquant {



//constructor

Termstructure::Termstructure(Vdoub days, Vdoub rate, int calendar, int daycount)
	{
	Vdoub qdays(2);
	Vdoub divident(2,0);
	qdays[0]=days.front();
	qdays[1]=days.back();

	Termstructure(days, rate, qdays, divident, calendar, daycount);
	}


Termstructure::Termstructure(Vdoub days, Vdoub rate, Vdoub qdays,Vdoub divident, int calendar, int daycount)
	{
	calendarSelector(calendar);
	selectCounter(daycount);
	input_forwards(days, rate);
	input_divident(qdays, divident);
	}
//destructor
Termstructure::~Termstructure()
	{
	kill_all(); 
	}

double Termstructure::discount(int date)
	{
	return curve_->discount(date);
	} // returns discount factor on a date

Vdoub Termstructure::discount(Vdate dates)
	{
	int n = dates.size();
	Vdoub discounts(n);
	for(int i = 0; i < n; i++)
		discounts[n] = curve_->discount(dates[i]);
	return discounts;
	} // returns array of days.

double Termstructure::rate(int date, int tenor)
	{
	return curve_->forwardRate(
		QuantLib::Date(date),  
		QuantLib::Date(date+tenor), 
		mcount_, 
		QuantLib::Simple);
	} // returns a forward rate on a date with given tenor

Vdoub Termstructure::rate(Vdate dates, int tenor)
	{
	int n = dates.size();
	Vdoub rates(n);
	for(int i = 0; i < n; i++)
		rates[n] = curve_->forwardRate(
		QuantLib::Date(dates[i]), 
		QuantLib::Date(dates[i]+tenor), 
		mcount_, 
		QuantLib::Simple);
	return rates;
	} //returns forward rates for array of days with given tenors.

double Termstructure::divident(int date, int tenor)
	{
	return qcurve_->forwardRate(
		QuantLib::Date(date),  
		QuantLib::Date(date+tenor), 
		mcount_, 
		QuantLib::Simple);
	} // returns a divident yield on a date with given tenor
Vdoub Termstructure::divident(Vdate dates, int tenor)
	{
	int n = dates.size();
	Vdoub dividents(n);
	for(int i = 0; i < n; i++)
		dividents[n] = qcurve_->forwardRate(
		QuantLib::Date(dates[i]), 
		QuantLib::Date(dates[i]+tenor), 
		mcount_, 
		QuantLib::Simple);
	return dividents;


	} //returns divident yields for array of days with given tenors.
}