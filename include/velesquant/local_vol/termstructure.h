#ifndef TERMSTRUCTURE_H
#define TERMSTRUCTURE_H

#include <ql/quantlib.hpp>

namespace velesquant {


#pragma warning (disable:4244)

typedef std::vector<double> Vdoub;
typedef std::vector<QuantLib::Date> Vdate;

class Termstructure{

public:
	Termstructure(Vdoub days, Vdoub rate, int calendar, int daycount);
	Termstructure(Vdoub days, Vdoub rate, Vdoub qdays, 
		Vdoub divident, int calendar, int daycount);
	~Termstructure();

	double discount(int date); // returns discount factor on a date
	Vdoub discount(Vdate dates); // returns array of days.
	double rate(int date, int tenor); // returns a forward rate on a date with given tenor
	Vdoub rate(Vdate dates, int tenor); //returns forward rates for array of days with given tenors.
	double divident(int date, int tenor); // returns a divident yield on a date with given tenor
	Vdoub divident(Vdate dates, int tenor); //returns divident yields for array of days with given tenors.


private:
	Vdate dates_;
	Vdoub rates_;
	Vdate qdates_;
	Vdoub dividents_;
	QuantLib::Calendar cal_;
	QuantLib::DayCounter mcount_;
	boost::shared_ptr<QuantLib::YieldTermStructure> curve_;
	boost::shared_ptr<QuantLib::YieldTermStructure> qcurve_;

	void input_forwards(Vdoub dates, Vdoub rates)
		{
		int n = dates.size();

		for (int i = 0 ; i < n ; i++ ) 
			{
			dates_.push_back(QuantLib::Date(dates[i]));
			rates_.push_back(rates[i]);			
			}
		boost::shared_ptr<QuantLib::YieldTermStructure> curve (new QuantLib::InterpolatedForwardCurve
			<QuantLib::Linear>(dates_,rates_,mcount_,cal_));
		curve_=curve;
		};
	void input_divident(Vdoub qdates, Vdoub dividents)
		{
		int n = qdates.size();

		for (int i = 0 ; i < n ; i++ ) 
			{
			qdates_.push_back(QuantLib::Date(qdates[i]));
			dividents_.push_back(dividents[i]);			
			}
		boost::shared_ptr<QuantLib::YieldTermStructure> qcurve (new QuantLib::InterpolatedForwardCurve
			<QuantLib::Linear>(qdates_,dividents_,mcount_,cal_));
		qcurve_=qcurve;
		qcurve_->enableExtrapolation();
		};
	void selectCounter(int counter){
		using namespace QuantLib;
		if (counter == 1){
			mcount_=Actual360();
			}
		else if (counter == 2){			
			mcount_=Actual365Fixed();
			}
		else if (counter == 3)	{	
			mcount_=ActualActual(ActualActual::ISDA);
			}
		else if (counter == 4){			
			mcount_=Actual365Fixed();
			}
		else if (counter == 5){	
			mcount_=Business252(cal_);
			}
		else if (counter == 6){		
			mcount_=OneDayCounter();
			}
		else if (counter == 7){				
			mcount_=SimpleDayCounter();
			}
		else if (counter == 8){			
			mcount_=Thirty360(Thirty360::USA);
			}
		else  mcount_=Thirty360(Thirty360::USA);
		}




	void calendarSelector(int calendar)
		{
		using namespace QuantLib;
		if (calendar == 1) cal_=UnitedKingdom();
		else if (calendar == 2) cal_=Argentina();
		else if (calendar == 3) cal_=Australia();
		else if (calendar == 4) cal_=BespokeCalendar();
		else if (calendar == 5) cal_=Brazil();
		else if (calendar == 6) cal_=Canada(); 
		else if (calendar == 7) cal_=China();
		else if (calendar == 8) cal_=CzechRepublic(); 
		else if (calendar == 9) cal_=Denmark();	
		else if (calendar == 10) cal_=Finland();
		else if (calendar == 11) cal_=Germany();
		else if (calendar == 12) cal_=HongKong();
		else if (calendar == 13) cal_=Hungary();
		else if (calendar == 14) cal_=Iceland(); 
		else if (calendar == 15) cal_=India();
		else if (calendar == 16) cal_=Indonesia(); 
		else if (calendar == 17) cal_=Italy();
		else if (calendar == 18) cal_=Japan ();
		else if (calendar == 19) cal_=Ukraine(); 
		else if (calendar == 20) cal_=Mexico();
		else if (calendar == 21) cal_=NewZealand(); 
		else if (calendar == 22) cal_=Norway();
		else if (calendar == 23) cal_=Poland();
		else if (calendar == 24) cal_=Russia();
		else if (calendar == 25) cal_=SaudiArabia();
		else if (calendar == 26) cal_=Singapore();
		else if (calendar == 27) cal_=Slovakia();
		else if (calendar == 28) cal_=SouthAfrica();
		else if (calendar == 29) cal_=SouthKorea(); 
		else if (calendar == 30) cal_=Sweden();
		else if (calendar == 31) cal_=Switzerland(); 
		else if (calendar == 32) cal_=Taiwan();
		else if (calendar == 33) cal_=TARGET(); 
		else if (calendar == 34) cal_=Turkey(); 
		else cal_=UnitedKingdom();
		}


	void kill_all()
		{
		// Removed dangerous and unnecessary deletes of stack members
		}
	};


}
#endif