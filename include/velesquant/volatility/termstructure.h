#ifndef TERMSTRUCTURE_H
#define TERMSTRUCTURE_H

#include <ql/quantlib.hpp>

namespace velesquant {

#ifdef _MSC_VER
#pragma warning(disable : 4244)
#endif

typedef std::vector<double> Vdoub;
typedef std::vector<QuantLib::Date> Vdate;

// Calendar types for Termstructure
enum class CalendarType {
  CAL_UnitedKingdom = 1,
  CAL_Argentina = 2,
  CAL_Australia = 3,
  CAL_BespokeCalendar = 4,
  CAL_Brazil = 5,
  CAL_Canada = 6,
  CAL_China = 7,
  CAL_CzechRepublic = 8,
  CAL_Denmark = 9,
  CAL_Finland = 10,
  CAL_Germany = 11,
  CAL_HongKong = 12,
  CAL_Hungary = 13,
  CAL_Iceland = 14,
  CAL_India = 15,
  CAL_Indonesia = 16,
  CAL_Italy = 17,
  CAL_Japan = 18,
  CAL_Ukraine = 19,
  CAL_Mexico = 20,
  CAL_NewZealand = 21,
  CAL_Norway = 22,
  CAL_Poland = 23,
  CAL_Russia = 24,
  CAL_SaudiArabia = 25,
  CAL_Singapore = 26,
  CAL_Slovakia = 27,
  CAL_SouthAfrica = 28,
  CAL_SouthKorea = 29,
  CAL_Sweden = 30,
  CAL_Switzerland = 31,
  CAL_Taiwan = 32,
  CAL_TARGET = 33,
  CAL_Turkey = 34
};

// Day counter types for Termstructure
enum class DayCounterType {
  DC_Actual360 = 1,
  DC_Actual365Fixed = 2,
  DC_ActualActual = 3,
  DC_Actual365NoLeap = 4,
  DC_Business252 = 5,
  DC_OneDayCounter = 6,
  DC_SimpleDayCounter = 7,
  DC_Thirty360 = 8
};

class Termstructure {

public:
  Termstructure(Vdoub days, Vdoub rate, CalendarType calendar,
                DayCounterType daycount);
  Termstructure(Vdoub days, Vdoub rate, Vdoub qdays, Vdoub divident,
                CalendarType calendar, DayCounterType daycount);
  ~Termstructure();

  double discount(int date);   // returns discount factor on a date
  Vdoub discount(Vdate dates); // returns array of days.
  double rate(int date,
              int tenor); // returns a forward rate on a date with given tenor
  Vdoub
  rate(Vdate dates,
       int tenor); // returns forward rates for array of days with given tenors.
  double
  divident(int date,
           int tenor); // returns a divident yield on a date with given tenor
  Vdoub divident(Vdate dates, int tenor); // returns divident yields for array
                                          // of days with given tenors.

private:
  Vdate dates_;
  Vdoub rates_;
  Vdate qdates_;
  Vdoub dividents_;
  QuantLib::Calendar cal_;
  QuantLib::DayCounter mcount_;
  boost::shared_ptr<QuantLib::YieldTermStructure> curve_;
  boost::shared_ptr<QuantLib::YieldTermStructure> qcurve_;

  void input_forwards(Vdoub dates, Vdoub rates) {
    int n = dates.size();

    for (int i = 0; i < n; i++) {
      dates_.push_back(QuantLib::Date(dates[i]));
      rates_.push_back(rates[i]);
    }
    boost::shared_ptr<QuantLib::YieldTermStructure> curve(
        new QuantLib::InterpolatedForwardCurve<QuantLib::Linear>(
            dates_, rates_, mcount_, cal_));
    curve_ = curve;
  };
  void input_divident(Vdoub qdates, Vdoub dividents) {
    int n = qdates.size();

    for (int i = 0; i < n; i++) {
      qdates_.push_back(QuantLib::Date(qdates[i]));
      dividents_.push_back(dividents[i]);
    }
    boost::shared_ptr<QuantLib::YieldTermStructure> qcurve(
        new QuantLib::InterpolatedForwardCurve<QuantLib::Linear>(
            qdates_, dividents_, mcount_, cal_));
    qcurve_ = qcurve;
    qcurve_->enableExtrapolation();
  };
  void selectCounter(DayCounterType counter) {
    using namespace QuantLib;
    switch (counter) {
    case DayCounterType::DC_Actual360:
      mcount_ = Actual360();
      break;
    case DayCounterType::DC_Actual365Fixed:
      mcount_ = Actual365Fixed();
      break;
    case DayCounterType::DC_ActualActual:
      mcount_ = ActualActual(ActualActual::ISDA);
      break;
    case DayCounterType::DC_Actual365NoLeap:
      mcount_ = Actual365Fixed();
      break;
    case DayCounterType::DC_Business252:
      mcount_ = Business252(cal_);
      break;
    case DayCounterType::DC_OneDayCounter:
      mcount_ = OneDayCounter();
      break;
    case DayCounterType::DC_SimpleDayCounter:
      mcount_ = SimpleDayCounter();
      break;
    case DayCounterType::DC_Thirty360:
      mcount_ = Thirty360(Thirty360::USA);
      break;
    default:
      mcount_ = Thirty360(Thirty360::USA);
      break;
    }
  }

  void calendarSelector(CalendarType calendar) {
    using namespace QuantLib;
    switch (calendar) {
    case CalendarType::CAL_UnitedKingdom:
      cal_ = UnitedKingdom();
      break;
    case CalendarType::CAL_Argentina:
      cal_ = Argentina();
      break;
    case CalendarType::CAL_Australia:
      cal_ = Australia();
      break;
    case CalendarType::CAL_BespokeCalendar:
      cal_ = BespokeCalendar();
      break;
    case CalendarType::CAL_Brazil:
      cal_ = Brazil();
      break;
    case CalendarType::CAL_Canada:
      cal_ = Canada();
      break;
    case CalendarType::CAL_China:
      cal_ = China();
      break;
    case CalendarType::CAL_CzechRepublic:
      cal_ = CzechRepublic();
      break;
    case CalendarType::CAL_Denmark:
      cal_ = Denmark();
      break;
    case CalendarType::CAL_Finland:
      cal_ = Finland();
      break;
    case CalendarType::CAL_Germany:
      cal_ = Germany();
      break;
    case CalendarType::CAL_HongKong:
      cal_ = HongKong();
      break;
    case CalendarType::CAL_Hungary:
      cal_ = Hungary();
      break;
    case CalendarType::CAL_Iceland:
      cal_ = Iceland();
      break;
    case CalendarType::CAL_India:
      cal_ = India();
      break;
    case CalendarType::CAL_Indonesia:
      cal_ = Indonesia();
      break;
    case CalendarType::CAL_Italy:
      cal_ = Italy();
      break;
    case CalendarType::CAL_Japan:
      cal_ = Japan();
      break;
    case CalendarType::CAL_Ukraine:
      cal_ = Ukraine();
      break;
    case CalendarType::CAL_Mexico:
      cal_ = Mexico();
      break;
    case CalendarType::CAL_NewZealand:
      cal_ = NewZealand();
      break;
    case CalendarType::CAL_Norway:
      cal_ = Norway();
      break;
    case CalendarType::CAL_Poland:
      cal_ = Poland();
      break;
    case CalendarType::CAL_Russia:
      cal_ = Russia();
      break;
    case CalendarType::CAL_SaudiArabia:
      cal_ = SaudiArabia();
      break;
    case CalendarType::CAL_Singapore:
      cal_ = Singapore();
      break;
    case CalendarType::CAL_Slovakia:
      cal_ = Slovakia();
      break;
    case CalendarType::CAL_SouthAfrica:
      cal_ = SouthAfrica();
      break;
    case CalendarType::CAL_SouthKorea:
      cal_ = SouthKorea();
      break;
    case CalendarType::CAL_Sweden:
      cal_ = Sweden();
      break;
    case CalendarType::CAL_Switzerland:
      cal_ = Switzerland();
      break;
    case CalendarType::CAL_Taiwan:
      cal_ = Taiwan();
      break;
    case CalendarType::CAL_TARGET:
      cal_ = TARGET();
      break;
    case CalendarType::CAL_Turkey:
      cal_ = Turkey();
      break;
    default:
      cal_ = UnitedKingdom();
      break;
    }
  }

  void kill_all() {
    // Removed dangerous and unnecessary deletes of stack members
  }
};

} // namespace velesquant
#endif