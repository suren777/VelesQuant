// bindings/bind_enums.cpp - Enum bindings
#include "bind_common.h"
#include <velesquant/models/utility.h>

#include <velesquant/volatility/c_tree.h>
#include <velesquant/volatility/termstructure.h>

namespace velesquant {
namespace bindings {

void bind_enums(py::module_ &m) {
  // Tree types
  py::enum_<tType>(m, "TreeType")
      .value("Recombining", tType::recomb)
      .value("NonRecombining", tType::nonrecomb)
      .export_values();

  py::enum_<exStyle>(m, "ExerciseStyle")
      .value("American", exStyle::American)
      .value("European", exStyle::European)
      .value("Bermudan", exStyle::Bermudan)
      .export_values();

  // Calendar types
  py::enum_<CalendarType>(m, "CalendarType")
      .value("UnitedKingdom", CalendarType::CAL_UnitedKingdom)
      .value("Argentina", CalendarType::CAL_Argentina)
      .value("Australia", CalendarType::CAL_Australia)
      .value("BespokeCalendar", CalendarType::CAL_BespokeCalendar)
      .value("Brazil", CalendarType::CAL_Brazil)
      .value("Canada", CalendarType::CAL_Canada)
      .value("China", CalendarType::CAL_China)
      .value("CzechRepublic", CalendarType::CAL_CzechRepublic)
      .value("Denmark", CalendarType::CAL_Denmark)
      .value("Finland", CalendarType::CAL_Finland)
      .value("Germany", CalendarType::CAL_Germany)
      .value("HongKong", CalendarType::CAL_HongKong)
      .value("Hungary", CalendarType::CAL_Hungary)
      .value("Iceland", CalendarType::CAL_Iceland)
      .value("India", CalendarType::CAL_India)
      .value("Indonesia", CalendarType::CAL_Indonesia)
      .value("Italy", CalendarType::CAL_Italy)
      .value("Japan", CalendarType::CAL_Japan)
      .value("Ukraine", CalendarType::CAL_Ukraine)
      .value("Mexico", CalendarType::CAL_Mexico)
      .value("NewZealand", CalendarType::CAL_NewZealand)
      .value("Norway", CalendarType::CAL_Norway)
      .value("Poland", CalendarType::CAL_Poland)
      .value("Russia", CalendarType::CAL_Russia)
      .value("SaudiArabia", CalendarType::CAL_SaudiArabia)
      .value("Singapore", CalendarType::CAL_Singapore)
      .value("Slovakia", CalendarType::CAL_Slovakia)
      .value("SouthAfrica", CalendarType::CAL_SouthAfrica)
      .value("SouthKorea", CalendarType::CAL_SouthKorea)
      .value("Sweden", CalendarType::CAL_Sweden)
      .value("Switzerland", CalendarType::CAL_Switzerland)
      .value("Taiwan", CalendarType::CAL_Taiwan)
      .value("TARGET", CalendarType::CAL_TARGET)
      .value("Turkey", CalendarType::CAL_Turkey)
      .export_values();

  py::enum_<DayCounterType>(m, "DayCounterType")
      .value("Actual360", DayCounterType::DC_Actual360)
      .value("Actual365Fixed", DayCounterType::DC_Actual365Fixed)
      .value("ActualActual", DayCounterType::DC_ActualActual)
      .value("Actual365NoLeap", DayCounterType::DC_Actual365NoLeap)
      .value("Business252", DayCounterType::DC_Business252)
      .value("OneDayCounter", DayCounterType::DC_OneDayCounter)
      .value("SimpleDayCounter", DayCounterType::DC_SimpleDayCounter)
      .value("Thirty360", DayCounterType::DC_Thirty360)
      .export_values();

  py::enum_<OptionType>(m, "OptionType")
      .value("Call", OptionType::Call)
      .value("Put", OptionType::Put)
      .export_values();

  py::enum_<CalibrationTarget>(m, "CalibrationTarget")
      .value("Price", CalibrationTarget::Price)
      .value("Volatility", CalibrationTarget::Volatility)
      .export_values();

  py::class_<defSwap>(m, "DefSwap")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("expiry"), py::arg("tenor"), py::arg("frequency"),
           py::arg("swap_rate"), py::arg("vol_atm"), py::arg("value"))
      .def(py::init<>())
      .def_readwrite("expiry", &defSwap::Expiry)
      .def_readwrite("tenor", &defSwap::Tenor)
      .def_readwrite("frequency", &defSwap::Frequency)
      .def_readwrite("swap_rate", &defSwap::SwapRate)
      .def_readwrite("vol_atm", &defSwap::VolATM)
      .def_readwrite("value", &defSwap::Value);
}

} // namespace bindings
} // namespace velesquant
