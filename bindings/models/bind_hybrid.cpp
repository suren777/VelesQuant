#include "../bind_common.h"
#include <velesquant/volatility/hhw.h>
#include <velesquant/volatility/schobzhu.h>
#include <velesquant/volatility/termstructure.h>

namespace velesquant {
namespace bindings {

void bind_hybrid(py::module_ &m) {
  // HHW - Hybrid Hull-White
  py::class_<HHW>(m, "HHW")
      .def(py::init<double, double, double, double, double, double, double,
                    double, double>(),
           py::arg("s0"), py::arg("v0"), py::arg("initial_rate"),
           py::arg("kappa"), py::arg("eta"), py::arg("rho"), py::arg("sigma1"),
           py::arg("sigma2"), py::arg("a"))
      .def("price",
           static_cast<double (HHW::*)(double, double) const>(&HHW::HHWPrice),
           py::arg("maturity"), py::arg("strike"));

  // Schobel-Zhu
  py::class_<SchobelZhu>(m, "SchobelZhu")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("spot"), py::arg("var0"), py::arg("kappa"), py::arg("theta"),
           py::arg("xi"), py::arg("rho"))
      .def("price", &SchobelZhu::SchobelPrice, py::arg("maturity"),
           py::arg("forward"), py::arg("strike"))
      .def("simulate", &SchobelZhu::simulation, py::arg("times"),
           py::arg("forwards"))
      .def("calibrate", &SchobelZhu::calibrator, py::arg("maturities"),
           py::arg("forwards"), py::arg("strikes"), py::arg("market_quotes"),
           py::arg("target") = CalibrationTarget::Price)
      .def_property("var0", &SchobelZhu::getParameterVar0,
                    &SchobelZhu::setParameterVar0)
      .def_property("kappa", &SchobelZhu::getParameterKappa,
                    &SchobelZhu::setParameterKappa)
      .def_property("theta", &SchobelZhu::getParameterTheta,
                    &SchobelZhu::setParameterTheta)
      .def_property("xi", &SchobelZhu::getParameterXi,
                    &SchobelZhu::setParameterXi)
      .def_property("rho", &SchobelZhu::getParameterRho,
                    &SchobelZhu::setParameterRho);

  // Termstructure
  py::class_<Termstructure>(m, "Termstructure")
      .def(py::init<std::vector<double>, std::vector<double>, CalendarType,
                    DayCounterType>(),
           py::arg("days"), py::arg("rate"), py::arg("calendar"),
           py::arg("daycount"))
      .def(py::init<std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>, CalendarType,
                    DayCounterType>(),
           py::arg("days"), py::arg("rate"), py::arg("qdays"),
           py::arg("dividend"), py::arg("calendar"), py::arg("daycount"))
      .def(
          "discount",
          static_cast<double (Termstructure::*)(int)>(&Termstructure::discount),
          py::arg("date"))
      .def("rate",
           static_cast<double (Termstructure::*)(int, int)>(
               &Termstructure::rate),
           py::arg("date"), py::arg("tenor"))
      .def("dividend",
           static_cast<double (Termstructure::*)(int, int)>(
               &Termstructure::divident),
           py::arg("date"), py::arg("tenor"));
}

} // namespace bindings
} // namespace velesquant
