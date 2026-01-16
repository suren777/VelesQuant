#include "../bind_common.h"
#include <velesquant/models/swaption.h>

namespace velesquant {
namespace bindings {

void bind_swaption(py::module_ &m) {
  py::class_<swaption>(m, "Swaption")
      .def(py::init<double, double, double, double, double, double, double,
                    double>(),
           py::arg("expiry"), py::arg("tenor"), py::arg("forward"),
           py::arg("annuity"), py::arg("beta") = 0.85, py::arg("alpha") = 0.5,
           py::arg("nu") = 0.25, py::arg("rho") = -0.75)
      .def("fair_value", &swaption::swaptionFairValue, py::arg("strike"),
           py::arg("call_or_put") = OptionType::Call)
      .def("swap_fair_value", &swaption::swapFairValue, py::arg("strike"))
      .def("get_implied_vol", &swaption::getImpliedVol, py::arg("strike"))
      .def_property_readonly("alpha", &swaption::getParameterAlpha)
      .def_property_readonly("nu", &swaption::getParameterNu)
      .def_property_readonly("rho", &swaption::getParameterRho);
}

} // namespace bindings
} // namespace velesquant
