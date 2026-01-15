// bindings/models/bind_hull_white.cpp - Hull-White and interest rate bindings
#include "../bind_common.h"
#include <velesquant/models/hw.h>
#include <velesquant/models/utility.h>

namespace velesquant {
namespace bindings {

void bind_hull_white(py::module_ &m) {
  py::class_<HullWhite, std::shared_ptr<HullWhite>>(m, "HullWhite")
      .def(py::init<double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>>(),
           py::arg("kappa"), py::arg("time_sigmas"), py::arg("sigmas"),
           py::arg("discount_factor_times"), py::arg("discount_factors"))
      // Methods with snake_case
      .def("option_bond", &HullWhite::optionBond, py::arg("expiry"),
           py::arg("maturity"), py::arg("strike"), py::arg("type"),
           "Price a bond option")
      .def("swaption", &HullWhite::swaption, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency") = 0.5,
           "Price a swaption")
      .def("simulate", &HullWhite::simulation, py::arg("times"),
           "Simulate Hull-White paths")
      .def("zero_coupon", &HullWhite::ZC, py::arg("expiry"),
           "Calculate zero-coupon bond price")
      .def("calibrate", &HullWhite::calibrator, py::arg("swap_quotes"),
           py::arg("target") = CalibrationTarget::Volatility,
           "Calibrate model to swaption quotes")
      // Getters
      .def("get_kappa", &HullWhite::getKappa)
      .def("get_sigmas", &HullWhite::getSigmas)
      .def("get_time_sigmas", &HullWhite::getTimeSigmas)
      // Vectorized methods
      .def("swaption_vec",
           py::vectorize([](HullWhite &self, double expiry, double tenor,
                            double strike, double pay_freq) {
             return self.swaption(expiry, tenor, strike, pay_freq);
           }),
           py::arg("expiry"), py::arg("tenor"), py::arg("strike"),
           py::arg("pay_frequency") = 0.5)
      .def("zero_coupon_vec", py::vectorize([](HullWhite &self, double expiry) {
             return self.ZC(expiry);
           }),
           py::arg("expiry"))
      // __repr__
      .def("__repr__",
           [](HullWhite &hw) {
             return "<HullWhite kappa=" + std::to_string(hw.getKappa()) + ">";
           })
      // NumPy array returns for zero-copy performance
      .def(
          "get_sigmas_array",
          [](HullWhite &hw) {
            auto vec = hw.getSigmas();
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::return_value_policy::move, "Get sigmas as NumPy array")
      .def(
          "get_time_sigmas_array",
          [](HullWhite &hw) {
            auto vec = hw.getTimeSigmas();
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::return_value_policy::move, "Get time sigmas as NumPy array")
      .def(
          "simulate_array",
          [](HullWhite &hw, const std::vector<double> &times) {
            auto vec = hw.simulation(times);
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::arg("times"), py::return_value_policy::move,
          "Simulate and return as NumPy array");
}

} // namespace bindings
} // namespace velesquant
