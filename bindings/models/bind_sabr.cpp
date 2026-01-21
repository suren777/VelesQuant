// bindings/models/bind_sabr.cpp - SABR model bindings
#include "../bind_common.h"
#include <velesquant/models/utility.h>
#include <velesquant/volatility/sabr.h>

namespace velesquant {
namespace bindings {

void bind_sabr(py::module_ &m) {
  py::class_<Sabr, std::shared_ptr<Sabr>>(m, "Sabr")
      .def(py::init<>())
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("maturity"), py::arg("forward"), py::arg("beta") = 0.85,
           py::arg("alpha") = 0.5, py::arg("nu") = 0.25, py::arg("rho") = -0.75)
      .def(py::init<double, double, double, double, double, double, double>(),
           py::arg("maturity"), py::arg("forward"), py::arg("beta"),
           py::arg("alpha"), py::arg("nu"), py::arg("rho"), py::arg("shift"))
      // Methods with snake_case names
      .def("implied_vol", &Sabr::impliedVol, py::arg("strike"),
           "Calculate implied volatility for a given strike")
      .def("normal_vol", &Sabr::normalVol, py::arg("K"),
           "Calculate normal volatility")
      .def("premium_bachelier", &Sabr::premiumBachelier, py::arg("strike"),
           py::arg("call_or_put") = OptionType::Call,
           "Calculate Bachelier premium")
      .def("premium_black_scholes", &Sabr::premiumBlackScholes,
           py::arg("strike"), py::arg("call_or_put") = OptionType::Call,
           "Calculate Black-Scholes premium")
      .def("local_vol", &Sabr::localVol, py::arg("spot"),
           "Calculate local volatility")
      // Properties
      .def_property("alpha", &Sabr::getParameterAlpha, &Sabr::setParameterAlpha)
      .def_property("nu", &Sabr::getParameterNu, &Sabr::setParameterNu)
      .def_property("rho", &Sabr::getParameterRho, &Sabr::setParameterRho)
      .def_property("maturity", &Sabr::getMaturity, &Sabr::setMaturity)
      .def_property("forward", &Sabr::getForward, &Sabr::setForward)
      .def_property("beta", &Sabr::getBeta, &Sabr::setBeta)
      // __repr__ for debugging
      .def("__repr__",
           [](const Sabr &s) {
             return "<Sabr alpha=" + std::to_string(s.getParameterAlpha()) +
                    " nu=" + std::to_string(s.getParameterNu()) +
                    " rho=" + std::to_string(s.getParameterRho()) + ">";
           })
      // Calibration
      .def(
          "calibrate",
          [](Sabr &self, const Eigen::MatrixXd &strikes,
             const Eigen::MatrixXd &quotes, CalibrationTarget quoteType) {
            std::vector<double> s_vec, q_vec;
            int rows = strikes.rows();
            int cols = strikes.cols();
            for (int i = 0; i < rows; i++) {
              for (int j = 0; j < cols; j++) {
                if (!std::isnan(strikes(i, j)) && !std::isnan(quotes(i, j))) {
                  s_vec.push_back(strikes(i, j));
                  q_vec.push_back(quotes(i, j));
                }
              }
            }
            if (s_vec.empty()) {
              throw std::runtime_error("No valid data points for calibration");
            }
            self.calibrator(s_vec, q_vec, quoteType);
            return py::dict("alpha"_a = self.getParameterAlpha(),
                            "nu"_a = self.getParameterNu(),
                            "rho"_a = self.getParameterRho());
          },
          py::arg("strikes"), py::arg("quotes"),
          py::arg("quote_type") = CalibrationTarget::Volatility,
          "Calibrate SABR model to market quotes");
}

} // namespace bindings
} // namespace velesquant
