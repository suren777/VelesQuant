// bindings/models/bind_heston.cpp - Heston model bindings
#include "../bind_common.h"
#include <velesquant/models/utility.h>
#include <velesquant/volatility/s_vol.h>

namespace velesquant {
namespace bindings {

void bind_heston(py::module_ &m) {
  py::class_<sVol, std::shared_ptr<sVol>>(m, "Heston")
      .def(py::init<>())
      .def(py::init<double, double, double, double, double, double, int>(),
           py::arg("spot"), py::arg("var0"), py::arg("kappa"), py::arg("theta"),
           py::arg("xi"), py::arg("rho"), py::arg("seed") = 42)
      // Methods with snake_case
      .def("price", &sVol::hestonPrice, py::arg("maturity"), py::arg("forward"),
           py::arg("strike"), py::arg("opt_type") = OptionType::Call,
           "Calculate Heston option price")
      .def("simulate", &sVol::simulationHeston, py::arg("times"),
           py::arg("forwards"), "Simulate Heston paths")
      // Properties (read-only)
      .def_property_readonly("var0", &sVol::getParameterVar0)
      .def_property_readonly("kappa", &sVol::getParameterKappa)
      .def_property_readonly("theta", &sVol::getParameterTheta)
      .def_property_readonly("xi", &sVol::getParameterXi)
      .def_property_readonly("rho", &sVol::getParameterRho)
      // __repr__
      .def("__repr__",
           [](const sVol &h) {
             return "<Heston var0=" + std::to_string(h.getParameterVar0()) +
                    " kappa=" + std::to_string(h.getParameterKappa()) +
                    " theta=" + std::to_string(h.getParameterTheta()) + ">";
           })
      // Calibration
      .def(
          "calibrate",
          [](sVol &self, const Eigen::MatrixXd &maturities,
             const Eigen::MatrixXd &forwards, const Eigen::MatrixXd &strikes,
             const Eigen::MatrixXd &quotes, CalibrationTarget quoteType) {
            std::vector<double> m_vec, f_vec, k_vec, q_vec;
            int m = strikes.rows();
            int k = quotes.rows();
            int l = quotes.cols();

            if ((m == k) && (l == 1)) {
              for (int i = 0; i < m; i++) {
                m_vec.push_back(maturities(i, 0));
                f_vec.push_back(forwards(i, 0));
                k_vec.push_back(strikes(i, 0));
                q_vec.push_back(quotes(i, 0));
              }
            } else {
              m = strikes.cols();
              if (m != l)
                throw std::runtime_error("Strikes don't match vols dimensions");

              for (int i = 0; i < k; i++)
                for (int j = 0; j < l; j++)
                  if (quotes(i, j) > 1e-10) {
                    m_vec.push_back(maturities(i, 0));
                    f_vec.push_back(forwards(i, 0));
                    k_vec.push_back(strikes(0, j));
                    q_vec.push_back(quotes(i, j));
                  }
            }

            self.calibrator(m_vec, f_vec, k_vec, q_vec, quoteType);
            return py::dict("var0"_a = self.getParameterVar0(),
                            "kappa"_a = self.getParameterKappa(),
                            "theta"_a = self.getParameterTheta(),
                            "xi"_a = self.getParameterXi(),
                            "rho"_a = self.getParameterRho());
          },
          py::arg("maturities"), py::arg("forwards"), py::arg("strikes"),
          py::arg("quotes"), py::arg("quote_type") = CalibrationTarget::Price,
          "Calibrate Heston model to market quotes");
}

} // namespace bindings
} // namespace velesquant
