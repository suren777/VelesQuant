#include "../bind_common.h"
#include <velesquant/volatility/l_vol.h>

namespace velesquant {
namespace bindings {

void bind_local_vol(py::module_ &m) {
  py::class_<lVol>(m, "LocalVol")
      .def(py::init<>())
      .def(py::init<std::vector<Sabr>>(), py::arg("sabr_models"))
      .def(py::init<std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>, double>(),
           py::arg("maturities"), py::arg("forwards"), py::arg("betas"),
           py::arg("alphas"), py::arg("nus"), py::arg("rhos"), py::arg("spot"))
      .def("call_pde", &lVol::callPDE, py::arg("maturity"), py::arg("strike"),
           py::arg("n") = 100)
      .def("put_pde", &lVol::putPDE, py::arg("maturity"), py::arg("strike"),
           py::arg("n") = 100)
      .def("dnt_pde", &lVol::dntPDE, py::arg("maturity"),
           py::arg("upper_barrier"), py::arg("lower_barrier"),
           py::arg("n") = 100)
      .def("density", &lVol::density, py::arg("maturity"), py::arg("nt"))
      .def("export_lv", &lVol::exportLV, py::arg("times"))
      .def("simulate",
           static_cast<std::vector<double> (lVol::*)(
               std::vector<double>, std::vector<double>) const>(
               &lVol::simulation),
           py::arg("times"), py::arg("rands"))
      .def(
          "simulate_paths",
          static_cast<std::vector<double> (lVol::*)(std::vector<double>) const>(
              &lVol::simulation),
          py::arg("times"))
      .def_property("spot", &lVol::getSpot, &lVol::setSpot);
}

} // namespace bindings
} // namespace velesquant
