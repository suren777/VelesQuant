// bindings/models/bind_hull_white.cpp - Hull-White and interest rate bindings
#include "../bind_common.h"
#include <velesquant/engines/hullwhite_analytic_engine.h>
#include <velesquant/models/hullwhite_model.h>

namespace velesquant {
namespace bindings {

void bind_hull_white(py::module_ &m) {

  // Bind HullWhiteModel
  py::class_<models::HullWhiteModel, std::shared_ptr<models::HullWhiteModel>>(
      m, "HullWhiteModel")
      .def(py::init<double, const std::vector<double> &,
                    const std::vector<double> &, const std::vector<double> &,
                    const std::vector<double> &>(),
           py::arg("kappa"), py::arg("time_sigmas"), py::arg("sigmas"),
           py::arg("time_dfs"), py::arg("dfs"))
      .def("get_kappa", &models::HullWhiteModel::getKappa)
      .def("get_sigmas", &models::HullWhiteModel::getSigmas)
      .def("get_time_sigmas", &models::HullWhiteModel::getTimeSigmas)
      .def("get_dfs", &models::HullWhiteModel::getDFs)
      .def("get_time_dfs", &models::HullWhiteModel::getTimeDFs)
      .def("get_discount_factor", &models::HullWhiteModel::getDiscountFactor,
           py::arg("t"))
      .def("simulate", &models::HullWhiteModel::simulation, py::arg("times"),
           "Simulate Hull-White paths")
      // Expose A and B for educational/debug purposes if needed
      .def("A", &models::HullWhiteModel::A, py::arg("t"), py::arg("T"))
      .def("B", &models::HullWhiteModel::B, py::arg("t"), py::arg("T"));

  // Bind HullWhiteAnalyticEngine
  using HWEngine = engines::HullWhiteAnalyticEngine<models::HullWhiteModel>;

  py::class_<HWEngine, std::shared_ptr<HWEngine>>(m, "HullWhiteAnalyticEngine")
      .def(py::init<std::shared_ptr<models::HullWhiteModel>>(),
           py::arg("model"))
      .def("option_bond", &HWEngine::optionBond, py::arg("expiry"),
           py::arg("maturity"), py::arg("strike"), py::arg("type"),
           "Price a bond option")
      .def("swaption", &HWEngine::swaption, py::arg("expiry"), py::arg("tenor"),
           py::arg("strike"), py::arg("pay_frequency") = 0.5,
           "Price a swaption");
}

} // namespace bindings
} // namespace velesquant
