#include "../bind_common.h"
#include <velesquant/volatility/c_tree.h>
#include <velesquant/volatility/skew_mc.h>

namespace velesquant {
namespace bindings {

void bind_trees(py::module_ &m) {
  py::class_<CTree>(m, "CTree")
      .def(py::init<>())
      .def(py::init<double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("S"), py::arg("T"), py::arg("F"), py::arg("IV"),
           py::arg("rate") = std::vector<double>(),
           py::arg("q") = std::vector<double>())
      .def("calculate_binomial",
           static_cast<double (CTree::*)(double, double, int, exStyle,
                                         OptionType, tType)>(
               &CTree::calculateBinomial),
           py::arg("strike"), py::arg("maturity"), py::arg("n_nodes"),
           py::arg("style"), py::arg("pay"), py::arg("tree"))
      .def("calculate_trinomial", &CTree::calculateTrinomial, py::arg("strike"),
           py::arg("maturity"), py::arg("n_nodes"), py::arg("style"),
           py::arg("pay"), py::arg("tree"));

  // skewMC - Skew Monte Carlo
  py::class_<skewMC>(m, "SkewMC")
      .def(py::init<>())
      .def(py::init<std::vector<Sabr>>(), py::arg("sabr_models"))
      .def("simulate", &skewMC::simulation, py::arg("times"), py::arg("spot"),
           py::arg("kappa"));
}

} // namespace bindings
} // namespace velesquant
