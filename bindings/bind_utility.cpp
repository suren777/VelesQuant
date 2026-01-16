#include "bind_common.h"
#include <velesquant/models/black_formula.h>
#include <velesquant/models/black_scholes.h>
#include <velesquant/models/utility.h>

namespace velesquant {
namespace bindings {

void bind_utility(py::module_ &m) {
  m.def("cdf_normal", &cdf_normal,
        "Cumulative distribution function for normal distribution");
  m.def("pdf_normal", &pdf_normal,
        "Probability density function for normal distribution");
  m.def("implied_vol", &implied_vol, py::arg("maturity"), py::arg("forward"),
        py::arg("strike"), py::arg("price"),
        "Calculate implied volatility for a Black-Scholes price");

  // Black-Scholes / Black Formula
  m.def("black_scholes_call", &BlackScholesCall, py::arg("spot"),
        py::arg("strike"), py::arg("rate"), py::arg("d"), py::arg("vol"),
        py::arg("expiry"));
  m.def("black_scholes_call_vega", &BlackScholesCallVega, py::arg("spot"),
        py::arg("strike"), py::arg("rate"), py::arg("d"), py::arg("vol"),
        py::arg("expiry"));
  m.def("black_formula_call", &BlackFormulaCall, py::arg("forward"),
        py::arg("strike"), py::arg("vol"), py::arg("expiry"));
  m.def("black_formula_call_vega", &BlackFormulaCallVega, py::arg("forward"),
        py::arg("strike"), py::arg("vol"), py::arg("expiry"));
}

} // namespace bindings
} // namespace velesquant
