
#include "bind_common.h"

namespace velesquant {
namespace bindings {

// Declarations from other files
void bind_enums(py::module_ &m);
void bind_utility(py::module_ &m);
void bind_sabr(py::module_ &m);
void bind_heston(py::module_ &m);
void bind_local_vol(py::module_ &m);
void bind_hull_white(py::module_ &m);
void bind_swaption(py::module_ &m);
void bind_cms(py::module_ &m);
void bind_basket(py::module_ &m);
void bind_hybrid(py::module_ &m);
void bind_pde_solvers(py::module_ &m);
void bind_trees(py::module_ &m);
void bind_exceptions(py::module_ &m);

} // namespace bindings
} // namespace velesquant

// Debug export to check if symbols are being stripped
extern "C" __attribute__((visibility("default"))) int
velesquant_debug_symbol() {
  return 42;
}

PYBIND11_MODULE(_core, m) {
  m.doc() = "VelesQuant Core Bindings - Quantitative Finance Library";

  using namespace velesquant::bindings;

  // Custom exceptions (must be first)
  bind_exceptions(m);

  // Enums and utility types
  bind_enums(m);
  bind_utility(m);

  // Volatility models
  bind_sabr(m);
  bind_heston(m);
  bind_local_vol(m);

  // Interest rate models
  bind_hull_white(m);
  bind_swaption(m);
  bind_cms(m);

  // Exotic models
  bind_basket(m);
  bind_hybrid(m);

  // Numerical methods
  bind_pde_solvers(m);
  bind_trees(m);
}
