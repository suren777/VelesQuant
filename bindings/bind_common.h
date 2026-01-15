// bindings/bind_common.h - Common pybind11 includes and forward declarations
#pragma once

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace velesquant {
namespace bindings {

// Forward declarations for modular binding functions
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

} // namespace bindings
} // namespace velesquant
