#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <velesquant/local_vol/sabr.h>

namespace py = pybind11;
using namespace velesquant;

PYBIND11_MODULE(_core, m) {
    m.doc() = "Valuations Library Core Bindings";

    m.def("add", [](int i, int j) {
        return i + j;
    }, "A dummy function");

    py::class_<sabr>(m, "Sabr")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double, double>(),
             py::arg("maturity"), py::arg("forward"), py::arg("beta")=0.85,
             py::arg("alpha")=0.5, py::arg("nu")=0.25, py::arg("rho")=-0.75)
        .def(py::init<double, double, double, double, double, double, double>(),
             py::arg("maturity"), py::arg("forward"), py::arg("beta"),
             py::arg("alpha"), py::arg("nu"), py::arg("rho"), py::arg("shift"))
        .def("impliedVol", &sabr::impliedVol, py::arg("strike"))
        .def("normalVol", &sabr::normalVol, py::arg("K"))
        .def("premiumBachelier", &sabr::premiumBachelier, py::arg("strike"), py::arg("callORput")="call")
        .def("premiumBlackScholes", &sabr::premiumBlackScholes, py::arg("strike"), py::arg("callORput")="call")
        .def("localVol", &sabr::localVol, py::arg("spot"))
        .def("setParameterAlpha", &sabr::setParameterAlpha)
        .def("setParameterNu", &sabr::setParameterNu)
        .def("setParameterRho", &sabr::setParameterRho)
        .def("getParameterAlpha", &sabr::getParameterAlpha)
        .def("getParameterNu", &sabr::getParameterNu)
        .def("getParameterRho", &sabr::getParameterRho)
        .def("getMaturity", &sabr::getMaturity)
        .def("getForward", &sabr::getForward)
        .def("getBeta", &sabr::getBeta)
        .def("getShift", &sabr::getShift)
        .def("setMaturity", &sabr::setMaturity)
        .def("setForward", &sabr::setForward)
        .def("setBeta", &sabr::setBeta)
        .def_property("alpha", &sabr::getParameterAlpha, &sabr::setParameterAlpha)
        .def_property("nu", &sabr::getParameterNu, &sabr::setParameterNu)
        .def_property("rho", &sabr::getParameterRho, &sabr::setParameterRho);
}
