#include "../bind_common.h"
#include <velesquant/models/log_basket.h>

namespace velesquant {
namespace bindings {

void bind_basket(py::module_ &m) {
  py::class_<lBasket>(m, "LogNormalBasket")
      .def(py::init<std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<std::vector<double>>,
                    std::vector<std::vector<double>>,
                    std::vector<std::vector<double>>>(),
           py::arg("spot"), py::arg("strike"), py::arg("maturities"),
           py::arg("forwards"), py::arg("iv"), py::arg("correlation"))
      .def("simulate", &lBasket::simulate_basket, py::arg("schedule"))
      .def("simulate_with_rebalancing", &lBasket::simulate_basketWR,
           py::arg("schedule"))
      .def("get_n_assets", &lBasket::get_nassets);
}

} // namespace bindings
} // namespace velesquant
