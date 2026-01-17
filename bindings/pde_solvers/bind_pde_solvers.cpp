// bindings/pde_solvers/bind_pde_solvers.cpp - PDE solver bindings
#include "../bind_common.h"
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/hw_pde.h>
#include <velesquant/pde_solvers/sabr_pde.h>
#include <velesquant/pde_solvers/short_rate_1f_pde.h>
#include <velesquant/pde_solvers/short_rate_2f_pde.h>

namespace velesquant {
namespace bindings {

void bind_pde_solvers(py::module_ &m) {
  // HWPDE - Hull-White PDE Solver
  py::class_<HWPDE, std::shared_ptr<HWPDE>>(m, "HWPDE")
      .def(py::init<double, double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>, int, double>(),
           py::arg("initial_rate"), py::arg("kappa"), py::arg("time_sigmas"),
           py::arg("sigmas"), py::arg("time_thetas"), py::arg("thetas"),
           py::arg("grid_points") = 512, py::arg("time_step") = 0.001)
      .def(py::init<double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>, int, double>(),
           py::arg("kappa"), py::arg("time_sigmas"), py::arg("sigmas"),
           py::arg("discount_factor_times"), py::arg("discount_factors"),
           py::arg("grid_points") = 512, py::arg("time_step") = 0.001)
      // Pricing methods (snake_case)
      .def("price_swaption", &HWPDE::pricingSwaption, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency"))
      .def("price_bermudan", &HWPDE::pricingBermudan, py::arg("expiry"),
           py::arg("tenor"), py::arg("exercises"), py::arg("strike"),
           py::arg("pay_frequency"))
      .def("price_callable_swap", &HWPDE::pricingCallableSwap,
           py::arg("expiry"), py::arg("tenor"), py::arg("exercises"),
           py::arg("coupon"), py::arg("strike"), py::arg("pay_frequency"),
           py::arg("type"))
      .def("price_swap", &HWPDE::pricingSwap, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency"))
      .def("price_zero_bond_option", &HWPDE::pricingZBO, py::arg("expiry"),
           py::arg("maturity"), py::arg("strike"), py::arg("type"))
      .def("price_zero_bond", &HWPDE::pricingZB, py::arg("maturity"))
      .def("price_coupon_bond_option", &HWPDE::pricingCBO, py::arg("expiry"),
           py::arg("tenor"), py::arg("coupon"), py::arg("strike"),
           py::arg("pay_frequency"), py::arg("type"))
      .def("price_coupon_bond", &HWPDE::pricingCouponBond, py::arg("expiry"),
           py::arg("tenor"), py::arg("coupon"), py::arg("pay_frequency"))
      // Analysis methods
      .def("get_swap_rate", &HWPDE::getSwapRate, py::arg("expiry"),
           py::arg("tenor"), py::arg("pay_frequency"))
      .def("get_implied_vol_atm", &HWPDE::getImpVolATM, py::arg("expiry"),
           py::arg("tenor"), py::arg("pay_frequency"))
      .def("get_discount_factors", &HWPDE::getDFs, py::arg("time_points"))
      .def("simulate", &HWPDE::simulationPDE, py::arg("times"))
      // Calibration
      .def("calibrate", &HWPDE::calibrator, py::arg("discount_factor_times"),
           py::arg("discount_factors"), py::arg("swap_quotes"))
      // Getters
      .def("get_kappa", &HWPDE::getKappa)
      .def("get_initial_rate", &HWPDE::getR0)
      .def("get_time_sigmas", &HWPDE::getTimeSigmas)
      .def("get_sigmas", &HWPDE::getSigmas)
      .def("get_time_thetas", &HWPDE::getTimeThetas)
      .def("get_thetas", &HWPDE::getThetas)
      // __repr__
      .def("__repr__",
           [](HWPDE &p) {
             return "<HWPDE kappa=" + std::to_string(p.getKappa()) + ">";
           })
      // NumPy array returns for zero-copy performance
      .def(
          "get_sigmas_array",
          [](HWPDE &p) {
            auto vec = p.getSigmas();
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::return_value_policy::move)
      .def(
          "get_thetas_array",
          [](HWPDE &p) {
            auto vec = p.getThetas();
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::return_value_policy::move)
      .def(
          "get_discount_factors_array",
          [](HWPDE &p, std::vector<double> time_points) {
            auto vec = p.getDFs(time_points);
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::arg("time_points"), py::return_value_policy::move)
      .def(
          "simulate_array",
          [](HWPDE &p, const std::vector<double> &times) {
            auto vec = p.simulationPDE(times);
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::arg("times"), py::return_value_policy::move);

  // SABR PDE
  py::class_<sabr_pde>(m, "SabrPDE")
      .def("get_alpha", &sabr_pde::getAlpha)
      .def("get_beta", &sabr_pde::getBeta)
      .def("get_nu", &sabr_pde::getNu)
      .def("get_rho", &sabr_pde::getRho)
      .def("get_density", &sabr_pde::getDensity)
      .def("get_f_grid", &sabr_pde::getFgrid);

  py::class_<afSABR, sabr_pde>(m, "AfSabr")
      .def(py::init<double, double, double, double, double, double, double, int,
                    int, double>(),
           py::arg("alpha"), py::arg("beta"), py::arg("nu"), py::arg("rho"),
           py::arg("shift"), py::arg("maturity"), py::arg("F"),
           py::arg("size_x"), py::arg("size_t"), py::arg("nd"))
      .def("get_shift", &afSABR::getShift);

  py::class_<antonovSABR, sabr_pde>(m, "AntonovSabr")
      .def(py::init<double, double, double, double, double, double, int, int,
                    double>(),
           py::arg("alpha"), py::arg("beta"), py::arg("nu"), py::arg("rho"),
           py::arg("maturity"), py::arg("F"), py::arg("size_x"),
           py::arg("size_t"), py::arg("nd"));

  // ShortRate1FPDE
  py::class_<ShortRate1FPDE>(m, "ShortRate1FPDE")
      .def(py::init<double, double, double, double, double, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("initial_rate"), py::arg("kappa"), py::arg("alpha"),
           py::arg("beta"), py::arg("gamma"), py::arg("time_sigmas"),
           py::arg("sigmas"), py::arg("time_thetas"), py::arg("thetas"))
      .def("price_swaption", &ShortRate1FPDE::pricingSwaption,
           py::arg("expiry"), py::arg("tenor"), py::arg("strike"),
           py::arg("pay_frequency") = 0.5)
      .def("price_zero_bond_option", &ShortRate1FPDE::pricingZBO,
           py::arg("expiry"), py::arg("maturity"), py::arg("strike"),
           py::arg("type"))
      .def("price_coupon_bond_option", &ShortRate1FPDE::pricingCBO,
           py::arg("expiry"), py::arg("tenor"), py::arg("coupon"),
           py::arg("strike"), py::arg("pay_frequency"), py::arg("type"));

  // ShortRate2FPDE
  py::class_<ShortRate2FPDE>(m, "ShortRate2FPDE")
      .def(py::init<double, double, double, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("kappa1"), py::arg("kappa2"), py::arg("lambda"),
           py::arg("time_sigma1s"), py::arg("sigma1s"), py::arg("time_sigma2s"),
           py::arg("sigma2s"), py::arg("time_alphas"), py::arg("alphas"))
      .def("price_swaption", &ShortRate2FPDE::pricingSwaption,
           py::arg("expiry"), py::arg("tenor"), py::arg("strike"),
           py::arg("pay_frequency") = 0.5)
      .def("calibrate", &ShortRate2FPDE::calibrator,
           py::arg("discount_factor_times"), py::arg("discount_factors"),
           py::arg("swap_quotes"));
}

} // namespace bindings
} // namespace velesquant
