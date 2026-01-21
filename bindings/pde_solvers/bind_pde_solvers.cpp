// bindings/pde_solvers/bind_pde_solvers.cpp - PDE solver bindings
#include "../bind_common.h"
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/short_rate_2f_model.h>
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/hw_pde.h>
#include <velesquant/pde_solvers/sabr_pde.h>
#include <velesquant/pde_solvers/short_rate_1f_pde.h>
#include <velesquant/pde_solvers/short_rate_2f_pde.h>
#include <velesquant/volatility/sabr.h>

namespace velesquant {
namespace bindings {

void bind_pde_solvers(py::module_ &m) {
  // HWPDE - Hull-White PDE Solver (Template Instantiation)
  using HWPDE_HW = HWPDE<models::HullWhiteModel>;

  py::class_<HWPDE_HW, std::shared_ptr<HWPDE_HW>>(m, "HWPDE")
      .def(py::init<std::shared_ptr<models::HullWhiteModel>, int, double>(),
           py::arg("model"), py::arg("grid_points") = 512,
           py::arg("time_step") = 0.001)
      // Pricing methods (snake_case)
      .def("price_swaption", &HWPDE_HW::pricingSwaption, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency"))
      .def("price_bermudan", &HWPDE_HW::pricingBermudan, py::arg("expiry"),
           py::arg("tenor"), py::arg("exercises"), py::arg("strike"),
           py::arg("pay_frequency"))
      .def("price_callable_swap", &HWPDE_HW::pricingCallableSwap,
           py::arg("expiry"), py::arg("tenor"), py::arg("exercises"),
           py::arg("coupon"), py::arg("strike"), py::arg("pay_frequency"),
           py::arg("type"))
      .def("price_swap", &HWPDE_HW::pricingSwap, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency"))
      .def("price_zero_bond_option", &HWPDE_HW::pricingZBO, py::arg("expiry"),
           py::arg("maturity"), py::arg("strike"), py::arg("type"))
      .def("price_zero_bond", &HWPDE_HW::pricingZB, py::arg("maturity"))
      .def("price_coupon_bond_option", &HWPDE_HW::pricingCBO, py::arg("expiry"),
           py::arg("tenor"), py::arg("coupon"), py::arg("strike"),
           py::arg("pay_frequency"), py::arg("type"))
      .def("price_coupon_bond", &HWPDE_HW::pricingCouponBond, py::arg("expiry"),
           py::arg("tenor"), py::arg("coupon"), py::arg("pay_frequency"))
      // Analysis methods
      .def("get_swap_rate", &HWPDE_HW::getSwapRate, py::arg("expiry"),
           py::arg("tenor"), py::arg("pay_frequency"))
      .def("get_implied_vol_atm", &HWPDE_HW::getImpVolATM, py::arg("expiry"),
           py::arg("tenor"), py::arg("pay_frequency"))
      .def("get_discount_factors", &HWPDE_HW::getDFs, py::arg("time_points"))
      .def("simulate", &HWPDE_HW::simulationPDE, py::arg("times"))
      // Calibration
      .def("calibrate", &HWPDE_HW::calibrator, py::arg("discount_factor_times"),
           py::arg("discount_factors"), py::arg("swap_quotes"),
           py::arg("optimizer_params") = std::map<std::string, double>())
      // Getters
      .def("get_kappa", &HWPDE_HW::getKappa)
      .def("get_initial_rate", &HWPDE_HW::getR0)
      .def("get_time_sigmas", &HWPDE_HW::getTimeSigmas)
      .def("get_sigmas", &HWPDE_HW::getSigmas)
      .def("get_time_thetas", &HWPDE_HW::getTimeThetas)
      .def("get_thetas", &HWPDE_HW::getThetas)
      // __repr__
      .def("__repr__",
           [](HWPDE_HW &p) {
             return "<HWPDE kappa=" + std::to_string(p.getKappa()) + ">";
           })
      // NumPy array returns for zero-copy performance
      .def(
          "get_sigmas_array",
          [](HWPDE_HW &p) {
            auto vec = p.getSigmas();
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::return_value_policy::move)
      .def(
          "get_thetas_array",
          [](HWPDE_HW &p) {
            auto vec = p.getThetas();
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::return_value_policy::move)
      .def(
          "get_discount_factors_array",
          [](HWPDE_HW &p, std::vector<double> time_points) {
            auto vec = p.getDFs(time_points);
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::arg("time_points"), py::return_value_policy::move)
      .def(
          "simulate_array",
          [](HWPDE_HW &p, const std::vector<double> &times) {
            auto vec = p.simulationPDE(times);
            return py::array_t<double>(vec.size(), vec.data());
          },
          py::arg("times"), py::return_value_policy::move);

  // SABR PDE

  using SabrPDE_S = SabrPDE<Sabr>;
  using AfSABR_S = AfSABR<Sabr>;
  using AntonovSABR_S = AntonovSABR<Sabr>;

  py::class_<SabrPDE_S, std::shared_ptr<SabrPDE_S>>(m, "SabrPDE")
      .def("get_alpha", &SabrPDE_S::getAlpha)
      .def("get_beta", &SabrPDE_S::getBeta)
      .def("get_nu", &SabrPDE_S::getNu)
      .def("get_rho", &SabrPDE_S::getRho)
      .def("get_density", &SabrPDE_S::getDensity)
      .def("get_f_grid", &SabrPDE_S::getFgrid);

  py::class_<AfSABR_S, SabrPDE_S, std::shared_ptr<AfSABR_S>>(m, "AfSabr")
      .def(py::init<std::shared_ptr<Sabr>, int, int, double>(),
           py::arg("model"), py::arg("size_x") = 100, py::arg("size_t") = 100,
           py::arg("nd") = 3.0)
      .def("get_shift", &AfSABR_S::getShift);

  py::class_<AntonovSABR_S, SabrPDE_S, std::shared_ptr<AntonovSABR_S>>(
      m, "AntonovSabr")
      .def(py::init<std::shared_ptr<Sabr>, int, int, double>(),
           py::arg("model"), py::arg("size_x") = 100, py::arg("size_t") = 100,
           py::arg("nd") = 3.0);

  // ShortRate1FModel
  py::class_<models::ShortRate1FModel,
             std::shared_ptr<models::ShortRate1FModel>>(m, "ShortRate1FModel")
      .def(py::init<double, double, double, double, double, std::vector<double>,
                    std::vector<double>>(),
           py::arg("R0"), py::arg("kappa"), py::arg("alpha"), py::arg("beta"),
           py::arg("gamma"), py::arg("time_sigmas"), py::arg("sigmas"))
      .def("get_r0", &models::ShortRate1FModel::getR0)
      .def("get_kappa", &models::ShortRate1FModel::getKappa)
      .def("get_alpha", &models::ShortRate1FModel::getAlpha)
      .def("get_beta", &models::ShortRate1FModel::getBeta)
      .def("get_gamma", &models::ShortRate1FModel::getGamma);

  // ShortRate1FPDE
  using SR1FPDE_M = ShortRate1FPDE<models::ShortRate1FModel>;
  py::class_<SR1FPDE_M, std::shared_ptr<SR1FPDE_M>>(m, "ShortRate1FPDE")
      .def(py::init<std::shared_ptr<models::ShortRate1FModel>, int, double>(),
           py::arg("model"), py::arg("grid_points") = 256,
           py::arg("time_step") = 0.0025)
      .def("price_swaption", &SR1FPDE_M::pricingSwaption, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency") = 0.5)
      .def("price_bermudan", &SR1FPDE_M::pricingBermudan, py::arg("expiry"),
           py::arg("tenor"), py::arg("exercises"), py::arg("strike"),
           py::arg("pay_frequency"))
      .def("price_callable_swap", &SR1FPDE_M::pricingCallableSwap,
           py::arg("expiry"), py::arg("tenor"), py::arg("exercises"),
           py::arg("coupon"), py::arg("strike"), py::arg("pay_frequency"),
           py::arg("type"))
      .def("price_swap", &SR1FPDE_M::pricingSwap, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency"))
      .def("price_zero_bond_option", &SR1FPDE_M::pricingZBO, py::arg("expiry"),
           py::arg("maturity"), py::arg("strike"), py::arg("type"))
      .def("price_zero_bond", &SR1FPDE_M::pricingZB, py::arg("maturity"))
      .def("price_coupon_bond_option", &SR1FPDE_M::pricingCBO,
           py::arg("expiry"), py::arg("tenor"), py::arg("coupon"),
           py::arg("strike"), py::arg("pay_frequency"), py::arg("type"))
      .def("price_coupon_bond", &SR1FPDE_M::pricingCouponBond,
           py::arg("expiry"), py::arg("tenor"), py::arg("coupon"),
           py::arg("pay_frequency"))
      .def("calibrate", &SR1FPDE_M::calibrator,
           py::arg("discount_factor_times"), py::arg("discount_factors"),
           py::arg("swap_quotes"),
           py::arg("optimizer_params") = std::map<std::string, double>());

  // ShortRate2FModel
  py::class_<models::ShortRate2FModel,
             std::shared_ptr<models::ShortRate2FModel>>(m, "ShortRate2FModel")
      .def(py::init<double, double, double, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("kappa1"), py::arg("kappa2"), py::arg("lambda"),
           py::arg("time_sigma1s"), py::arg("sigma1s"), py::arg("time_sigma2s"),
           py::arg("sigma2s"), py::arg("time_alphas"), py::arg("alphas"))
      .def("get_kappa1", &models::ShortRate2FModel::getKappa1)
      .def("get_kappa2", &models::ShortRate2FModel::getKappa2)
      .def("get_lambda", &models::ShortRate2FModel::getLambda);

  // ShortRate2FPDE
  using SR2FPDE_M = ShortRate2FPDE<models::ShortRate2FModel>;
  py::class_<SR2FPDE_M, std::shared_ptr<SR2FPDE_M>>(m, "ShortRate2FPDE")
      .def(py::init<std::shared_ptr<models::ShortRate2FModel>, double, int>(),
           py::arg("model"), py::arg("time_step") = 0.020833333333,
           py::arg("grid_points") = 21)
      .def("price_swaption", &SR2FPDE_M::pricingSwaption, py::arg("expiry"),
           py::arg("tenor"), py::arg("strike"), py::arg("pay_frequency") = 0.5)
      .def("price_zero_bond", &SR2FPDE_M::pricingZB, py::arg("maturity"))
      .def("calibrate", &SR2FPDE_M::calibrator,
           py::arg("discount_factor_times"), py::arg("discount_factors"),
           py::arg("swap_quotes"),
           py::arg("optimizer_params") = std::map<std::string, double>());
}

} // namespace bindings
} // namespace velesquant
