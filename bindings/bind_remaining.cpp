// bindings/bind_remaining.cpp - Remaining bindings (utility, CMS, trees, etc.)
#include "bind_common.h"
#include <velesquant/models/black_formula.h>
#include <velesquant/models/black_scholes.h>
#include <velesquant/models/cms.h>
#include <velesquant/models/cms_spread.h>
#include <velesquant/models/local_vol_facade.h>
#include <velesquant/models/log_basket.h>
#include <velesquant/models/quantoed_cms_spread.h>
#include <velesquant/models/simulation_facade.h>
#include <velesquant/models/swaption.h>
#include <velesquant/models/utility.h>
#include <velesquant/volatility/c_tree.h>
#include <velesquant/volatility/hhw.h>
#include <velesquant/volatility/l_vol.h>
#include <velesquant/volatility/schobzhu.h>
#include <velesquant/volatility/skew_mc.h>
#include <velesquant/volatility/termstructure.h>

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

void bind_local_vol(py::module_ &m) {
  py::class_<lVol>(m, "LocalVol")
      .def(py::init<>())
      .def(py::init<std::vector<Sabr>>(), py::arg("sabr_models"))
      .def("call_pde", &lVol::callPDE, py::arg("maturity"), py::arg("strike"),
           py::arg("n") = 100)
      .def("put_pde", &lVol::putPDE, py::arg("maturity"), py::arg("strike"),
           py::arg("n") = 100)
      .def("density", &lVol::density, py::arg("maturity"), py::arg("nt"))
      .def_property("spot", &lVol::getSpot, &lVol::setSpot);
}

void bind_swaption(py::module_ &m) {
  py::class_<swaption>(m, "Swaption")
      .def(py::init<double, double, double, double, double, double, double,
                    double>(),
           py::arg("expiry"), py::arg("tenor"), py::arg("forward"),
           py::arg("annuity"), py::arg("beta") = 0.85, py::arg("alpha") = 0.5,
           py::arg("nu") = 0.25, py::arg("rho") = -0.75)
      .def("fair_value", &swaption::swaptionFairValue, py::arg("strike"),
           py::arg("call_or_put") = OptionType::Call)
      .def("swap_fair_value", &swaption::swapFairValue, py::arg("strike"))
      .def("get_implied_vol", &swaption::getImpliedVol, py::arg("strike"))
      .def_property_readonly("alpha", &swaption::getParameterAlpha)
      .def_property_readonly("nu", &swaption::getParameterNu)
      .def_property_readonly("rho", &swaption::getParameterRho);
}

void bind_cms(py::module_ &m) {
  py::class_<cms>(m, "CMS")
      .def(py::init<double, double, double, double, double, double, double,
                    double, double, double>(),
           py::arg("expiry_sr"), py::arg("tenor_sr"), py::arg("forward_sr"),
           py::arg("annuity_sr"), py::arg("pay_cms"), py::arg("discount_cms"),
           py::arg("beta") = 0.85, py::arg("alpha") = 0.5, py::arg("nu") = 0.25,
           py::arg("rho") = -0.75)
      .def("fair_value", &cms::fairValue, py::arg("strike"),
           py::arg("call_or_put") = OptionType::Call)
      .def("get_forward", &cms::getForward)
      .def("get_discount_cms", &cms::getDiscountCMS)
      .def("get_maturity", &cms::getMaturity)
      .def("get_atm_vol", &cms::getATMvol)
      .def("get_implied_vol", &cms::getImpliedVol, py::arg("strike"))
      .def_property_readonly("alpha", &cms::getParameterAlpha)
      .def_property_readonly("nu", &cms::getParameterNu)
      .def_property_readonly("rho", &cms::getParameterRho);

  // CMS Spread
  py::class_<cms_spread>(m, "CmsSpread")
      .def(py::init([](double expiry1, double tenor1, double fwd1,
                       double annuity1, double pay1, double disc1, double beta1,
                       Eigen::MatrixXd &strikes1, Eigen::MatrixXd &quotes1,
                       CalibrationTarget type1, double expiry2, double tenor2,
                       double fwd2, double annuity2, double pay2, double disc2,
                       double beta2, Eigen::MatrixXd &strikes2,
                       Eigen::MatrixXd &quotes2, CalibrationTarget type2,
                       double corr) {
             CellMatrix sm1(strikes1), qm1(quotes1);
             CellMatrix sm2(strikes2), qm2(quotes2);
             return new cms_spread(expiry1, tenor1, fwd1, annuity1, pay1, disc1,
                                   beta1, sm1, qm1, type1, expiry2, tenor2,
                                   fwd2, annuity2, pay2, disc2, beta2, sm2, qm2,
                                   type2, corr);
           }),
           py::arg("expiry1"), py::arg("tenor1"), py::arg("fwd1"),
           py::arg("annuity1"), py::arg("pay1"), py::arg("disc1"),
           py::arg("beta1"), py::arg("strikes1"), py::arg("quotes1"),
           py::arg("type1") = CalibrationTarget::Price, py::arg("expiry2"),
           py::arg("tenor2"), py::arg("fwd2"), py::arg("annuity2"),
           py::arg("pay2"), py::arg("disc2"), py::arg("beta2"),
           py::arg("strikes2"), py::arg("quotes2"),
           py::arg("type2") = CalibrationTarget::Price, py::arg("corr"))
      .def("spread_option", &cms_spread::spreadOption, py::arg("K"),
           py::arg("a"), py::arg("b"))
      .def("simulate", static_cast<std::vector<double> (cms_spread::*)()>(
                           &cms_spread::simulationCMSs));

  // Quantoed CMS
  py::class_<quantoedCMS>(m, "QuantoedCMS")
      .def(py::init([](double expiry, double tenor, double fwd, double annuity,
                       double pay, double disc, double corFX, double atmVolFX,
                       double beta, Eigen::MatrixXd &strikes,
                       Eigen::MatrixXd &quotes, CalibrationTarget type) {
             CellMatrix sm(strikes), qm(quotes);
             return new quantoedCMS(expiry, tenor, fwd, annuity, pay, disc,
                                    corFX, atmVolFX, beta, sm, qm, type);
           }),
           py::arg("expiry"), py::arg("tenor"), py::arg("fwd"),
           py::arg("annuity"), py::arg("pay"), py::arg("disc"),
           py::arg("cor_fx"), py::arg("atm_vol_fx"), py::arg("beta"),
           py::arg("strikes"), py::arg("quotes"),
           py::arg("type") = CalibrationTarget::Price)
      .def("fair_value", &quantoedCMS::fairValue, py::arg("strike"),
           py::arg("call_or_put") = OptionType::Call)
      .def("get_forward", &quantoedCMS::getForward)
      .def("simulate", &quantoedCMS::simulation, py::arg("corr_rn"));

  // Quantoed CMS Spread
  py::class_<quantoedCMSspread>(m, "QuantoedCmsSpread")
      .def(py::init([](double expiry1, double tenor1, double fwd1,
                       double annuity1, double pay1, double disc1,
                       double corFX1, double atmVolFX1, double beta1,
                       Eigen::MatrixXd &strikes1, Eigen::MatrixXd &quotes1,
                       CalibrationTarget type1, double expiry2, double tenor2,
                       double fwd2, double annuity2, double pay2, double disc2,
                       double corFX2, double atmVolFX2, double beta2,
                       Eigen::MatrixXd &strikes2, Eigen::MatrixXd &quotes2,
                       CalibrationTarget type2, double corr) {
             CellMatrix sm1(strikes1), qm1(quotes1);
             CellMatrix sm2(strikes2), qm2(quotes2);
             return new quantoedCMSspread(
                 expiry1, tenor1, fwd1, annuity1, pay1, disc1, corFX1,
                 atmVolFX1, beta1, sm1, qm1, type1, expiry2, tenor2, fwd2,
                 annuity2, pay2, disc2, corFX2, atmVolFX2, beta2, sm2, qm2,
                 type2, corr);
           }),
           py::arg("expiry1"), py::arg("tenor1"), py::arg("fwd1"),
           py::arg("annuity1"), py::arg("pay1"), py::arg("disc1"),
           py::arg("cor_fx1"), py::arg("atm_vol_fx1"), py::arg("beta1"),
           py::arg("strikes1"), py::arg("quotes1"),
           py::arg("type1") = CalibrationTarget::Price, py::arg("expiry2"),
           py::arg("tenor2"), py::arg("fwd2"), py::arg("annuity2"),
           py::arg("pay2"), py::arg("disc2"), py::arg("cor_fx2"),
           py::arg("atm_vol_fx2"), py::arg("beta2"), py::arg("strikes2"),
           py::arg("quotes2"), py::arg("type2") = CalibrationTarget::Price,
           py::arg("corr"))
      .def("simulate", &quantoedCMSspread::simulationQuantoedCMSs,
           py::arg("cr1"), py::arg("cr2"));
}

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

void bind_hybrid(py::module_ &m) {
  // HHW - Hybrid Hull-White
  py::class_<HHW>(m, "HHW")
      .def(py::init<double, double, double, double, double, double, double,
                    double, double>(),
           py::arg("s0"), py::arg("v0"), py::arg("initial_rate"),
           py::arg("kappa"), py::arg("eta"), py::arg("rho"), py::arg("sigma1"),
           py::arg("sigma2"), py::arg("a"))
      .def("price",
           static_cast<double (HHW::*)(double, double) const>(&HHW::HHWPrice),
           py::arg("maturity"), py::arg("strike"));

  // Schobel-Zhu
  py::class_<SchobelZhu>(m, "SchobelZhu")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("spot"), py::arg("var0"), py::arg("kappa"), py::arg("theta"),
           py::arg("xi"), py::arg("rho"))
      .def("price", &SchobelZhu::SchobelPrice, py::arg("maturity"),
           py::arg("forward"), py::arg("strike"))
      .def("simulate", &SchobelZhu::simulation, py::arg("times"),
           py::arg("forwards"))
      .def("calibrate", &SchobelZhu::calibrator, py::arg("maturities"),
           py::arg("forwards"), py::arg("strikes"), py::arg("market_quotes"),
           py::arg("target") = CalibrationTarget::Price)
      .def_property("var0", &SchobelZhu::getParameterVar0,
                    &SchobelZhu::setParameterVar0)
      .def_property("kappa", &SchobelZhu::getParameterKappa,
                    &SchobelZhu::setParameterKappa)
      .def_property("theta", &SchobelZhu::getParameterTheta,
                    &SchobelZhu::setParameterTheta)
      .def_property("xi", &SchobelZhu::getParameterXi,
                    &SchobelZhu::setParameterXi)
      .def_property("rho", &SchobelZhu::getParameterRho,
                    &SchobelZhu::setParameterRho);

  // Termstructure
  py::class_<Termstructure>(m, "Termstructure")
      .def(py::init<std::vector<double>, std::vector<double>, CalendarType,
                    DayCounterType>(),
           py::arg("days"), py::arg("rate"), py::arg("calendar"),
           py::arg("daycount"))
      .def(py::init<std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>, CalendarType,
                    DayCounterType>(),
           py::arg("days"), py::arg("rate"), py::arg("qdays"),
           py::arg("dividend"), py::arg("calendar"), py::arg("daycount"))
      .def(
          "discount",
          static_cast<double (Termstructure::*)(int)>(&Termstructure::discount),
          py::arg("date"))
      .def("rate",
           static_cast<double (Termstructure::*)(int, int)>(
               &Termstructure::rate),
           py::arg("date"), py::arg("tenor"))
      .def("dividend",
           static_cast<double (Termstructure::*)(int, int)>(
               &Termstructure::divident),
           py::arg("date"), py::arg("tenor"));

  // Facade function
  m.def("lv_export", &lvExport, py::arg("maturities"), py::arg("forwards"),
        py::arg("betas"), py::arg("alphas"), py::arg("nus"), py::arg("rhos"),
        py::arg("spot"), py::arg("times"));
}

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
