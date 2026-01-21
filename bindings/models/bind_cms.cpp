#include "../bind_common.h"
#include <velesquant/models/cms.h>
#include <velesquant/models/cms_spread.h>
#include <velesquant/models/quantoed_cms_spread.h>

namespace velesquant {
namespace bindings {

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

} // namespace bindings
} // namespace velesquant
