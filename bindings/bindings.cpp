#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <velesquant/models/black_formula.h>
#include <velesquant/models/black_scholes.h>
#include <velesquant/models/cms.h>
#include <velesquant/models/cms_spread.h>
#include <velesquant/models/heston_facade.h>
#include <velesquant/models/hw.h>
#include <velesquant/models/local_vol_facade.h>
#include <velesquant/models/log_basket.h>
#include <velesquant/models/quantoed_cms_spread.h>
#include <velesquant/models/sabr_facade.h>
#include <velesquant/models/simulation_facade.h>
#include <velesquant/models/swaption.h>
#include <velesquant/models/utility.h>
#include <velesquant/pde_solvers/hw_pde.h>
#include <velesquant/pde_solvers/sabr_pde.h>
#include <velesquant/pde_solvers/short_rate_1f_pde.h>
#include <velesquant/pde_solvers/short_rate_2f_pde.h>
#include <velesquant/volatility/c_tree.h>
#include <velesquant/volatility/hhw.h>
#include <velesquant/volatility/l_vol.h>
#include <velesquant/volatility/s_vol.h>
#include <velesquant/volatility/sabr.h>
#include <velesquant/volatility/schobzhu.h>
#include <velesquant/volatility/skew_mc.h>
#include <velesquant/volatility/termstructure.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace velesquant;
// using namespace velesquant::xlw; // Removed xlw namespace usage

PYBIND11_MODULE(native, m) {
  m.doc() = "Valuations Library Core Bindings";

  // Enums for CTree

  py::enum_<tType>(m, "TreeType")
      .value("Recombining", tType::recomb)
      .value("NonRecombining", tType::nonrecomb)
      .export_values();

  py::enum_<exStyle>(m, "ExerciseStyle")
      .value("American", exStyle::American)
      .value("European", exStyle::European)
      .value("Bermudan", exStyle::Bermudan)
      .export_values();

  // Enums for Termstructure
  py::enum_<CalendarType>(m, "CalendarType")
      .value("UnitedKingdom", CalendarType::CAL_UnitedKingdom)
      .value("Argentina", CalendarType::CAL_Argentina)
      .value("Australia", CalendarType::CAL_Australia)
      .value("BespokeCalendar", CalendarType::CAL_BespokeCalendar)
      .value("Brazil", CalendarType::CAL_Brazil)
      .value("Canada", CalendarType::CAL_Canada)
      .value("China", CalendarType::CAL_China)
      .value("CzechRepublic", CalendarType::CAL_CzechRepublic)
      .value("Denmark", CalendarType::CAL_Denmark)
      .value("Finland", CalendarType::CAL_Finland)
      .value("Germany", CalendarType::CAL_Germany)
      .value("HongKong", CalendarType::CAL_HongKong)
      .value("Hungary", CalendarType::CAL_Hungary)
      .value("Iceland", CalendarType::CAL_Iceland)
      .value("India", CalendarType::CAL_India)
      .value("Indonesia", CalendarType::CAL_Indonesia)
      .value("Italy", CalendarType::CAL_Italy)
      .value("Japan", CalendarType::CAL_Japan)
      .value("Ukraine", CalendarType::CAL_Ukraine)
      .value("Mexico", CalendarType::CAL_Mexico)
      .value("NewZealand", CalendarType::CAL_NewZealand)
      .value("Norway", CalendarType::CAL_Norway)
      .value("Poland", CalendarType::CAL_Poland)
      .value("Russia", CalendarType::CAL_Russia)
      .value("SaudiArabia", CalendarType::CAL_SaudiArabia)
      .value("Singapore", CalendarType::CAL_Singapore)
      .value("Slovakia", CalendarType::CAL_Slovakia)
      .value("SouthAfrica", CalendarType::CAL_SouthAfrica)
      .value("SouthKorea", CalendarType::CAL_SouthKorea)
      .value("Sweden", CalendarType::CAL_Sweden)
      .value("Switzerland", CalendarType::CAL_Switzerland)
      .value("Taiwan", CalendarType::CAL_Taiwan)
      .value("TARGET", CalendarType::CAL_TARGET)
      .value("Turkey", CalendarType::CAL_Turkey)
      .export_values();

  py::enum_<DayCounterType>(m, "DayCounterType")
      .value("Actual360", DayCounterType::DC_Actual360)
      .value("Actual365Fixed", DayCounterType::DC_Actual365Fixed)
      .value("ActualActual", DayCounterType::DC_ActualActual)
      .value("Actual365NoLeap", DayCounterType::DC_Actual365NoLeap)
      .value("Business252", DayCounterType::DC_Business252)
      .value("OneDayCounter", DayCounterType::DC_OneDayCounter)
      .value("SimpleDayCounter", DayCounterType::DC_SimpleDayCounter)
      .value("Thirty360", DayCounterType::DC_Thirty360)
      .export_values();

  py::enum_<OptionType>(m, "OptionType")
      .value("Call", OptionType::Call)
      .value("Put", OptionType::Put)
      .export_values();

  py::enum_<CalibrationTarget>(m, "CalibrationTarget")
      .value("Price", CalibrationTarget::Price)
      .value("Volatility", CalibrationTarget::Volatility)
      .export_values();

  // defSwap struct binding for calibration inputs
  py::class_<defSwap>(m, "DefSwap")
      .def(py::init<>())
      .def_readwrite("Expiry", &defSwap::Expiry)
      .def_readwrite("Tenor", &defSwap::Tenor)
      .def_readwrite("Frequency", &defSwap::Frequency)
      .def_readwrite("SwapRate", &defSwap::SwapRate)
      .def_readwrite("VolATM", &defSwap::VolATM)
      .def_readwrite("Value", &defSwap::Value);

  // Utility functions
  m.def("cdf_normal", &cdf_normal,
        "Cumulative distribution function for normal distribution");
  m.def("pdf_normal", &pdf_normal,
        "Probability density function for normal distribution");
  m.def("implied_vol", &implied_vol, py::arg("Maturity"), py::arg("Forward"),
        py::arg("Strike"), py::arg("Price"),
        "Calculate implied volatility for a Black-Scholes price");

  // Equity & Volatility - Black-Scholes / Black Formula
  m.def("black_scholes_call", &BlackScholesCall, py::arg("Spot"),
        py::arg("Strike"), py::arg("r"), py::arg("d"), py::arg("Vol"),
        py::arg("Expiry"));
  m.def("black_scholes_call_vega", &BlackScholesCallVega, py::arg("Spot"),
        py::arg("Strike"), py::arg("r"), py::arg("d"), py::arg("Vol"),
        py::arg("Expiry"));
  m.def("black_formula_call", &BlackFormulaCall, py::arg("Forward"),
        py::arg("Strike"), py::arg("Vol"), py::arg("Expiry"));
  m.def("black_formula_call_vega", &BlackFormulaCallVega, py::arg("Forward"),
        py::arg("Strike"), py::arg("Vol"), py::arg("Expiry"));

  // SABR
  py::class_<Sabr>(m, "Sabr")
      .def(py::init<>())
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("maturity"), py::arg("forward"), py::arg("beta") = 0.85,
           py::arg("alpha") = 0.5, py::arg("nu") = 0.25, py::arg("rho") = -0.75)
      .def(py::init<double, double, double, double, double, double, double>(),
           py::arg("maturity"), py::arg("forward"), py::arg("beta"),
           py::arg("alpha"), py::arg("nu"), py::arg("rho"), py::arg("shift"))
      .def("impliedVol", &Sabr::impliedVol, py::arg("strike"))
      .def("normalVol", &Sabr::normalVol, py::arg("K"))
      .def("premiumBachelier", &Sabr::premiumBachelier, py::arg("strike"),
           py::arg("callORput") = velesquant::OptionType::Call)
      .def("premiumBlackScholes", &Sabr::premiumBlackScholes, py::arg("strike"),
           py::arg("callORput") = velesquant::OptionType::Call)
      .def("localVol", &Sabr::localVol, py::arg("spot"))
      .def_property("alpha", &Sabr::getParameterAlpha, &Sabr::setParameterAlpha)
      .def_property("nu", &Sabr::getParameterNu, &Sabr::setParameterNu)
      .def_property("rho", &Sabr::getParameterRho, &Sabr::setParameterRho)
      .def_property("maturity", &Sabr::getMaturity, &Sabr::setMaturity)
      .def_property("forward", &Sabr::getForward, &Sabr::setForward)
      .def_property("beta", &Sabr::getBeta, &Sabr::setBeta)
      .def(
          "calibrate",
          [](Sabr &self, const Eigen::MatrixXd &strikes,
             const Eigen::MatrixXd &quotes, CalibrationTarget quoteType) {
            std::vector<double> s_vec, q_vec;
            int rows = strikes.rows();
            int cols =
                strikes.cols(); // Should be 1 typically for strikes inputs
            // Logic adapted from model.cpp: handle 1D arrays (N x 1)
            for (int i = 0; i < rows; i++) {
              for (int j = 0; j < cols; j++) {
                if (!std::isnan(strikes(i, j)) && !std::isnan(quotes(i, j))) {
                  s_vec.push_back(strikes(i, j));
                  q_vec.push_back(quotes(i, j));
                }
              }
            }
            if (s_vec.empty()) {
              throw std::runtime_error("No valid data points for calibration");
            }
            self.calibrator(s_vec, q_vec, quoteType);
            // Return calibrated parameters convenient map/dict
            return py::dict("alpha"_a = self.getParameterAlpha(),
                            "nu"_a = self.getParameterNu(),
                            "rho"_a = self.getParameterRho());
          },
          py::arg("strikes"), py::arg("quotes"),
          py::arg("quote_type") = CalibrationTarget::Volatility);

  // Heston (sVol)
  py::class_<sVol>(m, "Heston")
      .def(py::init<>())
      .def(py::init<double, double, double, double, double, double, int>(),
           py::arg("spot"), py::arg("var0"), py::arg("kappa"), py::arg("theta"),
           py::arg("xi"), py::arg("rho"), py::arg("seed") = 42)
      .def("hestonPrice", &sVol::hestonPrice, py::arg("maturity"),
           py::arg("forward"), py::arg("strike"), py::arg("optType") = "call")
      .def("simulationHeston", &sVol::simulationHeston, py::arg("times"),
           py::arg("forwards"))
      .def_property_readonly("var0", &sVol::getParameterVar0)
      .def_property_readonly("kappa", &sVol::getParameterKappa)
      .def_property_readonly("theta", &sVol::getParameterTheta)
      .def_property_readonly("xi", &sVol::getParameterXi)
      .def_property_readonly("rho", &sVol::getParameterRho)
      .def(
          "calibrate",
          [](sVol &self, const Eigen::MatrixXd &maturities,
             const Eigen::MatrixXd &forwards, const Eigen::MatrixXd &strikes,
             const Eigen::MatrixXd &quotes, std::string quoteType) {
            std::vector<double> m_vec, f_vec, k_vec, q_vec;
            int m = strikes.rows();
            int k = quotes.rows();
            int l = quotes.cols();

            if ((m == k) && (l == 1)) {
              for (int i = 0; i < m; i++) {
                m_vec.push_back(maturities(i, 0));
                f_vec.push_back(forwards(i, 0));
                k_vec.push_back(strikes(i, 0));
                q_vec.push_back(quotes(i, 0));
              }
            } else {
              m = strikes.cols();
              if (m != l)
                throw std::runtime_error("Strikes don't match vols dimensions");

              for (int i = 0; i < k; i++)
                for (int j = 0; j < l; j++)
                  if (quotes(i, j) > 1e-10) {
                    m_vec.push_back(maturities(i, 0));
                    f_vec.push_back(forwards(i, 0));
                    k_vec.push_back(strikes(0, j));
                    q_vec.push_back(quotes(i, j));
                  }
            }

            self.calibrator(m_vec, f_vec, k_vec, q_vec, quoteType);
            return py::dict("var0"_a = self.getParameterVar0(),
                            "kappa"_a = self.getParameterKappa(),
                            "theta"_a = self.getParameterTheta(),
                            "xi"_a = self.getParameterXi(),
                            "rho"_a = self.getParameterRho());
          },
          py::arg("maturities"), py::arg("forwards"), py::arg("strikes"),
          py::arg("quotes"), py::arg("quote_type") = "Price");

  // Local Volatility (lVol)
  py::class_<lVol>(m, "LocalVol")
      .def(py::init<>())
      .def(py::init<std::vector<Sabr>>(), py::arg("sabrModels"))
      .def("callPDE", &lVol::callPDE, py::arg("maturity"), py::arg("strike"),
           py::arg("N") = 100)
      .def("putPDE", &lVol::putPDE, py::arg("maturity"), py::arg("strike"),
           py::arg("N") = 100)
      .def("density", &lVol::density, py::arg("maturity"), py::arg("Nt"))
      .def_property("spot", &lVol::getSpot, &lVol::setSpot);

  // Hull-White (HullWhite)
  py::class_<HullWhite>(m, "HullWhite")
      .def(py::init<double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>>(),
           py::arg("kappa"), py::arg("timeSigmas"), py::arg("sigmas"),
           py::arg("timeDFs"), py::arg("DFs"))
      .def("optionBond", &HullWhite::optionBond, py::arg("Expiry"),
           py::arg("Maturity"), py::arg("Strike"), py::arg("type"))
      .def("swaption", &HullWhite::swaption, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("Strike"), py::arg("PayFrequency") = 0.5)
      .def("simulation", &HullWhite::simulation, py::arg("times"))
      .def("getParameterKappa", &HullWhite::getKappa)
      .def("getSigmas", &HullWhite::getSigmas)
      .def("getTimeSigmas", &HullWhite::getTimeSigmas)
      .def("ZC", &HullWhite::ZC, py::arg("Expiry"))
      .def("calibrate", &HullWhite::calibrator, py::arg("swapQuotes"),
           py::arg("target") = CalibrationTarget::Volatility)
      .def("swaption_vec",
           py::vectorize([](HullWhite &self, double expiry, double tenor,
                            double strike, double pay_freq) {
             return self.swaption(expiry, tenor, strike, pay_freq);
           }),
           py::arg("Expiry"), py::arg("Tenor"), py::arg("Strike"),
           py::arg("PayFrequency") = 0.5)
      .def("ZC_vec", py::vectorize([](HullWhite &self, double expiry) {
             return self.ZC(expiry);
           }),
           py::arg("Expiry"));

  // Swaption
  py::class_<swaption>(m, "Swaption")
      .def(py::init<double, double, double, double, double, double, double,
                    double>(),
           py::arg("expiry"), py::arg("tenor"), py::arg("forward"),
           py::arg("annuity"), py::arg("beta") = 0.85, py::arg("alpha") = 0.5,
           py::arg("nu") = 0.25, py::arg("rho") = -0.75)
      .def("swaptionFairValue", &swaption::swaptionFairValue, py::arg("strike"),
           py::arg("callORput") = velesquant::OptionType::Call)
      .def("swapFairValue", &swaption::swapFairValue, py::arg("strike"))
      .def("getImpliedVol", &swaption::getImpliedVol, py::arg("strike"))
      .def_property_readonly("alpha", &swaption::getParameterAlpha)
      .def_property_readonly("nu", &swaption::getParameterNu)
      .def_property_readonly("rho", &swaption::getParameterRho);

  // CMS
  py::class_<cms>(m, "CMS")
      .def(py::init<double, double, double, double, double, double, double,
                    double, double, double>(),
           py::arg("expirySR"), py::arg("tenorSR"), py::arg("forwardSR"),
           py::arg("annuitySR"), py::arg("payCMS"), py::arg("discountCMS"),
           py::arg("beta") = 0.85, py::arg("alpha") = 0.5, py::arg("nu") = 0.25,
           py::arg("rho") = -0.75)
      .def("fairValue", &cms::fairValue, py::arg("strike"),
           py::arg("callORput") = velesquant::OptionType::Call)
      .def("getForward", &cms::getForward)
      .def("getDiscountCMS", &cms::getDiscountCMS)
      .def("getMaturity", &cms::getMaturity)
      .def("getATMvol", &cms::getATMvol)
      .def("getImpliedVol", &cms::getImpliedVol, py::arg("strike"))
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
           py::arg("type1") = velesquant::CalibrationTarget::Price,
           py::arg("expiry2"), py::arg("tenor2"), py::arg("fwd2"),
           py::arg("annuity2"), py::arg("pay2"), py::arg("disc2"),
           py::arg("beta2"), py::arg("strikes2"), py::arg("quotes2"),
           py::arg("type2") = velesquant::CalibrationTarget::Price,
           py::arg("corr"))
      .def("spreadOption", &cms_spread::spreadOption, py::arg("K"),
           py::arg("a"), py::arg("b"))
      .def("simulationCMSs", static_cast<std::vector<double> (cms_spread::*)()>(
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
           py::arg("corFX"), py::arg("atmVolFX"), py::arg("beta"),
           py::arg("strikes"), py::arg("quotes"),
           py::arg("type") = velesquant::CalibrationTarget::Price)
      .def("fairValue", &quantoedCMS::fairValue, py::arg("strike"),
           py::arg("callORput") = velesquant::OptionType::Call)
      .def("getForward", &quantoedCMS::getForward)
      .def("simulation", &quantoedCMS::simulation, py::arg("corrRN"));

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
           py::arg("corFX1"), py::arg("atmVolFX1"), py::arg("beta1"),
           py::arg("strikes1"), py::arg("quotes1"),
           py::arg("type1") = velesquant::CalibrationTarget::Price,
           py::arg("expiry2"), py::arg("tenor2"), py::arg("fwd2"),
           py::arg("annuity2"), py::arg("pay2"), py::arg("disc2"),
           py::arg("corFX2"), py::arg("atmVolFX2"), py::arg("beta2"),
           py::arg("strikes2"), py::arg("quotes2"),
           py::arg("type2") = velesquant::CalibrationTarget::Price,
           py::arg("corr"))
      .def("simulationQuantoedCMSs", &quantoedCMSspread::simulationQuantoedCMSs,
           py::arg("cr1"), py::arg("cr2"));

  // Log-Normal Basket
  py::class_<lBasket>(m, "LogNormalBasket")
      .def(py::init<std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<std::vector<double>>,
                    std::vector<std::vector<double>>,
                    std::vector<std::vector<double>>>(),
           py::arg("Spot"), py::arg("Strike"), py::arg("Maturities"),
           py::arg("Forwards"), py::arg("IV"), py::arg("correlation"))
      .def("simulate_basket", &lBasket::simulate_basket, py::arg("schedule"))
      .def("simulate_basketWR", &lBasket::simulate_basketWR,
           py::arg("schedule"))
      .def("get_nassets", &lBasket::get_nassets);

  // HHW - Hybrid Hull-White
  py::class_<HHW>(m, "HHW")
      .def(py::init<double, double, double, double, double, double, double,
                    double, double>(),
           py::arg("s0"), py::arg("v0"), py::arg("r0"), py::arg("kappa"),
           py::arg("eta"), py::arg("rho"), py::arg("sigma1"), py::arg("sigma2"),
           py::arg("a"))
      .def("HHWPrice",
           static_cast<double (HHW::*)(double, double) const>(&HHW::HHWPrice),
           py::arg("maturity"), py::arg("strike"));

  // Schobel-Zhu
  py::class_<SchobelZhu>(m, "SchobelZhu")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("spot"), py::arg("var0"), py::arg("kappa"), py::arg("theta"),
           py::arg("xi"), py::arg("rho"))
      .def("SchobelPrice", &SchobelZhu::SchobelPrice, py::arg("maturity"),
           py::arg("forward"), py::arg("strike"))
      .def("simulation", &SchobelZhu::simulation, py::arg("times"),
           py::arg("forwards"))
      .def("calibrator", &SchobelZhu::calibrator, py::arg("maturitys"),
           py::arg("forwards"), py::arg("strikes"), py::arg("marketQuotes"),
           py::arg("target") = velesquant::CalibrationTarget::Price)
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
           py::arg("divident"), py::arg("calendar"), py::arg("daycount"))
      .def(
          "discount",
          static_cast<double (Termstructure::*)(int)>(&Termstructure::discount),
          py::arg("date"))
      .def("rate",
           static_cast<double (Termstructure::*)(int, int)>(
               &Termstructure::rate),
           py::arg("date"), py::arg("tenor"))
      .def("divident",
           static_cast<double (Termstructure::*)(int, int)>(
               &Termstructure::divident),
           py::arg("date"), py::arg("tenor"));

  // SABR PDE
  py::class_<sabr_pde>(m, "SabrPDE")
      .def("getAlpha", &sabr_pde::getAlpha)
      .def("getBeta", &sabr_pde::getBeta)
      .def("getNu", &sabr_pde::getNu)
      .def("getRho", &sabr_pde::getRho)
      .def("getDensity", &sabr_pde::getDensity)
      .def("getFgrid", &sabr_pde::getFgrid);

  py::class_<afSABR, sabr_pde>(m, "AfSabr")
      .def(py::init<double, double, double, double, double, double, double, int,
                    int, double>(),
           py::arg("alpha"), py::arg("beta"), py::arg("nu"), py::arg("rho"),
           py::arg("shift"), py::arg("maturity"), py::arg("F"),
           py::arg("sizeX"), py::arg("sizeT"), py::arg("nd"))
      .def("getShift", &afSABR::getShift);

  py::class_<antonovSABR, sabr_pde>(m, "AntonovSabr")
      .def(py::init<double, double, double, double, double, double, int, int,
                    double>(),
           py::arg("alpha"), py::arg("beta"), py::arg("nu"), py::arg("rho"),
           py::arg("maturity"), py::arg("F"), py::arg("sizeX"),
           py::arg("sizeT"), py::arg("nd"));

  // CTree - Binomial/Trinomial Trees
  py::class_<CTree>(m, "CTree")
      .def(py::init<>())
      .def(py::init<double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("S"), py::arg("T"), py::arg("F"), py::arg("IV"),
           py::arg("r") = std::vector<double>(),
           py::arg("q") = std::vector<double>())
      .def("calculateBinomial",
           static_cast<double (CTree::*)(double, double, int, exStyle,
                                         OptionType, tType)>(
               &CTree::calculateBinomial),
           py::arg("strike"), py::arg("Maturity"), py::arg("Nnodes"),
           py::arg("style"), py::arg("pay"), py::arg("tree"))
      .def("calculateTrinomial", &CTree::calculateTrinomial, py::arg("strike"),
           py::arg("Maturity"), py::arg("Nnodes"), py::arg("style"),
           py::arg("pay"), py::arg("tree"));

  // skewMC - Skew Monte Carlo
  py::class_<skewMC>(m, "SkewMC")
      .def(py::init<>())
      .def(py::init<std::vector<Sabr>>(), py::arg("sabrModels"))
      .def("simulation", &skewMC::simulation, py::arg("times"), py::arg("spot"),
           py::arg("kappa"));

  // HWPDE - Hull-White PDE Solver
  py::class_<HWPDE>(m, "HWPDE")
      .def(py::init<double, double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>>(),
           py::arg("R0"), py::arg("kappa"), py::arg("timeSigmas"),
           py::arg("sigmas"), py::arg("timeThetas"), py::arg("thetas"))
      .def(py::init<double, std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>>(),
           py::arg("kappa"), py::arg("timeSigmas"), py::arg("sigmas"),
           py::arg("timeDF"), py::arg("DF"))
      // Pricing methods
      .def("pricingSwaption", &HWPDE::pricingSwaption, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("Strike"), py::arg("PayFrequency"))
      .def("pricingBermudan", &HWPDE::pricingBermudan, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("Exercises"), py::arg("Strike"),
           py::arg("PayFrequency"))
      .def("pricingCallableSwap", &HWPDE::pricingCallableSwap,
           py::arg("Expiry"), py::arg("Tenor"), py::arg("Exercises"),
           py::arg("Coupon"), py::arg("Strike"), py::arg("PayFrequency"),
           py::arg("type"))
      .def("pricingSwap", &HWPDE::pricingSwap, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("Strike"), py::arg("PayFrequency"))
      .def("pricingZBO", &HWPDE::pricingZBO, py::arg("Expiry"),
           py::arg("Maturity"), py::arg("Strike"), py::arg("type"))
      .def("pricingZB", &HWPDE::pricingZB, py::arg("Maturity"))
      .def("pricingCBO", &HWPDE::pricingCBO, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("Coupon"), py::arg("Strike"),
           py::arg("PayFrequency"), py::arg("type"))
      .def("pricingCouponBond", &HWPDE::pricingCouponBond, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("Coupon"), py::arg("PayFrequency"))
      // Analysis methods
      .def("getSwapRate", &HWPDE::getSwapRate, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("PayFrequency"))
      .def("getImpVolATM", &HWPDE::getImpVolATM, py::arg("Expiry"),
           py::arg("Tenor"), py::arg("PayFrequency"))
      .def("getDFs", &HWPDE::getDFs, py::arg("timePoints"))
      .def("simulationPDE", &HWPDE::simulationPDE, py::arg("times"))
      // Calibration
      .def("calibrate", &HWPDE::calibrator, py::arg("timeDFs"), py::arg("DFs"),
           py::arg("swapQuotes"))
      // Getters
      .def("getKappa", &HWPDE::getKappa)
      .def("getR0", &HWPDE::getR0)
      .def("getTimeSigmas", &HWPDE::getTimeSigmas)
      .def("getSigmas", &HWPDE::getSigmas)
      .def("getTimeThetas", &HWPDE::getTimeThetas)
      .def("getThetas", &HWPDE::getThetas);

  // ShortRate1FPDE
  py::class_<ShortRate1FPDE>(m, "ShortRate1FPDE")
      .def(py::init<double, double, double, double, double, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("R0"), py::arg("kappa"), py::arg("alpha"), py::arg("beta"),
           py::arg("gamma"), py::arg("timeSigmas"), py::arg("sigmas"),
           py::arg("timeThetas"), py::arg("thetas"))
      .def("pricingSwaption", &ShortRate1FPDE::pricingSwaption,
           py::arg("Expiry"), py::arg("Tenor"), py::arg("Strike"),
           py::arg("PayFrequency") = 0.5);

  // ShortRate2FPDE
  py::class_<ShortRate2FPDE>(m, "ShortRate2FPDE")
      .def(py::init<double, double, double, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>, std::vector<double>,
                    std::vector<double>>(),
           py::arg("kappa1"), py::arg("kappa2"), py::arg("lambda"),
           py::arg("timeSigma1s"), py::arg("sigma1s"), py::arg("timeSigma2s"),
           py::arg("sigma2s"), py::arg("timeAlphas"), py::arg("alphas"))
      .def("pricingSwaption", &ShortRate2FPDE::pricingSwaption,
           py::arg("Expiry"), py::arg("Tenor"), py::arg("Strike"),
           py::arg("PayFrequency") = 0.5)
      .def("calibrate", &ShortRate2FPDE::calibrator, py::arg("timeDFs"),
           py::arg("DFs"), py::arg("swapQuotes"));

  // Facade Functions (from model.cpp)
  m.def("lv_export", &lvExport, py::arg("maturities"), py::arg("forwards"),
        py::arg("betas"), py::arg("alphas"), py::arg("nus"), py::arg("rhos"),
        py::arg("spot"), py::arg("times"));
}
