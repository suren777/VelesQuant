#include <velesquant/models/swaption.h>

#include <algorithm>
#include <string>
#include <vector>
#include <velesquant/engines/swaption_engine.h>
#include <velesquant/instruments/swaption.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>
#include <velesquant/volatility/sabr.h>

// using namespace std;
// using namespace velesquant::xlw; // Not needed if types.h used correctly

namespace velesquant {

swaption::swaption(double expiry, double tenor, double forward, double annuity,
                   double beta, double alpha, double nu, double rho) {
  // Instrument initialized with forward as the strike since the constructor
  // does not provide a specific strike.
  instrument_ = std::make_shared<instruments::Swaption>(expiry, tenor, forward,
                                                        annuity, forward);
  model_ = std::make_shared<Sabr>(expiry, forward, beta, alpha, nu, rho);
  engine_ = std::make_shared<engines::SwaptionAnalyticEngine<Sabr>>(model_);
};

swaption::swaption(double expiry, double tenor, double forward, double annuity,
                   double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
                   CalibrationTarget quoteType) {
  instrument_ = std::make_shared<instruments::Swaption>(expiry, tenor, forward,
                                                        annuity, forward);
  model_ = std::make_shared<Sabr>(expiry, forward, beta);

  // Calibration logic from old constructor
  int m = strikes.RowsInStructure();
  std::vector<double> theStrikes(m), theQuotes(m);
  for (int i = 0; i < m; ++i) {
    theStrikes[i] = strikes(i, 0).NumericValue();
    theQuotes[i] = marketQuotes(i, 0).NumericValue();
  }
  model_->calibrator(theStrikes, theQuotes, quoteType);

  engine_ = std::make_shared<engines::SwaptionAnalyticEngine<Sabr>>(model_);
};

swaption::swaption(double expiry, double tenor, double forward, double annuity,
                   double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
                   CalibrationTarget quoteType,
                   [[maybe_unused]] CellMatrix &initialParams) {
  instrument_ = std::make_shared<instruments::Swaption>(expiry, tenor, forward,
                                                        annuity, forward);
  model_ = std::make_shared<Sabr>(expiry, forward, beta);

  int m = strikes.RowsInStructure();
  std::vector<double> theStrikes(m), theQuotes(m);
  for (int i = 0; i < m; ++i) {
    theStrikes[i] = strikes(i, 0).NumericValue();
    theQuotes[i] = marketQuotes(i, 0).NumericValue();
  }
  // Legacy code handled 'initialParams' variable but didn't use it in the
  // final calibrator call. We replicate that behavior here.
  model_->calibratorWithInitial(theStrikes, theQuotes, quoteType);
  engine_ = std::make_shared<engines::SwaptionAnalyticEngine<Sabr>>(model_);
};

double swaption::swaptionFairValue(double strike, OptionType callORput) const {
  // Create a temporary instrument with the specific strike and option type
  // for this calculation, reusing other properties from the stored instrument.
  instruments::Swaption tempInstr(
      instrument_->getExpiry(), instrument_->getTenor(),
      instrument_->getForward(), instrument_->getAnnuity(), strike, callORput);
  return engine_->calculate(tempInstr);
};

double swaption::swapFairValue(double strike) const {
  instruments::Swaption tempInstr(
      instrument_->getExpiry(), instrument_->getTenor(),
      instrument_->getForward(), instrument_->getAnnuity(), strike);
  return engine_->calculateSwap(tempInstr);
};

double swaption::simulation(double corrRN) {
  return model_->simulation(corrRN);
};

double swaption::getImpliedVol(double strike) const {
  instruments::Swaption tempInstr(
      instrument_->getExpiry(), instrument_->getTenor(),
      instrument_->getForward(), instrument_->getAnnuity(), strike);
  return engine_->getImpliedVol(tempInstr);
};

double swaption::getParameterAlpha() const {
  return model_->getParameterAlpha();
};

double swaption::getParameterNu() const { return model_->getParameterNu(); };

double swaption::getParameterRho() const { return model_->getParameterRho(); };
} // namespace velesquant