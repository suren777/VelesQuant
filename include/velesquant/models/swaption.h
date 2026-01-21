//		swaption.h

#ifndef SWAPTION_H
#define SWAPTION_H

#include <memory>
#include <string>
#include <velesquant/engines/swaption_engine.h>
#include <velesquant/instruments/swaption.h>
#include <velesquant/models/utility.h>
#include <velesquant/types.h>
#include <velesquant/volatility/sabr.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor

class swaption {
public:
  swaption(double expiry, double tenor, double forward, double annuity,
           double beta = 0.85, double alpha = 0.5, double nu = 0.25,
           double rho = -0.75);

  swaption(double expiry, double tenor, double forward, double annuity,
           double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
           CalibrationTarget quoteType = CalibrationTarget::Price);

  swaption(double expiry, double tenor, double forward, double annuity,
           double beta, CellMatrix &strikes, CellMatrix &marketQuotes,
           CalibrationTarget quoteType, CellMatrix &initialParams);

  ~swaption() = default;

  double swaptionFairValue(double strike,
                           OptionType callORput = OptionType::Call) const;
  double swapFairValue(double strike) const;

  double simulation(double corrRN);

  double getImpliedVol(double strike) const;

  double getParameterAlpha() const;
  double getParameterNu() const;
  double getParameterRho() const;

private:
  std::shared_ptr<instruments::Swaption> instrument_;
  std::shared_ptr<Sabr> model_;
  std::shared_ptr<engines::SwaptionAnalyticEngine<Sabr>> engine_;
};

} // namespace velesquant
#endif