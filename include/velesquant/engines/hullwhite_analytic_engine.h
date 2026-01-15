#ifndef VELESQUANT_ENGINES_HULLWHITE_ANALYTIC_ENGINE_H
#define VELESQUANT_ENGINES_HULLWHITE_ANALYTIC_ENGINE_H

#include <memory>
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/utility.h>

namespace velesquant {
namespace engines {

class HullWhiteAnalyticEngine {
public:
  explicit HullWhiteAnalyticEngine(
      std::shared_ptr<models::HullWhiteModel> model);

  [[nodiscard]] double optionBond(double expiry, double maturity, double strike,
                                  OptionType type) const;

  [[nodiscard]] double swaption(double expiry, double tenor, double strike,
                                double payFrequency = 0.5) const;

  // Helper exposed for testing or other uses
  [[nodiscard]] double criticalPoint(double expiry, double tenor, double strike,
                                     double payFrequency) const;

private:
  std::shared_ptr<models::HullWhiteModel> model_;

  // Private helpers
  double formulaBlack(double dfT0, double dfTN, double strike, double impVol,
                      double T0, OptionType type) const;
  double equalityCP(double x, double v0, double expiry, double tenor,
                    double strike, double payFrequency) const;
};

} // namespace engines
} // namespace velesquant

#endif // VELESQUANT_ENGINES_HULLWHITE_ANALYTIC_ENGINE_H
