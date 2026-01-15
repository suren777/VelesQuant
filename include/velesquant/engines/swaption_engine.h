#ifndef VELESQUANT_ENGINES_SWAPTION_ENGINE_H
#define VELESQUANT_ENGINES_SWAPTION_ENGINE_H

#include <memory>
#include <velesquant/instruments/swaption.h>
#include <velesquant/volatility/sabr.h>

namespace velesquant {
namespace engines {

/**
 * @class SwaptionAnalyticEngine
 * @brief Evaluation engine for Swaption instruments using the SABR model.
 */
class SwaptionAnalyticEngine {
public:
  /**
   * @brief Constructor
   * @param model Pointer to a SABR model. The model should be initialized with
   *              maturity and forward rate corresponding to the Swaption being
   * priced, or updated before calculation.
   */
  explicit SwaptionAnalyticEngine(std::shared_ptr<Sabr> model);

  /**
   * @brief Calculates the fair value of a Swaption.
   * @param instrument The swaption instrument to price.
   * @return The calculated price (premium).
   */
  double calculate(const instruments::Swaption &instrument) const;

  /**
   * @brief Calculates the fair value of the underlying Swap.
   * @param instrument The swaption instrument (containing underlying details).
   * @return The swap value.
   */
  double calculateSwap(const instruments::Swaption &instrument) const;

  /**
   * @brief Calculates implied volatility.
   * @param instrument The swaption instrument.
   * @return The implied volatility.
   */
  double getImpliedVol(const instruments::Swaption &instrument) const;

private:
  std::shared_ptr<Sabr> model_;
};

} // namespace engines
} // namespace velesquant

#endif // VELESQUANT_ENGINES_SWAPTION_ENGINE_H
