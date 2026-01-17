#ifndef VELESQUANT_ENGINES_SWAPTION_ENGINE_H
#define VELESQUANT_ENGINES_SWAPTION_ENGINE_H

#include <memory>
#include <stdexcept>
#include <velesquant/instruments/swaption.h>
#include <velesquant/models/concepts.h>
// Note: We no longer include sabr.h here to avoid loose coupling.
// The ModelType template parameter will handle the dependencies via duck
// typing.

namespace velesquant {
namespace engines {

/**
 * @class SwaptionAnalyticEngine
 * @brief Generic evaluation engine for Swaption instruments.
 *
 * This class is templated on the ModelType, allowing it to work with any model
 * that satisfies the PricingModel concept.
 */
template <velesquant::models::PricingModel ModelType>
class SwaptionAnalyticEngine {
public:
  /**
   * @brief Constructor
   * @param model Pointer to the pricing model (e.g., Sabr).
   */
  explicit SwaptionAnalyticEngine(std::shared_ptr<ModelType> model)
      : model_(std::move(model)) {
    if (!model_) {
      throw std::invalid_argument(
          "SwaptionAnalyticEngine: model cannot be null.");
    }
  }

  /**
   * @brief Calculates the fair value of a Swaption.
   * @param instrument The swaption instrument to price.
   * @return The calculated price (premium).
   */
  double calculate(const instruments::Swaption &instrument) const {
    // ModelType is expected to have premiumBlackScholes()
    return instrument.getAnnuity() *
           model_->premiumBlackScholes(instrument.getStrike(),
                                       instrument.getType());
  }

  /**
   * @brief Calculates the fair value of the underlying Swap.
   * @param instrument The swaption instrument (containing underlying details).
   * @return The swap value.
   */
  double calculateSwap(const instruments::Swaption &instrument) const {
    // ModelType is expected to have getForward()
    return instrument.getAnnuity() *
           (model_->getForward() - instrument.getStrike());
  }

  /**
   * @brief Calculates implied volatility.
   * @param instrument The swaption instrument.
   * @return The implied volatility.
   */
  double getImpliedVol(const instruments::Swaption &instrument) const {
    // ModelType is expected to have impliedVol()
    return model_->impliedVol(instrument.getStrike());
  }

private:
  std::shared_ptr<ModelType> model_;
};

} // namespace engines
} // namespace velesquant

#endif // VELESQUANT_ENGINES_SWAPTION_ENGINE_H
