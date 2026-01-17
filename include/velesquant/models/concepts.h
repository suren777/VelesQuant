#ifndef VELESQUANT_MODELS_CONCEPTS_H
#define VELESQUANT_MODELS_CONCEPTS_H

#include <concepts>
#include <velesquant/models/utility.h>

namespace velesquant {
namespace models {

/**
 * @concept PricingModel
 * @brief Concept defining the requirements for a model used in pricing engines.
 *
 * A PricingModel must provide:
 * - getForward(): returns the forward rate.
 * - premiumBlackScholes(strike, type): returns option premium.
 * - impliedVol(strike): returns implied volatility.
 */
template <typename T>
concept PricingModel = requires(T m, double k, OptionType type) {
  { m.getForward() } -> std::convertible_to<double>;
  { m.premiumBlackScholes(k, type) } -> std::convertible_to<double>;
  { m.impliedVol(k) } -> std::convertible_to<double>;
};

/**
 * @concept ShortRateModel
 * @brief Concept defining requirements for a short-rate model (e.g.
 * Hull-White).
 */
template <typename T>
concept ShortRateModel = requires(T m, double t, double T_end) {
  { m.getKappa() } -> std::convertible_to<double>;
  { m.variance(t, T_end) } -> std::convertible_to<double>;
  { m.getDiscountFactor(t) } -> std::convertible_to<double>;
};

} // namespace models
} // namespace velesquant

#endif // VELESQUANT_MODELS_CONCEPTS_H
