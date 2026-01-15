#include <stdexcept>
#include <velesquant/engines/swaption_engine.h>

namespace velesquant {
namespace engines {

SwaptionAnalyticEngine::SwaptionAnalyticEngine(std::shared_ptr<Sabr> model)
    : model_(std::move(model)) {
  if (!model_) {
    throw std::invalid_argument(
        "SwaptionAnalyticEngine: model cannot be null.");
  }
}

double SwaptionAnalyticEngine::calculate(
    const instruments::Swaption &instrument) const {
  // Note: This assumes the Sabr model in model_ is already configured with
  // the correct forward and maturity matching the instrument.
  // In a more advanced design, we might check or update the model here,
  // but the Sabr class combines model params and market state.

  // Check consistency (optional but recommended)
  // if (std::abs(model_->getForward() - instrument.getForward()) > 1e-8) ...

  return instrument.getAnnuity() *
         model_->premiumBlackScholes(instrument.getStrike(),
                                     instrument.getType());
}

double SwaptionAnalyticEngine::calculateSwap(
    const instruments::Swaption &instrument) const {
  return instrument.getAnnuity() *
         (model_->getForward() - instrument.getStrike());
}

double SwaptionAnalyticEngine::getImpliedVol(
    const instruments::Swaption &instrument) const {
  return model_->impliedVol(instrument.getStrike());
}

} // namespace engines
} // namespace velesquant
