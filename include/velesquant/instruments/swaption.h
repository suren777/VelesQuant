#ifndef VELESQUANT_INSTRUMENTS_SWAPTION_H
#define VELESQUANT_INSTRUMENTS_SWAPTION_H

#include <velesquant/instruments/instrument.h>
#include <velesquant/models/utility.h> // For OptionType

namespace velesquant {
namespace instruments {

/**
 * @class Swaption
 * @brief Represents a Swaption contract (option on a swap).
 * Contains only contract details, no pricing model or calibration logic.
 */
class Swaption : public Instrument {
public:
  Swaption(double expiry, double tenor, double forward, double annuity,
           double strike, OptionType type = OptionType::Call)
      : expiry_(expiry), tenor_(tenor), forward_(forward), annuity_(annuity),
        strike_(strike), type_(type) {}

  double getExpiry() const { return expiry_; }
  double getTenor() const { return tenor_; }
  double getForward() const { return forward_; }
  double getAnnuity() const { return annuity_; }
  double getStrike() const { return strike_; }
  OptionType getType() const { return type_; }

private:
  double expiry_;
  double tenor_;
  double forward_;
  double annuity_;
  double strike_;
  OptionType type_;
};

} // namespace instruments
} // namespace velesquant

#endif // VELESQUANT_INSTRUMENTS_SWAPTION_H
