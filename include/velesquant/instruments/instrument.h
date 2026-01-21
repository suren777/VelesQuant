#ifndef VELESQUANT_INSTRUMENTS_INSTRUMENT_H
#define VELESQUANT_INSTRUMENTS_INSTRUMENT_H

namespace velesquant {
namespace instruments {

/**
 * @class Instrument
 * @brief Base class for all financial instruments.
 */
class Instrument {
public:
  virtual ~Instrument() = default;
};

} // namespace instruments
} // namespace velesquant

#endif // VELESQUANT_INSTRUMENTS_INSTRUMENT_H
