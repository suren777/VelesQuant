
namespace velesquant {

namespace Intepolators {

class LogLinear {
public:
  template <class I1, class I2>
  Interpolation interpolate(const I1 &xBegin, const I1 &xEnd,
                            const I2 &yBegin) const {
    return LogLinearInterpolation(xBegin, xEnd, yBegin);
  }
  static const int requiredPoints = 2;
};

} // namespace Intepolators
} // namespace velesquant