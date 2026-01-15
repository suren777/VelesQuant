#ifndef VELESQUANT_CONSTANTS_H
#define VELESQUANT_CONSTANTS_H

namespace velesquant {
namespace constants {

namespace calibration {
constexpr double TOLERANCE = 1e-10;
constexpr double EPSILON = 1e-10;
constexpr int MAX_ITERATIONS = 5000;
} // namespace calibration

namespace numerical {
constexpr double TIME_STEP = 0.001;
constexpr double HULL_WHITE_DT = 0.001;
} // namespace numerical

namespace sabr {
// SABR formula coefficients (Hagan et al.)
constexpr double COEFF_24 = 1.0 / 24.0;
constexpr double COEFF_1920 = 1.0 / 1920.0;
constexpr double COEFF_4 = 1.0 / 4.0;
constexpr double COEFF_12 = 1.0 / 12.0;

// Bump factors
constexpr double BUMP_UP = 1.001;
constexpr double BUMP_DOWN = 0.999;
constexpr double ATM_BUMP_UP = 1.00001;
constexpr double ATM_BUMP_DOWN = 0.99999;
constexpr double SPOT_BUMP_UP = 1.001;
constexpr double SPOT_BUMP_DOWN = 0.999;
} // namespace sabr

namespace calibration {
constexpr double PENALTY_HUGE = 1e32;
constexpr double PENALTY_LARGE = 1e8;
constexpr double PENALTY_LAMBDA = 1e6;
constexpr double TOL_COARSE = 1e-8;
constexpr double TOL_FINE = 1e-12; // 1e-12
constexpr int MAX_ITER_LOW = 20;
constexpr int MAX_ITER_HIGH = 100;
} // namespace calibration

namespace numerical {
constexpr double FD_BUMP_SIZE = 0.001;
} // namespace numerical

} // namespace constants
} // namespace velesquant

#endif // VELESQUANT_CONSTANTS_H
