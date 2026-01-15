#include <algorithm>
#include <boost/math/distributions.hpp>
#include <cmath>
#include <functional>
#include <memory>
#include <ql/quantlib.hpp>
#include <vector>

#include <velesquant/engines/hullwhite_analytic_engine.h>
#include <velesquant/models/hullwhite_calibrator.h>
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/hw.h>
#include <velesquant/models/utility.h>
// #include <velesquant/volatility/lm.h>

#include <velesquant/constants.h>

const double DT = velesquant::constants::numerical::HULL_WHITE_DT;
const double SDT = std::sqrt(DT);