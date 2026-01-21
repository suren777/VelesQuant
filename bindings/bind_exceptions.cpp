// bindings/bind_exceptions.cpp - Custom exception hierarchy for quant errors
#include "bind_common.h"
#include <stdexcept>

namespace velesquant {

// Custom exception classes
class VelesQuantError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class CalibrationError : public VelesQuantError {
public:
  using VelesQuantError::VelesQuantError;
};

class PricingError : public VelesQuantError {
public:
  using VelesQuantError::VelesQuantError;
};

class NumericalError : public VelesQuantError {
public:
  using VelesQuantError::VelesQuantError;
};

namespace bindings {

void bind_exceptions(py::module_ &m) {
  // Register base exception
  static py::exception<VelesQuantError> velesquant_error(m, "VelesQuantError");

  // Register derived exceptions with inheritance
  static py::exception<CalibrationError> calibration_error(
      m, "CalibrationError", velesquant_error.ptr());
  static py::exception<PricingError> pricing_error(m, "PricingError",
                                                   velesquant_error.ptr());
  static py::exception<NumericalError> numerical_error(m, "NumericalError",
                                                       velesquant_error.ptr());

  // Register exception translators using py::set_error (pybind11 3.x)
  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p)
        std::rethrow_exception(p);
    } catch (const CalibrationError &e) {
      py::set_error(calibration_error, e.what());
    } catch (const PricingError &e) {
      py::set_error(pricing_error, e.what());
    } catch (const NumericalError &e) {
      py::set_error(numerical_error, e.what());
    } catch (const VelesQuantError &e) {
      py::set_error(velesquant_error, e.what());
    }
  });
}

} // namespace bindings
} // namespace velesquant
