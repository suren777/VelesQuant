#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

namespace velesquant {

class VelesException : public std::runtime_error {
public:
  explicit VelesException(const std::string &message)
      : std::runtime_error(message) {}
};

// Helper macro to throw exception with file and line info
#define VEL_RAISE(msg)                                                         \
  do {                                                                         \
    std::ostringstream oss;                                                    \
    oss << msg << " [" << __FILE__ << ":" << __LINE__ << "]";                  \
    throw ::velesquant::VelesException(oss.str());                             \
  } while (0)

// Helper macro to check condition and throw if false
#define VEL_CHECK(cond, msg)                                                   \
  do {                                                                         \
    if (!(cond)) {                                                             \
      VEL_RAISE(msg);                                                          \
    }                                                                          \
  } while (0)

} // namespace velesquant
