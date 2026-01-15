#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace velesquant {

// Primary types using Eigen
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using RowVector = Eigen::RowVectorXd;

namespace xlw {

// CellMatrix - a wrapper around Eigen::MatrixXd that provides backward
// compatibility with the legacy xlw API
class CellMatrix {
public:
  CellMatrix() : data_() {}

  CellMatrix(size_t rows, size_t cols)
      : data_(Eigen::MatrixXd::Zero(rows, cols)) {}

  CellMatrix(double v) : data_(Eigen::MatrixXd::Constant(1, 1, v)) {}

  CellMatrix(const std::vector<double> &v)
      : data_(Eigen::Map<const Eigen::VectorXd>(v.data(), v.size())) {}

  CellMatrix(const std::string & /*s*/) : data_(Eigen::MatrixXd::Zero(1, 1)) {}

  CellMatrix(const Eigen::MatrixXd &m) : data_(m) {}

  // Proxy class for element access that mimics the old Cell behavior
  class CellProxy {
  public:
    CellProxy(double &ref) : ref_(ref) {}

    double NumericValue() const { return ref_; }
    bool IsANumber() const { return true; }
    operator double() const { return ref_; }

    CellProxy &operator=(double v) {
      ref_ = v;
      return *this;
    }
    CellProxy &operator+=(double v) {
      ref_ += v;
      return *this;
    }

  private:
    double &ref_;
  };

  class ConstCellProxy {
  public:
    ConstCellProxy(const double &ref) : ref_(ref) {}

    double NumericValue() const { return ref_; }
    bool IsANumber() const { return true; }
    operator double() const { return ref_; }

  private:
    const double &ref_;
  };

  CellProxy operator()(size_t i, size_t j) {
    // Auto-resize if needed (legacy behavior)
    if (static_cast<Eigen::Index>(i) >= data_.rows() ||
        static_cast<Eigen::Index>(j) >= data_.cols()) {
      Eigen::Index newRows =
          std::max(data_.rows(), static_cast<Eigen::Index>(i + 1));
      Eigen::Index newCols =
          std::max(data_.cols(), static_cast<Eigen::Index>(j + 1));
      data_.conservativeResize(newRows, newCols);
    }
    return CellProxy(data_(i, j));
  }

  ConstCellProxy operator()(size_t i, size_t j) const {
    return ConstCellProxy(data_(i, j));
  }

  // Legacy API compatibility
  size_t RowsInStructure() const { return data_.rows(); }
  size_t ColumnsInStructure() const { return data_.cols(); }

  // Modern API aliases
  size_t rows() const { return data_.rows(); }
  size_t cols() const { return data_.cols(); }

  // Access to underlying Eigen matrix
  Eigen::MatrixXd &eigen() { return data_; }
  const Eigen::MatrixXd &eigen() const { return data_; }

private:
  Eigen::MatrixXd data_;
};

using MyMatrix = CellMatrix;
using MyArray = std::vector<double>; // Keep as std::vector for compatibility

} // namespace xlw

// Bring xlw types into velesquant namespace for convenience
using xlw::CellMatrix;
using xlw::MyArray;
using xlw::MyMatrix;

} // namespace velesquant
