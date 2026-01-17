#ifndef VELESQUANT_MODELS_HULLWHITE_CALIBRATOR_H
#define VELESQUANT_MODELS_HULLWHITE_CALIBRATOR_H

#include <memory>
#include <vector>
#include <velesquant/engines/hullwhite_analytic_engine.h>
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/utility.h>

namespace velesquant {
namespace models {

class HullWhiteCalibrator {
public:
  HullWhiteCalibrator(
      std::shared_ptr<HullWhiteModel> model,
      std::shared_ptr<engines::HullWhiteAnalyticEngine<HullWhiteModel>> engine);

  void calibrate(const std::vector<defSwap> &swapQuotes,
                 CalibrationTarget target);
  void calibrateBootstrap(const std::vector<defSwap> &swapQuotes,
                          CalibrationTarget target);

private:
  std::shared_ptr<HullWhiteModel> model_;
  std::shared_ptr<engines::HullWhiteAnalyticEngine<HullWhiteModel>> engine_;

  // Internal state for optimization callbacks
  std::vector<defSwap> quoteSwap_;
  std::vector<double> marketSwaption_;

  // Objective functions
  void objFcnPrice(int m, int n, double *x, double *fvec, int *iflag,
                   double *lb, double *ub);
  void objFcnIV(int m, int n, double *x, double *fvec, int *iflag, double *lb,
                double *ub);
  double pen_fun(double *x, double *lb, double *ub, int n);

  // Helper for bootstrap
  void calibrateKappa();
  double objKappa(double x, const std::vector<double> &ivr,
                  const std::vector<double> &br, const std::vector<int> &ind);

  // Helpers refactored from HW
  double Bratio(double a, double M, double t1, double t2);
};

} // namespace models
} // namespace velesquant

#endif // VELESQUANT_MODELS_HULLWHITE_CALIBRATOR_H
