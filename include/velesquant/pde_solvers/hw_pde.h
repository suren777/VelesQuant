//		HWPDE.h

#ifndef HWPDE_H
#define HWPDE_H

#include <list>
#include <vector>
#include <velesquant/models/utility.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
#ifdef _MSC_VER
#pragma warning(disable : 4018)
#endif

class HWPDE {
public:
  HWPDE(double R0, double kappa, std::vector<double> timeSigmas,
        std::vector<double> sigmas, std::vector<double> timeThetas,
        std::vector<double> thetas)
      : R0_(R0), kappa_(kappa), timeSigmas_(timeSigmas), sigmas_(sigmas),
        timeThetas_(timeThetas), thetas_(thetas) {
    buildGrid();
    sigmas0_ = sigmas_;
    kappa0_ = kappa_;
  };

  HWPDE(double kappa, std::vector<double> timeSigmas,
        std::vector<double> sigmas, std::vector<double> timeDF,
        std::vector<double> DF)
      : kappa_(kappa), timeSigmas_(timeSigmas), sigmas_(sigmas),
        timeDFs_(timeDF), DFs_(DF) {
    timeThetas_.push_back(0.5);
    double dt = (timeDFs_[timeDFs_.size() - 1] - timeThetas_[0]) / 20;
    for (double i = timeThetas_[0] + dt; i <= timeDFs_[timeDFs_.size() - 1];
         i += dt)
      timeThetas_.push_back(i);
    R0_ = -log(DFs_[1]) / timeDFs_[1];
    buildGrid();
    for (size_t i = 0; i < timeThetas_.size(); i++)
      thetas_.push_back(.001);
    termStructureCalibrator();
    sigmas0_ = sigmas_;
    kappa0_ = kappa_;
  };

  HWPDE(std::vector<double> timeDFs, std::vector<double> DFs,
        std::vector<defSwap> quoteSwap)
      : timeDFs_(timeDFs), DFs_(DFs), quoteSwap_(quoteSwap) {
    R0_ = -log(DFs_[0]) / timeDFs_[0];
    kappa_ = 0.03;
    timeThetas_ = timeDFs_;
    thetas_.resize(timeThetas_.size(), 0.005);
    std::list<double> Times;
    for (size_t i = 0; i < quoteSwap_.size(); i++) {
      double expiry = quoteSwap_[i].Expiry;
      Times.push_back(expiry);
    };
    Times.sort();
    Times.unique(); // sorting in order and removing duplicates
    for (std::list<double>::iterator it = Times.begin(); it != Times.end();
         ++it)
      timeSigmas_.push_back(*it);
    sigmas_.resize(timeSigmas_.size(), 0.003);
    sigmas0_ = sigmas_;
    kappa0_ = kappa_;
    buildGrid();
    calibrator(timeDFs_, DFs_, quoteSwap_);
  };

  ~HWPDE() {};
  double get_cal_time() { return cal_time_; }
  void calibrator(std::vector<double> timeDFs, std::vector<double> DFs,
                  std::vector<defSwap> swapQuotes);
  void calibratorBootStrap(std::vector<double> timeDFs, std::vector<double> DFs,
                           std::vector<defSwap> swapQuotes);

  void discountBack(double Expiry, double Maturity, std::vector<double> &f);
  double pricingZB(double Maturity);
  double pricingZBO(double Expiry, double Maturity, double Strike,
                    const std::string &type);
  double pricingCouponBond(double Expiry, double Tenor, double Coupon,
                           double PayFrequency);
  void pricingCouponBondt(double Expiry, double Tenor, double Coupon,
                          double PayFrequency, std::vector<double> &f);
  double pricingCBO(double Expiry, double Tenor, double Coupon, double Strike,
                    double PayFrequency, const std::string &type);
  double pricingSwap(double Expiry, double Tenor, double Strike,
                     double PayFrequency);
  double pricingCallableSwap(double Expiry, double Tenor,
                             std::vector<double> Exercises, double Coupon,
                             double Strike, double PayFrequency,
                             const std::string &type);
  double pricingSwaption(double Expiry, double Tenor, double Strike,
                         double PayFrequency);
  double pricingBermudan(double Expiry, double Tenor,
                         std::vector<double> Exercises, double Strike,
                         double PayFrequency);

  std::vector<double> getDFs(std::vector<double> &timePoints);
  double getImpVolATM(double Expiry, double Tenor, double PayFrequency);
  double getSwapRate(double Expiry, double Tenor, double PayFrequency);

  double getDGinterp(double t);

  std::vector<double> simulationPDE(std::vector<double> timePoints) const;

  double getKappa() { return kappa_; };
  double getR0() { return R0_; };
  std::vector<double> getTimeSigmas() { return timeSigmas_; };
  std::vector<double> getSigmas() { return sigmas_; };
  std::vector<double> getTimeThetas() { return timeThetas_; };
  std::vector<double> getThetas() { return thetas_; };
  std::vector<double> parallel_check();

private:
  double R0_, kappa_, kappa0_, cal_time_;
  double vleft_, vright_;
  mutable std::vector<double> timeSigmas_, sigmas_, sigmas0_;
  mutable std::vector<double> timeThetas_, thetas_; // time_[0] first point
  mutable std::vector<double> a_, b_, c_, e_, f_, g_;
  std::vector<double> timeDFs_, DFs_;
  std::vector<defSwap> quoteSwap_;
  std::vector<defSwap> qSB_;
  int iter_, iR0_;

  mutable std::vector<double> marketSwaption_;
  mutable std::vector<double> gridR_;
  void buildGrid(double Rmax = 0.5, double factor = 0.7);
  void oneStepBackward(const int t, std::vector<double> &inV);
  void oneStepForward(const int t, std::vector<double> &inV);

  void termStructureCalibrator();
  void termStructureCalibratorBtstrp(double);
  void objFcnCalibration(int m, int n, double *x, double *fvec, int *iflag,
                         double *lb, double *ub);
  double objFcnCalibrationB(double x, int m);
  double pen_fun(double *x, double *lb, double *ub, int n);

  double getAnnuity(double Expiry, double Tenor, double PayFrequency);
  double swaptionValueATM(double Expiry, double Tenor, double PayFrequency,
                          double SwapRate, double VolATM);
  double trapezoidal(std::vector<double> &inV);
  double swaptionIVblack(double Expiry, double Tenor, double swaption_price);
  double swaptionATM(double Expiry, double Tenor, double VolATM) {
    return (getDFinterp(Expiry) - getDFinterp(Expiry + Tenor)) *
           (2.0 * cdf_normal(0.5 * VolATM * sqrt(Expiry)) - 1.0);
  };
  double whichSigma(double t) const;
  double whichTheta(double t) const;
  double Bratio(double a, double Mi, double Tj, double Tk);
  double dBratio(double a, double Mi, double Tk, double Tj);
  void calibrateKappa();
  double objKappa(double a, const std::vector<double> &ivr,
                  const std::vector<double> &br, const std::vector<int> &ind);
  double dobjKappa(double a, const std::vector<double> &ivr,
                   const std::vector<double> &br, const std::vector<int> &ind);
  // Mesh initialisation
  std::vector<double> grMesh(int, double Rmax = .5, double factor = 0.7);
  double getDFinterp(double t);
};

} // namespace velesquant
#endif