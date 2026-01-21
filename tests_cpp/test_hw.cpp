#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <vector>
#include <velesquant/engines/hullwhite_analytic_engine.h>
#include <velesquant/models/hullwhite_model.h>
#include <velesquant/models/utility.h> // for annuity

using namespace velesquant;

// Helper to generate flat rate discount factors
namespace {
std::vector<double> generateDFs(const std::vector<double> &times, double rate) {
  std::vector<double> dfs;
  dfs.reserve(times.size());
  for (double t : times) {
    dfs.push_back(std::exp(-rate * t));
  }
  return dfs;
}

// Calculate swap rate using model
double calculateSwapRate(const models::HullWhiteModel &model, double expiry,
                         double tenor, double payFreq) {
  double startTime = expiry;
  double endTime = expiry + tenor;

  double pStart = model.getDiscountFactor(startTime);
  double pEnd = model.getDiscountFactor(endTime);

  double annuity = 0.0;
  int nPeriods = static_cast<int>(std::round(tenor / payFreq));

  for (int i = 1; i <= nPeriods; ++i) {
    double t = startTime + i * payFreq;
    annuity += payFreq * model.getDiscountFactor(t);
  }

  if (annuity < 1e-12)
    return 0.0;
  return (pStart - pEnd) / annuity;
}

} // namespace

class HullWhiteTestFixture : public ::testing::Test {
protected:
  double r0 = 0.05;
  double kappa = 0.1;
  double sigma = 0.01;
  std::vector<double> timeSigmas = {0.0, 30.0};
  std::vector<double> sigmas = {sigma, sigma};
  std::vector<double> timeDFs;
  std::vector<double> DFs;

  void SetUp() override {
    // Setup times 0 to 10 years
    for (int i = 0; i <= 10; ++i) {
      timeDFs.push_back(static_cast<double>(i));
    }
    DFs = generateDFs(timeDFs, r0);
  }
};

TEST_F(HullWhiteTestFixture, Construction) {
  auto model = std::make_shared<models::HullWhiteModel>(kappa, timeSigmas,
                                                        sigmas, timeDFs, DFs);
  EXPECT_NEAR(model->getKappa(), kappa, 1e-12);
  auto gotSigmas = model->getSigmas();
  EXPECT_EQ(gotSigmas.size(), 2);
  EXPECT_DOUBLE_EQ(gotSigmas[0], sigma);
}

TEST_F(HullWhiteTestFixture, ZCBondPricing) {
  auto model = std::make_shared<models::HullWhiteModel>(kappa, timeSigmas,
                                                        sigmas, timeDFs, DFs);

  // Price a bond maturing in 1 year
  // For T=1, ZC bond price should be exactly DF(1) = exp(-0.05) ~ 0.9512
  double price = model->getDiscountFactor(1.0);
  double expected = DFs[1]; // timeDFs[1] is 1.0

  // The ZC(T) should be exactly the discount factor from the curve
  EXPECT_NEAR(price, expected, 1e-12);
}

TEST_F(HullWhiteTestFixture, OptionOnBond) {
  auto model = std::make_shared<models::HullWhiteModel>(kappa, timeSigmas,
                                                        sigmas, timeDFs, DFs);
  auto engine = std::make_shared<
      engines::HullWhiteAnalyticEngine<models::HullWhiteModel>>(model);

  double expiry = 1.0;
  double maturity = 2.0;

  // Strike = Forward Price = P(0,2)/P(0,1)
  double forwardBondPrice = DFs[2] / DFs[1];
  double strike = forwardBondPrice;

  double callPrice =
      engine->optionBond(expiry, maturity, strike, OptionType::Call);
  double putPrice =
      engine->optionBond(expiry, maturity, strike, OptionType::Put);

  // Put-Call Parity: C - P = P(0, M) - K * P(0, E)
  double parityTarget = DFs[2] - strike * DFs[1];

  // Since strike = DFs[2]/DFs[1], parityTarget should be 0.
  double parityDiff = callPrice - putPrice;

  EXPECT_NEAR(parityDiff, parityTarget, 1e-10);
}

TEST_F(HullWhiteTestFixture, SwaptionPricing) {
  auto model = std::make_shared<models::HullWhiteModel>(kappa, timeSigmas,
                                                        sigmas, timeDFs, DFs);
  auto engine = std::make_shared<
      engines::HullWhiteAnalyticEngine<models::HullWhiteModel>>(model);

  double expiry = 1.0;
  double tenor = 2.0;
  double payFreq = 0.5; // Semi-annual

  // Get forward swap rate
  double swapRate = calculateSwapRate(*model, expiry, tenor, payFreq);
  EXPECT_GT(swapRate, 0.0);
  EXPECT_LT(swapRate, 0.2); // Reasonable swap rate

  // Price ATM swaption
  double strike = swapRate;
  double price = engine->swaption(expiry, tenor, strike, payFreq);

  EXPECT_GT(price, 0.0);
  EXPECT_LT(price, 1.0); // Swaption price should be reasonable
}

TEST_F(HullWhiteTestFixture, Simulation) {
  auto model = std::make_shared<models::HullWhiteModel>(kappa, timeSigmas,
                                                        sigmas, timeDFs, DFs);

  std::vector<double> times = {0.25, 0.5, 0.75, 1.0};
  std::vector<double> rates = model->simulation(times);

  EXPECT_EQ(rates.size(), times.size());

  // All simulated rates should be finite
  for (double r : rates) {
    EXPECT_TRUE(std::isfinite(r));
  }
}
