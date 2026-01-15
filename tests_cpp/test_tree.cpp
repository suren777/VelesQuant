#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <velesquant/volatility/c_tree.h>
#include <velesquant/models/utility.h> // For OptionType

using namespace velesquant;

// Helper to generate flat curves
std::vector<double> generateFlatCurve(const std::vector<double> &times,
                                      double value) {
  std::vector<double> vals;
  vals.reserve(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    vals.push_back(value);
  }
  return vals;
}

// Helper to calculate Forwards
std::vector<double> calculateForwards(double spot, double r, double q,
                                      const std::vector<double> &times) {
  std::vector<double> fwds;
  fwds.reserve(times.size());
  for (double t : times) {
    fwds.push_back(spot * std::exp((r - q) * t));
  }
  return fwds;
}

class CTreeTestFixture : public ::testing::Test {
protected:
  double spot = 100.0;
  double r = 0.05;
  double q = 0.02;
  double sigma = 0.2;
  double T = 1.0;

  std::vector<double> timeGrid;
  std::vector<double> vols;
  std::vector<double> forwards;
  std::vector<double> rates;
  std::vector<double> divs;

  void SetUp() override {
    // Grid should cover the option maturity T=1.0
    // Let's define points at 0, 0.5, 1.0, 1.5, 2.0
    timeGrid = {0.0, 0.5, 1.0, 1.5, 2.0};
    vols = generateFlatCurve(timeGrid, sigma);
    forwards = calculateForwards(spot, r, q, timeGrid);

    // Pass r and q as flat curves (as CTree currently only uses first element
    // anyway)
    rates = generateFlatCurve(timeGrid, r);
    divs = generateFlatCurve(timeGrid, q);
  }
};

TEST_F(CTreeTestFixture, Construction) {
  CTree tree(spot, timeGrid, forwards, vols, rates, divs);
  SUCCEED();
}

TEST_F(CTreeTestFixture, EuropeanOptionPricing) {
  CTree tree(spot, timeGrid, forwards, vols, rates, divs);
  double strike = 100.0;
  int steps = 50;

  // calculateBinomial(strike, Maturity, Nnodes, style, pay, treeType)

  double priceBinomial = tree.calculateBinomial(
      strike, T, steps, exStyle::European, OptionType::Call, tType::recomb);

  // Analytical Black-Scholes estimate: ~9.23
  EXPECT_NEAR(priceBinomial, 9.23, 0.5);

  double priceTrinomial = tree.calculateTrinomial(
      strike, T, steps, exStyle::European, OptionType::Call, tType::recomb);

  EXPECT_NEAR(priceTrinomial, 9.23, 0.5);
}

TEST_F(CTreeTestFixture, AmericanOptionPricing) {
  CTree tree(spot, timeGrid, forwards, vols, rates, divs);
  double strike = 100.0;
  int steps = 50;

  double putEuro = tree.calculateBinomial(strike, T, steps, exStyle::European,
                                          OptionType::Put, tType::recomb);

  double putAmer = tree.calculateBinomial(strike, T, steps, exStyle::American,
                                          OptionType::Put, tType::recomb);

  // American Put >= European Put
  EXPECT_GE(putAmer, putEuro - 1e-10);
  EXPECT_GT(putAmer, 0.0);
}
