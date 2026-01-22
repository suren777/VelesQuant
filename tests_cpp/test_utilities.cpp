#include <cmath>
#include <gtest/gtest.h>
#include <vector>
#include <velesquant/models/utility.h>
#include <velesquant/numerics/interpolation.h>

using namespace velesquant;

TEST(UtilityTest, CumulativeNormal) {
  EXPECT_NEAR(cdf_normal(0.0), 0.5, 1e-12);
  EXPECT_NEAR(cdf_normal(1.96), 0.975, 1e-3);
  EXPECT_NEAR(cdf_normal(-1.96), 0.025, 1e-3);
}

TEST(UtilityTest, Interpolation) {
  std::vector<double> X = {0.0, 1.0, 2.0};
  std::vector<double> Y = {0.0, 1.0, 4.0}; // x^2

  // Test linear interpolation
  double val = interpolation("Linear", X, Y, 0.5);
  EXPECT_DOUBLE_EQ(val, 0.5); // (0+1)/2 = 0.5

  val = interpolation("Linear", X, Y, 1.5);
  EXPECT_DOUBLE_EQ(val, 2.5); // (1+4)/2 = 2.5
}

TEST(UtilityTest, PdfNormal) {
  // PDF of standard normal at x=0 should be 1/sqrt(2*pi) ≈ 0.3989
  EXPECT_NEAR(pdf_normal(0.0), 0.3989422804, 1e-6);
  EXPECT_GT(pdf_normal(1.0), 0.0);
  EXPECT_LT(pdf_normal(1.0), pdf_normal(0.0)); // PDF decreases away from mean
}

TEST(UtilityTest, OptionVega) {
  double F = 100.0;
  double K = 100.0;
  double sigma = 0.2;
  double T = 1.0;

  double vega = option_vega(F, K, sigma, T);

  EXPECT_GT(vega, 0.0); // Vega should be positive
  EXPECT_LT(vega, F);   // Vega should be less than forward price
}

TEST(UtilityTest, DFtoR) {
  std::vector<double> DF = {0.9512, 0.9048}; // exp(-0.05*1), exp(-0.05*2)
  std::vector<double> T = {1.0, 2.0};

  std::vector<double> rates = DFtoR(DF, T, 0); // Continuous compounding

  EXPECT_NEAR(rates[0], 0.05, 0.01);
  EXPECT_NEAR(rates[1], 0.05, 0.01);
}

TEST(UtilityTest, FwdPrice) {
  double spot = 100.0;
  double r = 0.05;
  double div = 0.02;
  double T = 1.0;

  double fwd = FwdPrice(spot, r, div, T);
  double expected = spot * std::exp((r - div) * T);

  EXPECT_DOUBLE_EQ(fwd, expected);
}

TEST(UtilityTest, Cholesky) {
  // Test 2x2 identity matrix (Cholesky should be identity)
  std::vector<std::vector<double>> corrMatrix = {{1.0, 0.0}, {0.0, 1.0}};

  std::vector<std::vector<double>> chol = cholesky(corrMatrix);

  EXPECT_NEAR(chol[0][0], 1.0, 1e-10);
  EXPECT_NEAR(chol[1][1], 1.0, 1e-10);
  EXPECT_NEAR(chol[1][0], 0.0, 1e-10);
}
