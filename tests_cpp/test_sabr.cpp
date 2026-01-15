#include <gtest/gtest.h>
#include <iostream>
#include <velesquant/volatility/sabr.h>

using namespace velesquant;

TEST(SabrTest, Construction) {
  // Basic construction check
  // sabr(double maturity, double forward, double beta, double alpha, double nu,
  // double rho)
  double maturity = 1.0;
  double forward = 0.05;
  double beta = 0.5;
  double alpha = 0.2;
  double nu = 0.4;
  double rho = -0.3;

  Sabr model(maturity, forward, beta, alpha, nu, rho);

  EXPECT_DOUBLE_EQ(model.getMaturity(), maturity);
  EXPECT_DOUBLE_EQ(model.getForward(), forward);
  EXPECT_DOUBLE_EQ(model.getBeta(), beta);
  EXPECT_DOUBLE_EQ(model.getParameterAlpha(), alpha);
}

TEST(SabrTest, ImpliedVolCheck) {
  // Check that implied vol returns a reasonable number (non-negative, not NaN)
  double maturity = 1.0;
  double forward = 100.0;
  double beta = 1.0; // Lognormal
  double alpha = 0.2;
  double nu = 0.3;
  double rho = -0.2;

  Sabr model(maturity, forward, beta, alpha, nu, rho);

  // Use proper public API call
  // Note: Depends on actual implementation in sabr.cpp
  // Assuming impliedVol(strike) exists based on header view earlier
  double strike = 100.0; // ATM
  double vol = model.impliedVol(strike);

  EXPECT_GT(vol, 0.0);
  std::cout << "ATM Implied Vol: " << vol << std::endl;
}

TEST(SabrTest, NormalVolCheck) {
  double maturity = 1.0;
  double forward = 0.05;
  double beta = 0.5;
  double alpha = 0.2;
  double nu = 0.4;
  double rho = -0.3;

  Sabr model(maturity, forward, beta, alpha, nu, rho);

  double strike = 0.05; // ATM
  double normalVol = model.normalVol(strike);

  EXPECT_GT(normalVol, 0.0);
  EXPECT_LT(normalVol, 1.0); // Normal vol should be reasonable
}

TEST(SabrTest, PremiumBlackScholes) {
  double maturity = 1.0;
  double forward = 100.0;
  double beta = 1.0;
  double alpha = 0.2;
  double nu = 0.3;
  double rho = -0.2;

  Sabr model(maturity, forward, beta, alpha, nu, rho);

  double strike = 100.0;
  double callPremium = model.premiumBlackScholes(strike, OptionType::Call);
  double putPremium = model.premiumBlackScholes(strike, OptionType::Put);

  EXPECT_GT(callPremium, 0.0);
  EXPECT_GT(putPremium, 0.0);

  // Put-call parity (approximately: C - P = F - K for undiscounted)
  // Since F = K for ATM, C and P should be similar
  EXPECT_NEAR(callPremium, putPremium, callPremium * 0.1);
}

TEST(SabrTest, LocalVolCheck) {
  double maturity = 1.0;
  double forward = 100.0;
  double beta = 0.5;
  double alpha = 0.2;
  double nu = 0.3;
  double rho = -0.2;

  Sabr model(maturity, forward, beta, alpha, nu, rho);

  double spot = 100.0;
  double localVol = model.localVolIV(spot);

  EXPECT_GT(localVol, 0.0);
  EXPECT_LT(localVol, 2.0); // Local vol should be reasonable
}
