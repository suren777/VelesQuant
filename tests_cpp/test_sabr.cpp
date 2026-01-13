#include <gtest/gtest.h>
#include <velesquant/local_vol/sabr.h>
#include <iostream>

using namespace velesquant;

TEST(SabrTest, Construction) {
    // Basic construction check
    // sabr(double maturity, double forward, double beta, double alpha, double nu, double rho)
    double maturity = 1.0;
    double forward = 0.05;
    double beta = 0.5;
    double alpha = 0.2;
    double nu = 0.4;
    double rho = -0.3;

    sabr model(maturity, forward, beta, alpha, nu, rho);

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

    sabr model(maturity, forward, beta, alpha, nu, rho);
    
    // Use proper public API call
    // Note: Depends on actual implementation in sabr.cpp
    // Assuming impliedVol(strike) exists based on header view earlier
    double strike = 100.0; // ATM
    double vol = model.impliedVol(strike);

    EXPECT_GT(vol, 0.0);
    std::cout << "ATM Implied Vol: " << vol << std::endl;
}
