import velesquant._core as m
import math


def test_sabr():
    print("Testing SABR model binding...")
    # sabr(double maturity, double forward, double beta, double alpha, double nu, double rho)
    # Using defaults: beta=0.85, alpha=0.5, nu=0.25, rho=-0.75
    maturity = 1.0
    forward = 0.05
    beta = 0.85
    alpha = 0.2
    nu = 0.3
    rho = -0.2

    sabr = m.Sabr(maturity, forward, beta, alpha, nu, rho)

    print(f"SABR Created: Maturity={sabr.getMaturity()}, Forward={sabr.getForward()}")
    print(
        f"Alpha={sabr.getParameterAlpha()}, Nu={sabr.getParameterNu()}, Rho={sabr.getParameterRho()}"
    )

    strike = 0.05
    iv = sabr.impliedVol(strike)
    print(f"Implied Vol at strike {strike}: {iv}")

    # Check modification
    sabr.alpha = 0.25
    print(f"New Alpha: {sabr.alpha}")
    iv2 = sabr.impliedVol(strike)
    print(f"Implied Vol at strike {strike} with new alpha: {iv2}")

    assert iv > 0
    assert iv != iv2
    print("Verification Successful!")


if __name__ == "__main__":
    test_sabr()
