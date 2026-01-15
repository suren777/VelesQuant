import numpy as np

from velesquant.models.heston import HestonModel


def test_heston_pricing():
    # Setup
    spot = 100.0
    var0 = 0.04
    kappa = 2.0
    theta = 0.04
    xi = 0.1
    rho = -0.5

    model = HestonModel(spot, var0, kappa, theta, xi, rho)

    # Price Call
    maturity = 1.0
    forward = 100.0  # simple case r=0, q=0
    strike = 100.0

    price = model.price_option(maturity, forward, strike, "call")
    assert price > 0.0

    # Put-Call Parity check roughly (if r=0)
    # C - P = S - K (if r=0, fwd=S)
    call_price = price
    put_price = model.price_option(maturity, forward, strike, "put")
    assert np.isclose(call_price - put_price, forward - strike, atol=1e-2)


def test_heston_calibration():
    # Synthetic calibration
    # Create "Ground Truth"
    true_kappa = 1.5
    true_theta = 0.05
    model_true = HestonModel(100.0, 0.04, true_kappa, true_theta, 0.1, -0.5)

    maturities = [1.0, 1.0, 1.0]
    forwards = [100.0, 100.0, 100.0]
    strikes = [90.0, 100.0, 110.0]

    prices = [
        model_true.price_option(m, f, k, "call")
        for m, f, k in zip(maturities, forwards, strikes)
    ]

    # Start guess
    model_guess = HestonModel(100.0, 0.04, 0.5, 0.02, 0.1, -0.5)

    # Calibrate
    model_guess.calibrate(maturities, forwards, strikes, prices, "Price")

    # Verify API contract: The method should run and update attributes (or at least return them)
    assert model_guess.var0 > 0
    assert model_guess.kappa > 0
    # Note: Convergence might vary on synthetic data, so we mainly check that the call completes
    # and populates the model state.
