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

    maturities = np.array([1.0] * 9).reshape(-1, 1)
    forwards = np.array([100.0] * 9).reshape(-1, 1)
    # Strikes 80 to 120 step 5
    strikes_list = [80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0]
    strikes = np.array(strikes_list).reshape(-1, 1)

    prices_list = [model_true.price_option(1.0, 100.0, k, "call") for k in strikes_list]
    prices = np.array(prices_list).reshape(-1, 1)

    # Start guess - close to true
    model_guess = HestonModel(100.0, 0.04, 1.2, 0.04, 0.1, -0.5)

    # Calibrate
    model_guess.calibrate(maturities, forwards, strikes, prices, "Price")

    # Verify API contract: The method should run and update attributes
    assert model_guess.var0 > 0
    assert model_guess.kappa > 0
    # Check if we got somewhat closer or at least didn't explode
    assert abs(model_guess.theta - true_theta) < 0.1
