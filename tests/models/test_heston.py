import numpy as np

from velesquant.models.heston import HestonModel


def test_heston():
    # spot, var0, kappa, theta, xi, rho, seed
    h = HestonModel(
        spot=100.0, var0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7, seed=42
    )
    # price_option(maturity, forward, strike, option_type)
    # Original test: hestonPrice(1.0, 100.0, 100.0, "call")
    price = h.price_option(1.0, 100.0, 100.0, "call")
    assert price > 0
    assert h.var0 == 0.04
    assert h.kappa == 2.0


def test_heston_extended():
    h = HestonModel(
        spot=100.0, var0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7, seed=42
    )

    # Test properties
    assert h.theta == 0.04
    assert h.xi == 0.3
    assert h.rho == -0.7

    # Test simulation
    times = [0.0, 0.5, 1.0]
    # forwards must match times size or logic
    forwards = [100.0, 100.0, 100.0]

    paths = h.simulate(times, forwards)
    # Returns vector<double> containing spot path
    assert isinstance(paths, list)
    assert len(paths) == 3
    assert paths[0] > 0.0  # Start at spot


def test_heston_calibration_example():
    """
    Example of using the new Object-Oriented Heston API for calibration.
    """
    # Create "Ground Truth" to generate consistent data
    true_kappa = 1.5
    true_theta = 0.05
    model_true = HestonModel(100.0, 0.04, true_kappa, true_theta, 0.1, -0.5)

    # 1. Define Market Surface Data
    maturities = np.array([1.0] * 9).reshape(-1, 1)
    forwards = np.array([100.0] * 9).reshape(-1, 1)
    strikes_list = [80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0]
    strikes = np.array(strikes_list).reshape(-1, 1)

    prices_list = [model_true.price_option(1.0, 100.0, k, "call") for k in strikes_list]
    # Quotes are prices here because target="Price"
    quotes = np.array(prices_list).reshape(-1, 1)

    # 2. Instantiate the Heston Model (Initial Guess)
    heston = HestonModel(
        spot=100.0, var0=0.04, kappa=1.2, theta=0.04, xi=0.1, rho=-0.5, seed=42
    )

    # 3. Running Calibration
    result_model = heston.calibrate(
        maturities=maturities,
        forwards=forwards,
        strikes=strikes,
        quotes=quotes,
        calibration_target="Price",
    )

    # 4. Verify Results
    # The wrapper updates attributes

    # Check that the object was updated
    assert heston.var0 == result_model.var0
    assert abs(heston.rho) <= 1.0
