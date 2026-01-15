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
    # 1. Define Market Surface Data
    # Heston requires a volatility surface (maturity x strike)
    maturities = np.array([[1.0]])
    forwards = np.array([[100.0]])
    strikes = np.array([[100.0]])
    quotes = np.array([[0.20]])

    # 2. Instantiate the Heston Model
    heston = HestonModel(
        spot=100.0, var0=0.04, kappa=1.0, theta=0.04, xi=0.1, rho=-0.5, seed=42
    )

    # 3. Running Calibration

    # Wrapper returns updated self OR params dict?
    # Wrapper implementation returns self.

    result_model = heston.calibrate(
        maturities=maturities,
        forwards=forwards,
        strikes=strikes,
        quotes=quotes,
        calibration_target="Price",  # Or "Volatility"
    )

    # 4. Verify Results
    # The wrapper updates attributes

    # Check that the object was updated
    assert heston.var0 == result_model.var0
    assert abs(heston.rho) <= 1.0
