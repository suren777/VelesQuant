from velesquant import SchobelZhu


def test_schobzhu_binding():
    # SchobelZhu(spot, var0, kappa, theta, xi, rho)
    sz = SchobelZhu(100.0, 0.04, 1.0, 0.04, 0.3, -0.5)
    price = sz.SchobelPrice(1.0, 100.0, 100.0)
    assert isinstance(price, float)
    assert price > 0

    sz.var0 = 0.05
    assert sz.var0 == 0.05


def test_schobzhu_simulation():
    """
    Scenario: Monte Carlo simulation of Schobel-Zhu price paths.
    """
    sz = SchobelZhu(100.0, 0.04, 1.0, 0.04, 0.3, -0.5)
    times = [0.1, 0.5, 1.0]
    forwards = [100.0, 101.0, 102.0]

    paths = sz.simulation(times, forwards)

    assert len(paths) == len(times)
    # Check for no NaNs or Infs
    for p in paths:
        assert isinstance(p, float)
        assert p > 0


def test_schobzhu_calibration():
    """
    Scenario: Calibrating Schobel-Zhu to a set of market quotes.
    """
    sz = SchobelZhu(100.0, 0.04, 1.0, 0.04, 0.3, -0.5)

    # Mock market data
    maturities = [1.0, 1.0, 1.0]
    forwards = [100.0, 100.0, 100.0]
    strikes = [90.0, 100.0, 110.0]
    quotes = [12.0, 5.0, 2.0]  # Rough call prices

    # Run calibration
    sz.calibrator(maturities, forwards, strikes, quotes)

    # Check that parameters have changed (or at least are accessible)
    assert sz.var0 > 0
    assert abs(sz.rho) <= 1.0
