import pytest
from velesquant import native
from velesquant.models import SchobelZhuModel


def test_schobzhu_binding():
    # SchobelZhu(spot, var0, kappa, theta, xi, rho)
    sz = native.SchobelZhu(100.0, 0.04, 1.0, 0.04, 0.3, -0.5)
    price = sz.SchobelPrice(1.0, 100.0, 100.0)
    assert isinstance(price, float)
    assert price > 0

    sz.var0 = 0.05
    assert sz.var0 == 0.05


def test_schobzhu_simulation():
    """
    Scenario: Monte Carlo simulation of Schobel-Zhu price paths.
    """
    sz = native.SchobelZhu(100.0, 0.04, 1.0, 0.04, 0.3, -0.5)
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
    sz = native.SchobelZhu(100.0, 0.04, 1.0, 0.04, 0.3, -0.5)

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

    model = SchobelZhuModel(100.0, 0.04, 2.0, 0.2, 0.3, -0.5)
    assert model.spot == 100.0
    assert model.kappa == 2.0
    assert model.price_european(1.0, 100.0, 100.0) > 0

    # Test calibration in wrapper
    model.calibrate_to_quotes(
        maturities=[1.0, 1.0, 1.0],
        forwards=[100.0, 100.0, 100.0],
        strikes=[90.0, 100.0, 110.0],
        quotes=[12.0, 5.0, 2.0],
        target="Price",
    )
    assert model.var0 > 0

    # Test Properties
    model.spot = 105.0
    assert model.spot == 105.0
    model.var0 = 0.05
    assert model.var0 == 0.05
    model.kappa = 2.5
    assert model.kappa == 2.5
    model.theta = 0.25
    assert model.theta == 0.25
    model.xi = 0.35
    assert model.xi == 0.35
    model.rho = -0.45
    assert model.rho == -0.45

    # Test to_dict
    d = model.to_dict()
    assert d["type"] == "SchobelZhuModel"
    assert d["spot"] == 105.0
    assert d["rho"] == -0.45

    # Test simulate
    sim = model.simulate([1.0], [105.0])
    assert len(sim) > 0

    # Test Not Implemented
    with pytest.raises(NotImplementedError):
        model.calibrate([], None)
    with pytest.raises(NotImplementedError):
        model.price(None, None)
