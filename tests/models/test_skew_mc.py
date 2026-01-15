"""Tests for SkewMC (Skew Monte Carlo) model."""

from velesquant import native
from velesquant.models import SkewMCModel, SabrModel


def test_skewmc_binding():
    """Test that SkewMC can be instantiated with SABR models."""
    # Create a list of SABR models for different maturities
    sabr1 = native.Sabr(0.5, 0.05, 0.5, 0.2, 0.3, -0.2)
    sabr2 = native.Sabr(1.0, 0.05, 0.5, 0.2, 0.3, -0.2)
    sabr3 = native.Sabr(2.0, 0.05, 0.5, 0.2, 0.3, -0.2)

    mc = native.SkewMC([sabr1, sabr2, sabr3])
    assert mc is not None


def test_skewmc_simulation():
    """Test Monte Carlo simulation execution."""
    sabr1 = native.Sabr(0.5, 100.0, 0.5, 0.2, 0.3, -0.2)
    sabr2 = native.Sabr(1.0, 100.0, 0.5, 0.2, 0.3, -0.2)

    mc = native.SkewMC([sabr1, sabr2])

    # Simulate at times 0.25, 0.5, 0.75, 1.0
    times = [0.25, 0.5, 0.75, 1.0]
    spot = 100.0
    kappa = 0.1

    path = mc.simulate(times, spot, kappa)
    assert len(path) == len(times)
    # Spot values should be positive
    assert all(s > 0 for s in path)


def test_skew_mc_wrapper_init():
    sabr = SabrModel(1.0, 0.05, 0.5, 0.2, 0.3, -0.2)
    model = SkewMCModel([sabr])
    assert model._cpp_model is not None
    path = model.simulate([0.5, 1.0], 100.0, 0.1)
    assert len(path) == 2

    # Test to_dict
    d = model.to_dict()
    assert d["type"] == "SkewMCModel"
    assert len(d["sabr_models"]) == 1
