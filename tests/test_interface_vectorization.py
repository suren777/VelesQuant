import numpy as np
import pytest

from velesquant.instruments.rates import Swaption
from velesquant.market.curves import DiscountCurve
from velesquant.models.hullwhite import HullWhiteModel


def test_swaption_vectorization_prototype():
    """
    Explore if we can pass numpy arrays to the existing Swaption class
    and if the Model can handle it (or be made to handle it).
    """
    # 1. Market
    curve = DiscountCurve([0, 10], [1.0, 0.6])

    # 2. Vectorized Instrument Inputs
    n_scenarios = 1000
    expiries = np.linspace(1.0, 5.0, n_scenarios)
    tenors = np.full(n_scenarios, 5.0)
    strikes = np.full(n_scenarios, 0.03)

    # Can we construct a "Vectorized Swaption"?
    # The current dataclass type hints say float, but at runtime it handles object/any.
    # However, the helper methods in the Model might fail if they expect scalar floats for C++ calls.

    inst = Swaption(expiry=expiries, tenor=tenors, strike=strikes)

    model = HullWhiteModel(kappa=0.03, sigma=0.01)

    # This call should now SUCCEED and return a numpy array
    try:
        prices = model.price(inst, curve)
        assert isinstance(prices, np.ndarray)
        assert len(prices) == n_scenarios
        assert np.all(prices > 0)
    except Exception as e:
        pytest.fail(f"Vectorization failed: {e}")


def test_bond_vectorization():
    from velesquant.instruments.bonds import ZeroCouponBond

    curve = DiscountCurve([0, 10], [1.0, 0.6])
    model = HullWhiteModel(kappa=0.03, sigma=0.01)

    maturities = np.array([1.0, 2.0, 3.0])
    # Create vectorized bond
    zcb_vec = ZeroCouponBond(maturity=maturities)

    prices = model.price(zcb_vec, curve)

    assert isinstance(prices, np.ndarray)
    assert len(prices) == 3
    assert np.all(prices > 0)
