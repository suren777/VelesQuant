import numpy as np

from velesquant.instruments.bonds import ZeroCouponBond
from velesquant.instruments.portfolio import Portfolio
from velesquant.instruments.rates import Swaption
from velesquant.market.curves import DiscountCurve
from velesquant.models.hullwhite import HullWhiteModel


def test_portfolio_pricing():
    # Setup
    curve = DiscountCurve([0, 10], [1.0, 0.61])
    model = HullWhiteModel(kappa=0.03, sigma=0.01)

    # Portfolio
    pf = Portfolio()

    # Add Swaptions
    s1 = Swaption(expiry=1.0, tenor=5.0, strike=0.03)
    s2 = Swaption(expiry=2.0, tenor=5.0, strike=0.03)
    pf.add(s1)
    pf.add(s2)

    # Add Bonds
    z1 = ZeroCouponBond(maturity=1.0)
    z2 = ZeroCouponBond(maturity=2.0)
    pf.add(z1)
    pf.add(z2)

    # Price
    # This should trigger the _price_portfolio batching logic
    total_npv = model.price(pf, curve)

    # Verify against individual pricing
    p_s1 = model.price(s1, curve)
    p_s2 = model.price(s2, curve)
    p_z1 = model.price(z1, curve)
    p_z2 = model.price(z2, curve)

    expected_total = p_s1 + p_s2 + p_z1 + p_z2

    assert np.isclose(total_npv, expected_total)
