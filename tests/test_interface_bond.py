import numpy as np

from velesquant.instruments.bonds import CouponBond, ZeroCouponBond
from velesquant.market.curves import DiscountCurve
from velesquant.models.hullwhite import HullWhiteModel


def test_bond_pricing():
    # Setup Market
    r = 0.05
    # Use denser grid to avoid Cubic Spline artifacts on DFs
    times = list(range(11))  # 0 to 10
    dfs = [np.exp(-r * t) for t in times]
    curve = DiscountCurve(times, dfs)

    # Setup Model
    model = HullWhiteModel(kappa=0.01, sigma=0.01)

    # 1. Zero Coupon Bond
    zcb = ZeroCouponBond(maturity=1.0)
    p_zcb = model.price(zcb, curve)

    # Analytic check P(0,T) = exp(-rT)
    expected_zcb = np.exp(-r * 1.0)  # ~ 0.9512
    assert np.isclose(p_zcb, expected_zcb, atol=1e-4)

    # 2. Coupon Bond
    cb = CouponBond(maturity=2.0, pay_frequency=1.0, coupon_rate=0.05)
    p_cb = model.price(cb, curve)

    # Cashflows:
    # T=1: 0.05
    # T=2: 1.05
    # Price = 0.05 * P(0,1) + 1.05 * P(0,2)
    # We can check against manual calc using model.price(zcb) logic logic

    zcb1 = ZeroCouponBond(maturity=1.0)
    zcb2 = ZeroCouponBond(maturity=2.0)

    p1 = model.price(zcb1, curve)
    p2 = model.price(zcb2, curve)

    # Consistency check
    assert np.isclose(p_cb, 0.05 * p1 + 1.05 * p2)
