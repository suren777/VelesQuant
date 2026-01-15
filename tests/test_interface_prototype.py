import numpy as np

from velesquant.instruments.rates import Swaption
from velesquant.market.curves import DiscountCurve
from velesquant.models.hullwhite import HullWhiteModel


def test_interface_swaption_pricing():
    # 1. Setup Market Data
    # Flat curve at 5%
    times = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0]
    r = 0.05
    dfs = [np.exp(-r * t) for t in times]
    curve = DiscountCurve(times, dfs)

    # 2. Define Instrument
    # 5Y x 10Y Swaption, ATMish (roughly)
    swaption = Swaption(expiry=5.0, tenor=10.0, strike=0.05)

    # 3. Define Model
    model = HullWhiteModel(kappa=0.03, sigma=0.01)

    # 4. Price
    price = model.price(swaption, curve)

    # 5. Assertions
    assert price > 0.0

    # Check sensitivity (higher vol -> higher price)
    model_high_vol = HullWhiteModel(kappa=0.03, sigma=0.02)
    price_high_vol = model_high_vol.price(swaption, curve)
    assert price_high_vol > price
