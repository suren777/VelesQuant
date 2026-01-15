import pytest

from velesquant.instruments.rates import Swaption
from velesquant.market.container import Market
from velesquant.market.curves import DiscountCurve
from velesquant.models.hullwhite import HullWhiteModel


def test_hullwhite_with_market_container():
    # Setup
    curve = DiscountCurve(times=[0.0, 10.0], dfs=[1.0, 0.6])
    market = Market()
    market.add("USD", curve)

    model = HullWhiteModel(kappa=0.03, sigma=0.01)
    inst = Swaption(expiry=1.0, tenor=5.0, strike=0.03)

    # Price using Market container
    price_via_market = model.price(inst, market_data=market, curve_name="USD")

    # Price using direct Curve
    price_direct = model.price(inst, market_data=curve)

    assert price_via_market == price_direct
    assert price_via_market > 0.0


def test_hullwhite_with_market_container_missing_key():
    market = Market()
    model = HullWhiteModel(kappa=0.03, sigma=0.01)
    inst = Swaption(expiry=1.0, tenor=5.0, strike=0.03)

    with pytest.raises(KeyError):
        model.price(inst, market_data=market, curve_name="EUR")


def test_hullwhite_invalid_market_type():
    model = HullWhiteModel(kappa=0.03, sigma=0.01)
    inst = Swaption(expiry=1.0, tenor=5.0, strike=0.03)

    with pytest.raises(TypeError):
        model.price(inst, market_data="InvalidString")
