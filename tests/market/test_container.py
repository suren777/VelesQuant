import pytest

from velesquant.market.container import Market
from velesquant.market.curves import DiscountCurve


def test_market_container_basic():
    market = Market()

    # Create some dummy data
    curve = DiscountCurve(times=[0.0, 1.0], dfs=[1.0, 0.95])

    # Add to market
    market.add("USD_Curve", curve)

    # Retrieve
    retrieved = market.get("USD_Curve", DiscountCurve)
    assert retrieved is curve
    assert retrieved.discount(0.5) < 1.0


def test_market_container_missing():
    market = Market()
    with pytest.raises(KeyError):
        market.get("NonExistent", DiscountCurve)


def test_market_container_collisions():
    """Test that same name but different types are stored separately."""
    market = Market()

    data1 = "StringData"
    data2 = 123

    market.add("SharedName", data1)
    market.add("SharedName", data2)

    assert market.get("SharedName", str) == "StringData"
    assert market.get("SharedName", int) == 123
