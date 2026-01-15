import json


from velesquant.market.container import Market
from velesquant.market.curves import DiscountCurve


def test_market_serialization():
    market = Market()
    curve = DiscountCurve(times=[0.0, 1.0], dfs=[1.0, 0.95])
    market.add("USD", curve)

    # Serialize
    data_dict = market.to_dict()
    assert "USD:DiscountCurve" in data_dict
    assert data_dict["USD:DiscountCurve"]["type"] == "DiscountCurve"

    # Test JSON compatibility
    json_str = json.dumps(data_dict)
    loaded_dict = json.loads(json_str)

    # Deserialize
    new_market = Market.from_dict(loaded_dict)
    retrieved = new_market.get("USD", DiscountCurve)

    # Verify content
    assert retrieved.times == [0.0, 1.0]
    assert retrieved.dfs == [1.0, 0.95]
    assert abs(retrieved.discount(0.5) - curve.discount(0.5)) < 1e-9
