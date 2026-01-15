import json

from velesquant.models.hullwhite import HullWhiteModel


def test_hullwhite_serialization():
    model = HullWhiteModel(kappa=0.03, sigma=0.01)

    # Serialize
    data_dict = model.to_dict()
    assert data_dict["type"] == "HullWhiteModel"
    assert data_dict["kappa"] == 0.03
    assert data_dict["sigma"] == 0.01

    # JSON roundtrip
    json_str = json.dumps(data_dict)
    loaded = json.loads(json_str)

    # Deserialize
    restored = HullWhiteModel.from_dict(loaded)
    assert restored.kappa == 0.03
    assert restored.sigma == 0.01
