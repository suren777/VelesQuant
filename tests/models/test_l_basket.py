import pytest

from velesquant import LogNormalBasket
from velesquant.models import LogNormalBasketModel


def test_log_basket_init_and_sim():
    # 2 Assets
    spots = [100.0, 100.0]
    strikes = [100.0, 100.0]
    maturities = [1.0, 1.0]

    forwards = [[100.0], [100.0]]
    iv = [[0.2], [0.2]]
    correlation = [[1.0, 0.5], [0.5, 1.0]]

    lb = LogNormalBasket(spots, strikes, maturities, forwards, iv, correlation)

    assert lb.get_n_assets() == 2

    schedule = [1.0]
    paths = lb.simulate(schedule)
    assert isinstance(paths, list)
    assert len(paths) > 0


def test_basket_wrapper_init():
    spots = [100.0, 110.0]
    strikes = [100.0, 110.0]
    maturities = [1.0, 2.0]
    fwds = [[100.0, 100.0], [110.0, 110.0]]
    ivs = [[0.2, 0.2], [0.22, 0.22]]
    corr = [[1.0, 0.5], [0.5, 1.0]]
    model = LogNormalBasketModel(spots, strikes, maturities, fwds, ivs, corr)
    assert model.get_num_assets() == 2

    # Test simulation variants
    sim1 = model.simulate([1.0])
    assert len(sim1) > 0
    sim2 = model.simulate([1.0], with_rate_curves=True)
    assert len(sim2) > 0

    # Calibration stub test
    with pytest.raises(NotImplementedError):
        model.calibrate([], None)

    with pytest.raises(NotImplementedError):
        model.price(None, None)

    # Test to_dict
    d = model.to_dict()
    assert d["type"] == "LogNormalBasketModel"
    assert d["spots"] == spots
