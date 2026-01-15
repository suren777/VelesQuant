from velesquant.models.localvol import LocalVolModel
from velesquant.models.sabr import SabrModel


def test_local_vol():
    s1 = SabrModel(maturity=1.0, forward=100.0, beta=0.85, alpha=0.2, nu=0.3, rho=-0.5)
    lv = LocalVolModel([s1], spot=100.0)
    price = lv.call_pde(0.5, 100.0, 50)
    assert price > 0


def test_local_vol_extended():
    s1 = SabrModel(maturity=1.0, forward=100.0, beta=0.85, alpha=0.2, nu=0.3, rho=-0.5)
    lv = LocalVolModel([s1])

    # Test spot property
    lv.spot = 105.0
    assert lv.spot == 105.0

    # Test put_pde
    put_price = lv.put_pde(0.5, 100.0, 50)
    assert put_price >= 0

    # Test density
    d = lv.density(1.0, 50)
    assert isinstance(d, list)
    assert len(d) > 0
    assert isinstance(d[0], float)
