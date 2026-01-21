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


def test_local_vol_basic():
    # Setup inputs mocking what facade expected
    # Use milder params (LogNormal beta=1.0)
    maturities = [1.0, 5.0]
    forwards = [100.0, 100.0]
    betas = [1.0, 1.0]
    alphas = [0.2, 0.2]
    nus = [0.1, 0.1]
    rhos = [-0.3, -0.3]
    spot = 100.0

    lv = LocalVolModel(
        maturities=maturities,
        forwards=forwards,
        betas=betas,
        alphas=alphas,
        nus=nus,
        rhos=rhos,
        spot=spot,
    )

    # Test Density
    # Increase n_points for stability
    density = lv.get_density(maturity=1.0, n_points=100)
    assert len(density) > 0
    # Density should integrate somewhat to 1 (ignoring grid range cutoffs)?
    # Just check values are non-negative, allowing for small numerical noise
    assert all(d >= -1e-5 for d in density)

    # Test Pricing
    call_price = lv.price_european(maturity=1.0, strike=100.0, option_type="Call")
    assert call_price > 0

    put_price = lv.price_european(maturity=1.0, strike=100.0, option_type="Put")
    assert put_price > 0

    # Test Barrier
    dnt_price = lv.price_barrier(maturity=1.0, upper_barrier=120.0, lower_barrier=80.0)
    assert dnt_price >= 0
    # Naive check: DNT price should be less than sum of vanilla calls/puts (roughly) or positive
    # DNT is a touch option? No, Double No Touch pays if never touches.
    # So it pays 1 if spot stays inside.
    assert dnt_price <= 1.0

    # Test Export Surface
    times = [0.5, 1.0, 1.5]
    surface = lv.export_local_vol_surface(times)
    assert len(surface) == 3
    for slice_vol in surface:
        assert len(slice_vol) > 0


def test_local_vol_simulation():
    maturities = [1.0]
    forwards = [100.0]
    betas = [1.0]
    alphas = [0.2]
    nus = [0.1]
    rhos = [-0.3]
    spot = 100.0

    lv = LocalVolModel(
        maturities=maturities,
        forwards=forwards,
        betas=betas,
        alphas=alphas,
        nus=nus,
        rhos=rhos,
        spot=spot,
    )
    times = [0.0, 0.5, 1.0]

    path = lv.simulate(times)
    assert len(path) == 3
    assert path[0] == spot
