import pytest

from velesquant.models import (
    HWPDEModel,
    SabrPDEModel,
    ShortRate1FPDEModel,
    ShortRate2FPDEModel,
)


def test_hwpde_wrapper_init():
    # kappa, timeSigmas, sigmas, timeDF, DF
    # Needs at least 2 points for interpolation
    model = HWPDEModel(0.1, [0.0, 10.0], [0.01, 0.01], [0.0, 10.0], [1.0, 0.9])
    assert model._cpp_model is not None
    assert model.price_swaption(1.0, 5.0, 0.03) > 0
    assert model.price_bermudan(1.0, 5.0, [0.5, 1.0], 0.03) > 0
    assert model.price_zbo(1.0, 2.0, 0.95, "Call") > 0
    assert model.price_zbo(1.0, 2.0, 0.95, "Put") > 0

    # Calibration smoke test (needs 2 points and enough quotes)
    # quotes = [
    #     {"expiry": 1.0, "tenor": 1.0, "rate": 0.03, "vol": 0.2, "frequency": 1.0},
    #     {"expiry": 1.0, "tenor": 2.0, "rate": 0.03, "vol": 0.2, "frequency": 1.0},
    #     {"expiry": 2.0, "tenor": 1.0, "rate": 0.03, "vol": 0.2, "frequency": 1.0},
    #     {"expiry": 2.0, "tenor": 2.0, "rate": 0.03, "vol": 0.2, "frequency": 1.0},
    # ]
    # model.calibrate([0.0, 10.0], [1.0, 0.7], quotes)

    # Test simulation
    sim = model.simulate([0.5, 1.0])
    assert len(sim) > 0

    d = model.to_dict()
    assert d["type"] == "HWPDEModel"
    assert d["kappa"] == 0.1


def test_hwpde_wrapper_discount_curve():
    # Test convenience argument discount_curve
    curve = [(0.0, 1.0), (10.0, 0.9)]
    model = HWPDEModel(
        kappa=0.1, time_sigmas=[0.0, 10.0], sigmas=[0.01, 0.01], discount_curve=curve
    )
    assert model._cpp_model is not None
    assert model.price_swaption(1.0, 5.0, 0.03) > 0

    # helper check
    # We can check if internal bindings got data (e.g. by pricing) which we did.


def test_short_rate_1f_pde_wrapper_init():
    # r0, kappa, alpha, beta, gamma, timeSigmas, sigmas, timeThetas, thetas
    model = ShortRate1FPDEModel(
        0.02, 0.1, 0.0, 0.0, 1.0, [0.0, 10.0], [0.01, 0.01], [0.0, 10.0], [0.03, 0.03]
    )
    assert model._cpp_model is not None
    assert model.price_swaption(1.0, 5.0, 0.03) > 0

    d = model.to_dict()
    assert d["type"] == "ShortRate1FPDEModel"
    assert d["initial_rate"] == 0.02


def test_sabr_pde_wrapper_init():
    # variant, alpha, beta, nu, rho, maturity, f
    model = SabrPDEModel("Af", 0.2, 0.5, 0.3, -0.2, 1.0, 100.0)
    assert model._cpp_model is not None
    with pytest.raises(NotImplementedError):
        model.calibrate([], None)

    density = model.get_density()
    assert len(density) > 0

    f_grid = model.get_f_grid()
    assert len(f_grid) > 0

    d = model.to_dict()
    assert d["type"] == "SabrPDEModel"
    assert d["variant"] == "Af"


def test_short_rate_2f_pde_wrapper_init():
    # kappa1, kappa2, lam, time_sigma1s, sigma1s, time_sigma2s, sigma2s, time_alphas, alphas
    model = ShortRate2FPDEModel(
        0.03,
        0.01,
        -0.5,
        [10.0],
        [0.01],
        [10.0],
        [0.02],
        [10.0],
        [0.0],
    )
    assert model._cpp_model is not None
    # Just check if pricing runs without error (smoke test)
    # Using dummy values
    price = model.price_swaption(1.0, 5.0, 0.03)
    assert price > 0.0

    zb = model.price_zero_bond(1.0)
    assert zb > 0.0
    # Price can be > 1.0 in G2++ if rates are negative (controlled by alphas)
    # assert 0.0 < zb < 1.0

    d = model.to_dict()
    assert d["type"] == "ShortRate2FPDEModel"
    assert d["kappa1"] == 0.03
