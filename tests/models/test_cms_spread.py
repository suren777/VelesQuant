import numpy as np
import pytest

from velesquant import CalibrationTarget, CmsSpread
from velesquant.models import CMSSpreadModel


def test_cms_spread_init_and_pricing():
    # Construct Inputs for Leg 1
    exp1, ten1, fwd1, ann1 = 5.0, 10.0, 0.03, 1.0
    pay1, disc1, beta1 = 0.0, 0.8, 0.5
    strikes1 = np.array([[0.02], [0.03], [0.04]])
    quotes1 = np.array([[0.20], [0.20], [0.20]])  # Flat vol 20%

    # Construct Inputs for Leg 2
    exp2, ten2, fwd2, ann2 = 5.0, 2.0, 0.02, 1.0
    pay2, disc2, beta2 = 0.0, 0.8, 0.5
    strikes2 = np.array([[0.01], [0.02], [0.03]])
    quotes2 = np.array([[0.25], [0.25], [0.25]])  # Flat vol 25%

    corr = 0.6

    cs = CmsSpread(
        exp1,
        ten1,
        fwd1,
        ann1,
        pay1,
        disc1,
        beta1,
        strikes1,
        quotes1,
        CalibrationTarget.Price,
        exp2,
        ten2,
        fwd2,
        ann2,
        pay2,
        disc2,
        beta2,
        strikes2,
        quotes2,
        CalibrationTarget.Price,
        corr,
    )

    # Test Spread Option
    val = cs.spread_option(0.01, 1.0, -1.0)
    assert val > 0.0


def test_cms_spread_simulation():
    # Simplified inputs for simulation test
    exp1, ten1, fwd1, ann1 = 1.0, 1.0, 0.03, 1.0
    pay1, disc1, beta1 = 0.0, 0.9, 0.5
    strikes1 = np.array([[0.03]])
    quotes1 = np.array([[0.20]])

    exp2, ten2, fwd2, ann2 = 1.0, 1.0, 0.03, 1.0
    pay2, disc2, beta2 = 0.0, 0.9, 0.5
    strikes2 = np.array([[0.03]])
    quotes2 = np.array([[0.20]])

    cs = CmsSpread(
        exp1,
        ten1,
        fwd1,
        ann1,
        pay1,
        disc1,
        beta1,
        strikes1,
        quotes1,
        CalibrationTarget.Price,
        exp2,
        ten2,
        fwd2,
        ann2,
        pay2,
        disc2,
        beta2,
        strikes2,
        quotes2,
        CalibrationTarget.Price,
        0.0,
    )

    sim_res = cs.simulate()
    assert isinstance(sim_res, list)
    assert len(sim_res) == 2


def test_cms_spread_wrapper_init():
    strikes = np.array([0.01, 0.02, 0.03]).reshape(-1, 1)
    quotes = np.array([0.2, 0.2, 0.2]).reshape(-1, 1)
    model = CMSSpreadModel(
        1.0,
        10.0,
        0.03,
        8.5,
        1.0,
        0.97,
        0.85,
        strikes,
        quotes,
        "Price",
        1.0,
        10.0,
        0.02,
        8.5,
        1.0,
        0.97,
        0.85,
        strikes,
        quotes,
        "Price",
        0.5,
    )
    assert model._cpp_model is not None
    assert model.spread_option(0.01) > 0

    # Test to_dict
    d = model.to_dict()
    assert d["type"] == "CMSSpreadModel"
    assert d["corr"] == 0.5

    # Test simulate (no args in wrapper)
    sim = model.simulate()
    assert len(sim) == 2

    # Test price stub
    with pytest.raises(NotImplementedError):
        model.price(None, None)

    # Test Vol target
    model_vol = CMSSpreadModel(
        1.0,
        10.0,
        0.03,
        1.0,
        0.04,
        1.0,
        0.97,
        strikes,
        quotes,
        "Vol",
        1.0,
        10.0,
        0.02,
        1.0,
        0.04,
        1.0,
        0.97,
        strikes,
        quotes,
        "Vol",
        0.5,
    )
    assert model_vol._cpp_model is not None

    with pytest.raises(NotImplementedError):
        model.calibrate([], None)
