import pytest
import numpy as np

import velesquant.native as m
from velesquant.models import QuantoedCMSSpreadModel


def test_quantoed_cms_spread_init_and_sim():
    # Leg 1
    exp1, ten1, fwd1, ann1 = 5.0, 10.0, 0.03, 1.0
    pay1, disc1 = 0.0, 0.8
    corFX1, atmVolFX1 = 0.3, 0.15
    beta1 = 0.5
    strikes1 = np.array([[0.02]])
    quotes1 = np.array([[0.20]])

    # Leg 2
    exp2, ten2, fwd2, ann2 = 5.0, 2.0, 0.02, 1.0
    pay2, disc2 = 0.0, 0.8
    corFX2, atmVolFX2 = 0.2, 0.10
    beta2 = 0.5
    strikes2 = np.array([[0.02]])
    quotes2 = np.array([[0.25]])

    corr = 0.6

    qcs = m.QuantoedCmsSpread(
        exp1,
        ten1,
        fwd1,
        ann1,
        pay1,
        disc1,
        corFX1,
        atmVolFX1,
        beta1,
        strikes1,
        quotes1,
        m.CalibrationTarget.Price,
        exp2,
        ten2,
        fwd2,
        ann2,
        pay2,
        disc2,
        corFX2,
        atmVolFX2,
        beta2,
        strikes2,
        quotes2,
        m.CalibrationTarget.Price,
        corr,
    )

    sim = qcs.simulate(0.0, 0.0)
    assert isinstance(sim, list)
    assert len(sim) == 2


def test_quantoed_cms_spread_wrapper_init():
    strikes = np.array([0.01, 0.02, 0.03]).reshape(-1, 1)
    quotes = np.array([0.2, 0.2, 0.2]).reshape(-1, 1)
    model = QuantoedCMSSpreadModel(
        1.0,
        10.0,
        0.03,
        8.5,
        1.0,
        0.97,
        0.3,
        0.15,
        0.85,
        strikes,
        quotes,
        "Price",
        1.0,
        10.0,
        0.03,
        8.5,
        1.0,
        0.97,
        0.3,
        0.15,
        0.85,
        strikes,
        quotes,
        "Price",
        0.5,
    )
    assert model._cpp_model is not None

    # Test to_dict
    d = model.to_dict()
    assert d["type"] == "QuantoedCMSSpreadModel"
    assert d["corr"] == 0.5

    # Test simulate (2 args)
    sim = model.simulate(0.5, 0.5)
    assert len(sim) > 0

    with pytest.raises(NotImplementedError):
        model.price(None, None)

    with pytest.raises(NotImplementedError):
        model.calibrate([], None)
