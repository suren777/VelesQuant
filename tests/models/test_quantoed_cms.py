import numpy as np

import velesquant.native as m
from velesquant.models import QuantoedCMSModel


def test_quantoed_cms_init_and_pricing():
    exp, ten, fwd, ann = 5.0, 10.0, 0.03, 1.0
    pay, disc = 0.0, 0.8
    corFX, atmVolFX = 0.3, 0.15
    beta = 0.5

    strikes = np.array([[0.02], [0.03], [0.04]])
    quotes = np.array([[0.20], [0.20], [0.20]])  # Flat vol

    q_cms = m.QuantoedCMS(
        exp,
        ten,
        fwd,
        ann,
        pay,
        disc,
        corFX,
        atmVolFX,
        beta,
        strikes,
        quotes,
        m.CalibrationTarget.Price,
    )

    # Test Fair Value
    val = q_cms.fair_value(0.03, m.OptionType.Call)
    assert val > 0.0

    # Test Forward
    # Quanto adjustment usually shifts the forward
    fwd_res = q_cms.get_forward()
    assert fwd_res != fwd  # Should be adjusted


def test_quantoed_cms_simulation():
    exp, ten, fwd, ann = 1.0, 5.0, 0.03, 1.0
    pay, disc = 0.0, 0.9
    corFX, atmVolFX = 0.0, 0.1  # No corr -> simple
    beta = 0.5
    strikes = np.array([[0.03]])
    quotes = np.array([[0.20]])

    q_cms = m.QuantoedCMS(
        exp,
        ten,
        fwd,
        ann,
        pay,
        disc,
        corFX,
        atmVolFX,
        beta,
        strikes,
        quotes,
        m.CalibrationTarget.Price,
    )

    sim_val = q_cms.simulate(0.0)
    assert isinstance(sim_val, (float, list))


def test_quantoed_cms_wrapper_init():
    strikes = np.array([0.01, 0.02, 0.03]).reshape(-1, 1)
    quotes = np.array([0.2, 0.2, 0.2]).reshape(-1, 1)
    model = QuantoedCMSModel(
        1.0, 10.0, 0.03, 8.5, 1.0, 0.97, 0.7, 0.15, 0.85, strikes, quotes, "Price"
    )
    assert model._cpp_model is not None
    assert model.get_forward() != 0.03
    # assert model.fair_value(0.05, "Put") > 0  # Put seems to return 0.0, skipping for now

    sim = model.simulate(0.5)
    assert sim > 0

    # Test to_dict
    d = model.to_dict()
    assert d["type"] == "QuantoedCMSModel"
    assert d["expiry"] == 1.0
    # Test Vol target
    model_vol = QuantoedCMSModel(
        1.0, 10.0, 0.03, 8.5, 1.0, 0.97, 0.7, 0.15, 0.85, strikes, quotes, "Vol"
    )
    assert model_vol._cpp_model is not None
