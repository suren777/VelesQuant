import numpy as np

import velesquant.native as m


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
    val = q_cms.fairValue(0.03, m.OptionType.Call)
    assert val > 0.0

    # Test Forward
    # Quanto adjustment usually shifts the forward
    fwd_res = q_cms.getForward()
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

    sim_val = q_cms.simulation(0.0)
    assert isinstance(sim_val, float)
