import numpy as np

import velesquant.native as m


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

    sim = qcs.simulationQuantoedCMSs(0.0, 0.0)
    assert isinstance(sim, list)
    assert len(sim) == 2
