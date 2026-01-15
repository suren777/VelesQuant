import numpy as np

import velesquant.native as m


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

    cs = m.CmsSpread(
        exp1,
        ten1,
        fwd1,
        ann1,
        pay1,
        disc1,
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
        beta2,
        strikes2,
        quotes2,
        m.CalibrationTarget.Price,
        corr,
    )

    # Test Spread Option
    val = cs.spreadOption(0.01, 1.0, -1.0)
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

    cs = m.CmsSpread(
        exp1,
        ten1,
        fwd1,
        ann1,
        pay1,
        disc1,
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
        beta2,
        strikes2,
        quotes2,
        m.CalibrationTarget.Price,
        0.0,
    )

    sim_res = cs.simulationCMSs()
    assert isinstance(sim_res, list)
    assert len(sim_res) == 2
