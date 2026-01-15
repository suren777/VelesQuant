import numpy as np

from velesquant.market.curves import DiscountCurve


def test_discount_curve_interpolation():
    times = [0.0, 1.0, 2.0]
    dfs = [1.0, 0.9, 0.81]  # Roughly 10% rate

    curve = DiscountCurve(times, dfs)

    # Check exact points
    assert np.isclose(curve.discount(0.0), 1.0)
    assert np.isclose(curve.discount(1.0), 0.9)
    assert np.isclose(curve.discount(2.0), 0.81)

    # Check interpolation (t=1.5 should be geometric mean of 0.9 and 0.81 => sqrt(0.9*0.81) = 0.9*0.9 = 0.81?? No sqrt(0.729) approx 0.8538)
    # log(0.9) approx -0.1053, log(0.81) approx -0.2107
    # mid log = -0.158
    # exp(-0.158) = 0.8538..
    # 0.9 * 0.9 = 0.81? Wait
    # 0.9 = e^-r*1
    # 0.81 = e^-r*2 = (e^-r)^2 = 0.9^2. Correct.
    # discount(1.5) should be e^-r*1.5 = (0.9)^1.5 = 0.8538149

    expected_1_5 = 0.9**1.5
    assert np.isclose(curve.discount(1.5), expected_1_5)

    # Check simple extrapolation (just bounds check)
    assert np.isclose(curve.discount(-0.1), 1.0)  # Left bound
    assert np.isclose(curve.discount(2.1), 0.81)  # Right bound
