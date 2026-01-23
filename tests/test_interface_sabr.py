import numpy as np

from velesquant import Sabr
from velesquant.models.sabr import SabrModel


def test_sabr_initialization():
    sabr = SabrModel(maturity=1.0, forward=0.03, alpha=0.2, beta=0.5, nu=0.4, rho=-0.3)
    assert sabr.alpha == 0.2
    assert sabr.rho == -0.3

    # Check C++ object valid
    vol = sabr.implied_vol(0.03)
    assert vol > 0.0


def test_sabr_calibration():
    # Synthetic data
    fwd = 0.03
    tenor = 1.0
    true_alpha = 0.2
    true_nu = 0.3
    true_rho = -0.2

    # Generate "market" quotes using a ground truth model
    ground_truth = Sabr(tenor, fwd, 0.5, true_alpha, true_nu, true_rho)

    strikes = [0.02, 0.025, 0.03, 0.035, 0.04]
    vols = [ground_truth.implied_vol(k) for k in strikes]

    # Initialize a perturbed model
    model = SabrModel(maturity=tenor, forward=fwd, alpha=0.5, beta=0.5, nu=0.5, rho=0.0)

    # Calibrate
    model.calibrate(strikes, vols, calibration_target="Volatility")

    # Check if params moved closer to truth
    assert np.isclose(model.alpha, true_alpha, atol=1e-2)
    assert np.isclose(model.nu, true_nu, atol=1e-2)
    assert np.isclose(model.rho, true_rho, atol=1e-2)
