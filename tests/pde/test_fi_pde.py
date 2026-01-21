"""Tests for Fixed Income PDE Solvers (HWPDE, ShortRate1FPDE, ShortRate2FPDE)."""

from velesquant import HWPDE, ShortRate1FPDE, ShortRate2FPDE
from velesquant.native import HullWhiteModel, ShortRate1FModel, ShortRate2FModel


def test_hwpde_binding():
    """Test that HWPDE can be instantiated."""
    # Parameters for Hull-White PDE
    kappa = 0.05
    timeSigmas = [0.5, 1.0, 2.0]
    sigmas = [0.01, 0.012, 0.015]
    time_dfs = [0.5, 1.0, 2.0]
    dfs = [0.98, 0.95, 0.9]

    model = HullWhiteModel(kappa, timeSigmas, sigmas, time_dfs, dfs)
    hwpde = HWPDE(model)
    assert hwpde is not None


def test_hwpde_swaption():
    """Test swaption pricing with HWPDE."""
    R0 = 0.03
    kappa = 0.05
    timeSigmas = [1.0, 2.0, 5.0]
    sigmas = [0.01, 0.012, 0.015]
    time_dfs = [0.0, 0.5, 1.0, 5.0]
    dfs = [1.0, 0.99, 0.95, 0.80]

    model = HullWhiteModel(kappa, timeSigmas, sigmas, time_dfs, dfs)
    hwpde = HWPDE(model)

    # Price a 1Y into 5Y swaption
    price = hwpde.price_swaption(1.0, 5.0, 0.03, 0.5)
    assert isinstance(price, float)
    # Swaption price should be positive
    assert price >= 0.0


def test_short_rate_1f_pde_binding():
    """Test ShortRate1FPDE instantiation."""
    R0 = 0.03
    kappa = 0.05
    alpha = 1.0
    beta = 0.0001
    gamma = 1.0
    timeSigmas = [1.0, 2.0]
    sigmas = [0.01, 0.012]
    # timeThetas = [1.0, 2.0]
    # thetas = [0.03, 0.035]

    model = ShortRate1FModel(R0, kappa, alpha, beta, gamma, timeSigmas, sigmas)
    sr1f = ShortRate1FPDE(model)
    assert sr1f is not None


def test_short_rate_2f_pde_binding():
    """Test ShortRate2FPDE instantiation."""
    kappa1 = 0.03
    kappa2 = 0.05
    lam = 0.2
    timeSigma1s = [1.0, 2.0]
    sigma1s = [0.01, 0.012]
    timeSigma2s = [1.0, 2.0]
    sigma2s = [0.005, 0.006]
    timeAlphas = [1.0, 2.0]
    alphas = [0.03, 0.035]

    # No R0
    model = ShortRate2FModel(
        kappa1,
        kappa2,
        lam,
        timeSigma1s,
        sigma1s,
        timeSigma2s,
        sigma2s,
        timeAlphas,
        alphas,
    )
    sr2f = ShortRate2FPDE(model)
    assert sr2f is not None
