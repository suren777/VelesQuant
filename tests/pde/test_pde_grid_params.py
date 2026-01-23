from velesquant.models.pde_solvers import (
    HWPDEModel,
    ShortRate1FPDEModel,
    ShortRate2FPDEModel,
)


def test_hwpde_grid_params():
    kappa = 0.1
    time_sigmas = [1.0, 2.0]
    sigmas = [0.01, 0.012]
    discount_curve = [(0.0, 1.0), (1.0, 0.95)]

    # Test with default params
    model_default = HWPDEModel(
        kappa=kappa,
        time_sigmas=time_sigmas,
        sigmas=sigmas,
        discount_curve=discount_curve,
    )
    assert model_default is not None

    # Test with custom grid params
    model_custom = HWPDEModel(
        kappa=kappa,
        time_sigmas=time_sigmas,
        sigmas=sigmas,
        discount_curve=discount_curve,
        grid_points=100,
        time_step=0.01,
    )
    assert model_custom is not None

    # Simple pricing check to ensure grid setup works
    price = model_custom.price_zero_bond(maturity=1.0)
    assert price > 0.0


def test_shortrate1f_grid_params():
    initial_rate = 0.03
    kappa = 0.1
    alpha = 0.01
    beta = 0.5
    gamma = 0.0
    time_sigmas = [1.0]
    sigmas = [0.01]

    # Test with default params
    model_default = ShortRate1FPDEModel(
        initial_rate, kappa, alpha, beta, gamma, time_sigmas, sigmas
    )
    assert model_default is not None

    # Test with custom grid params
    model_custom = ShortRate1FPDEModel(
        initial_rate,
        kappa,
        alpha,
        beta,
        gamma,
        time_sigmas,
        sigmas,
        grid_points=100,
        time_step=0.01,
    )
    assert model_custom is not None

    # Simple pricing check
    price = model_custom.price_zbo(expiry=1.0, maturity=2.0, strike=0.9)
    assert price > 0.0


def test_shortrate2f_grid_params():
    kappa1 = 0.03
    kappa2 = 0.01
    lam = -0.7
    time_sigma1s = [1.0]
    sigma1s = [0.01]
    time_sigma2s = [1.0]
    sigma2s = [0.02]
    time_alphas = [1.0]
    alphas = [0.0]

    # Test with default params
    model_default = ShortRate2FPDEModel(
        kappa1,
        kappa2,
        lam,
        time_sigma1s,
        sigma1s,
        time_sigma2s,
        sigma2s,
        time_alphas,
        alphas,
    )
    assert model_default is not None

    # Test with custom grid params
    model_custom = ShortRate2FPDEModel(
        kappa1,
        kappa2,
        lam,
        time_sigma1s,
        sigma1s,
        time_sigma2s,
        sigma2s,
        time_alphas,
        alphas,
        time_step=0.1,
    )
    assert model_custom is not None

    # Simple pricing check
    price = model_custom.price_zero_bond(maturity=1.0)
    assert price > 0.0


def test_calibration_optimizer_params():
    # Setup HWPDE for calibration
    kappa = 0.1
    time_sigmas = [1.0]
    sigmas = [0.01]
    discount_curve = [(0.0, 1.0), (1.0, 0.95)]
    model = HWPDEModel(kappa, time_sigmas, sigmas, discount_curve=discount_curve)

    time_dfs = [0.0, 1.0, 2.0, 4.0]
    dfs = [1.0, 0.95, 0.90, 0.85]
    swap_quotes = [
        {"expiry": 1.0, "tenor": 1.0, "rate": 0.05, "vol": 0.2},
        {"expiry": 1.0, "tenor": 0.5, "rate": 0.05, "vol": 0.2},
        {"expiry": 2.0, "tenor": 1.0, "rate": 0.05, "vol": 0.2},
    ]

    optimizer_params = {"maxfev": 10000, "ftol": 1e-5}

    # Should run without error
    model.calibrate(time_dfs, dfs, swap_quotes, optimizer_params)
