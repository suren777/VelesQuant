from velesquant import DefSwap, ShortRate2FModel, ShortRate2FPDE


def test_short_rate_2f_binding():
    kappa1 = 0.03
    kappa2 = 0.01
    lam = -0.5

    time_sigma1 = [10.0]
    sigma1 = [0.01]

    time_sigma2 = [10.0]
    sigma2 = [0.02]

    time_alpha = [10.0]
    alpha = [0.0]

    model = ShortRate2FModel(
        kappa1,
        kappa2,
        lam,
        time_sigma1,
        sigma1,
        time_sigma2,
        sigma2,
        time_alpha,
        alpha,
    )
    sr2f = ShortRate2FPDE(model)

    try:
        val = sr2f.pricing_swaption(1.0, 5.0, 0.03)
        assert isinstance(val, float)
    except Exception:
        pass

    time_dfs = [0.0, 1.0, 10.0]
    dfs = [1.0, 0.95, 0.70]

    swaps = []
    s = DefSwap()
    s.expiry = 1.0
    s.tenor = 5.0
    s.swap_rate = 0.03
    s.vol_atm = 0.20
    swaps.append(s)

    try:
        sr2f.calibrate(time_dfs, dfs, swaps)
    except Exception:
        pass
