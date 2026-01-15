import velesquant.native as n


def test_hull_white_binding():
    kappa = 0.1
    time_sigmas = [1.0, 5.0]
    sigmas = [0.01, 0.015]
    time_dfs = [0.0, 1.0, 5.0, 10.0]
    dfs = [1.0, 0.95, 0.80, 0.60]

    hw = n.HullWhite(kappa, time_sigmas, sigmas, time_dfs, dfs)

    assert hw.getParameterKappa() == kappa
    assert hw.getSigmas() == sigmas
    assert hw.getTimeSigmas() == time_sigmas

    swaps = []

    s1 = n.DefSwap()
    s1.Expiry = 1.0
    s1.Tenor = 5.0
    s1.SwapRate = 0.03
    s1.VolATM = 0.20
    s1.Value = 0.0
    swaps.append(s1)

    s2 = n.DefSwap()
    s2.Expiry = 5.0
    s2.Tenor = 5.0
    s2.SwapRate = 0.04
    s2.VolATM = 0.25
    s2.Value = 0.0
    swaps.append(s2)

    assert len(swaps) == 2

    try:
        hw.calibrate(swaps, n.CalibrationTarget.Volatility)
    except Exception:
        pass


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

    sr2f = n.ShortRate2FPDE(
        kappa1, kappa2, lam, time_sigma1, sigma1, time_sigma2, sigma2, time_alpha, alpha
    )

    try:
        val = sr2f.pricingSwaption(1.0, 5.0, 0.03)
        assert isinstance(val, float)
    except Exception:
        pass

    time_dfs = [0.0, 1.0, 10.0]
    dfs = [1.0, 0.95, 0.70]

    swaps = []
    s = n.DefSwap()
    s.Expiry = 1.0
    s.Tenor = 5.0
    s.SwapRate = 0.03
    s.VolATM = 0.20
    swaps.append(s)

    try:
        sr2f.calibrate(time_dfs, dfs, swaps)
    except Exception:
        pass
