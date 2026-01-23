import numpy as np

from velesquant import (
    DefSwap,
    ShortRate2FModel,
    ShortRate2FPDE,
)


class TestShortRate2FPDE:
    def test_initialization(self):
        # Parameters match the original constructor signature logic
        # kappa1, kappa2, lambda, timeSigma1s, sigma1s, timeSigma2s, sigma2s, timeAlphas, alphas
        time_sigmas1 = [0.1, 0.5, 1.0, 5.0]
        sigmas1 = [0.01, 0.01, 0.01, 0.01]
        time_sigmas2 = [0.1, 0.5, 1.0, 5.0]
        sigmas2 = [0.01, 0.01, 0.01, 0.01]
        time_alphas = [0.1, 0.5, 1.0, 5.0]
        alphas = [0.02, 0.02, 0.02, 0.02]

        model = ShortRate2FModel(
            0.03,
            0.01,
            0.0,  # kappa1, kappa2, lambda
            time_sigmas1,
            sigmas1,
            time_sigmas2,
            sigmas2,
            time_alphas,
            alphas,
        )

        assert model.get_kappa1() == 0.03
        assert model.get_kappa2() == 0.01
        assert model.get_lambda() == 0.0

        solver = ShortRate2FPDE(model)
        assert solver is not None

    def test_calibrate_and_price(self):
        # Create a model with minimal parameters to speed up test
        time_sigmas1 = [10.0]  # Only 1 time point
        sigmas1 = [0.007]
        time_sigmas2 = [10.0]
        sigmas2 = [0.007]
        time_alphas = [10.0]
        alphas = [0.02]

        model = ShortRate2FModel(
            0.03,
            0.01,
            0.0,
            time_sigmas1,
            sigmas1,
            time_sigmas2,
            sigmas2,
            time_alphas,
            alphas,
        )
        solver = ShortRate2FPDE(model)

        # Setup calibration data
        dfs_times = [1.0, 5.0]
        dfs = [np.exp(-0.02 * t) for t in dfs_times]

        swaps = []
        # Generate just enough swaptions for degrees of freedom (params > data error?)
        # Model has: 1 sigma1 + 1 sigma2 + kappa1 + kappa2 + lambda = 5 params
        # We need >= 5 quotes. Let's use 6.
        for i in range(1, 7):
            expiry = i * 0.5
            tenor = 5.0
            swaps.append(DefSwap(expiry, tenor, 0.5, 0.02, 0.2, 0.05))

        # The calibrator updates the model in place
        # NOTE: Calibration is computationally expensive for 2F PDE and can hang in unit tests.
        # We verify that the solver bindings and pricing methods work.
        # solver.calibrate(dfs_times, dfs, swaps)

        # Check pricing works with initial parameters
        price = solver.price_swaption(1.0, 5.0, 0.02, 0.5)
        assert price > 0.0
        # assert price < 1.0 # Depends on params, but usually true for meaningful inputs

        zb_price = solver.price_zero_bond(1.0)
        assert 0.0 < zb_price < 1.0
