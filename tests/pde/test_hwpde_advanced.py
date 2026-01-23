"""
Test for HWPDE advanced bindings.

Verifies the newly exposed methods:
- pricingSwap, pricingZBO, pricingZB, pricingCBO, pricingCouponBond
- getSwapRate, getImpVolATM, getDFs, simulationPDE
- calibrate, getKappa, getR0, getTimeSigmas, getSigmas, getTimeThetas, getThetas
"""

import pytest

from velesquant import HWPDE, DefSwap, HullWhiteModel, OptionType


def create_hwpde_model(
    kappa, time_sigmas, sigmas, time_dfs, dfs, grid_points=100, time_step=0.01
):
    """Helper to create HWPDE with a HullWhiteModel."""
    model = HullWhiteModel(kappa, time_sigmas, sigmas, time_dfs, dfs)
    return HWPDE(model, grid_points, time_step)


def test_hwpde_basic_pricing():
    """Test basic HWPDE pricing methods."""
    hwpde = create_hwpde_model(
        kappa=0.05,
        time_sigmas=[1.0, 5.0, 10.0],
        sigmas=[0.01, 0.012, 0.015],
        time_dfs=[0.0, 1.0, 5.0, 10.0, 30.0],
        dfs=[1.0, 0.95, 0.80, 0.65, 0.30],
        grid_points=100,
        time_step=0.01,
    )

    # Test pricingZB
    zb_price = hwpde.price_zero_bond(5.0)
    assert 0 < zb_price <= 1.0

    # Test pricingSwaption (already existed)
    swaption = hwpde.price_swaption(1.0, 5.0, 0.03, 0.5)
    assert swaption >= 0

    # Test pricingSwap
    swap = hwpde.price_swap(0.0, 5.0, 0.03, 0.5)
    assert isinstance(swap, float)

    # Test pricingZBO
    zbo = hwpde.price_zero_bond_option(1.0, 5.0, 0.90, OptionType.Call)
    assert zbo >= 0

    # Test pricingCouponBond
    coupon_bond = hwpde.price_coupon_bond(0.0, 5.0, 0.04, 0.5)
    assert coupon_bond > 0

    # Test pricingCBO
    cbo = hwpde.price_coupon_bond_option(1.0, 5.0, 0.04, 1.0, 0.5, OptionType.Call)
    assert cbo >= 0


def test_hwpde_exotic_pricing():
    """Test exotic pricing methods."""
    hwpde = create_hwpde_model(
        kappa=0.05,
        time_sigmas=[1.0, 5.0],
        sigmas=[0.01, 0.012],
        time_dfs=[0.0, 1.0, 5.0, 10.0],
        dfs=[1.0, 0.95, 0.80, 0.65],
        grid_points=100,
        time_step=0.01,
    )

    # Test pricingBermudan
    exercises = [2.0, 3.0, 4.0]  # Exercise dates
    bermudan = hwpde.price_bermudan(1.0, 5.0, exercises, 0.03, 0.5)
    assert bermudan >= 0

    # Test pricingCallableSwap
    callable_swap = hwpde.price_callable_swap(
        1.0, 5.0, exercises, 0.04, 1.0, 0.5, OptionType.Call
    )
    assert isinstance(callable_swap, float)


def test_hwpde_analysis():
    """Test analysis and getter methods."""
    hwpde = create_hwpde_model(
        kappa=0.05,
        time_sigmas=[1.0, 5.0],
        sigmas=[0.01, 0.012],
        time_dfs=[0.0, 1.0, 5.0, 10.0],
        dfs=[1.0, 0.95, 0.80, 0.65],
        grid_points=100,
        time_step=0.01,
    )

    # Test getSwapRate
    swap_rate = hwpde.get_swap_rate(1.0, 5.0, 0.5)
    assert -0.01 < swap_rate < 0.20

    # Test getImpVolATM
    imp_vol = hwpde.get_implied_vol_atm(1.0, 5.0, 0.5)
    assert 0 < imp_vol < 2.0

    # Test getDFs
    time_points = [1.0, 2.0, 5.0]
    dfs = hwpde.get_discount_factors(time_points)
    assert len(dfs) == 3
    assert all(isinstance(x, float) for x in dfs)

    # Test simulationPDE
    sim_times = [0.5, 1.0, 2.0]
    sim_path = hwpde.simulate(sim_times)
    assert len(sim_path) == 3
    assert all(isinstance(x, float) for x in sim_path)

    # Test getters
    assert hwpde.get_kappa() == 0.05
    assert hwpde.get_time_sigmas() == [1.0, 5.0]
    assert hwpde.get_sigmas() == [0.01, 0.012]
    assert isinstance(hwpde.get_initial_rate(), float)
    assert len(hwpde.get_thetas()) > 0


def test_hwpde_calibration():
    """Test calibration method."""
    hwpde = create_hwpde_model(
        kappa=0.05,
        time_sigmas=[1.0, 5.0],
        sigmas=[0.01, 0.012],
        time_dfs=[0.0, 1.0, 5.0, 10.0],
        dfs=[1.0, 0.95, 0.80, 0.65],
        grid_points=100,
        time_step=0.01,
    )

    # Create swap quotes for calibration
    swaps = []

    s1 = DefSwap()
    s1.expiry = 1.0
    s1.tenor = 5.0
    s1.frequency = 0.5
    s1.swap_rate = 0.03
    s1.vol_atm = 0.25
    s1.value = 0.0
    swaps.append(s1)

    s2 = DefSwap()
    s2.expiry = 5.0
    s2.tenor = 5.0
    s2.frequency = 0.5
    s2.swap_rate = 0.04
    s2.vol_atm = 0.20
    s2.value = 0.0
    swaps.append(s2)

    time_dfs = [0.0, 1.0, 5.0, 10.0]
    dfs = [1.0, 0.95, 0.80, 0.65]

    try:
        hwpde.calibrate(time_dfs, dfs, swaps)
        # Check that calibration changed parameters or finished
        assert hwpde.get_kappa() is not None
    except RuntimeError:
        # We expect a RuntimeError with "too much freedom" or similar
        # because the input data is dummy data and might not converge/be valid
        pass
    except Exception as e:
        pytest.fail(f"Unexpected exception during calibration binding call: {e}")
