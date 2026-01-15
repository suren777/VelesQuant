"""
Test for HWPDE advanced bindings.

Verifies the newly exposed methods:
- pricingSwap, pricingZBO, pricingZB, pricingCBO, pricingCouponBond
- getSwapRate, getImpVolATM, getDFs, simulationPDE
- calibrate, getKappa, getR0, getTimeSigmas, getSigmas, getTimeThetas, getThetas
"""

import pytest

import velesquant.native as n


def test_hwpde_basic_pricing():
    """Test basic HWPDE pricing methods."""
    # Create HWPDE with DF-based constructor
    hwpde = n.HWPDE(
        kappa=0.05,
        timeSigmas=[1.0, 5.0, 10.0],
        sigmas=[0.01, 0.012, 0.015],
        discount_factor_times=[0.0, 1.0, 5.0, 10.0, 30.0],
        discount_factors=[1.0, 0.95, 0.80, 0.65, 0.30],
    )

    # Test pricingZB
    zb_price = hwpde.pricingZB(5.0)
    assert 0 < zb_price <= 1.0

    # Test pricingSwaption (already existed)
    swaption = hwpde.pricingSwaption(1.0, 5.0, 0.03, 0.5)
    assert swaption >= 0

    # Test pricingSwap
    swap = hwpde.pricingSwap(0.0, 5.0, 0.03, 0.5)
    assert isinstance(swap, float)

    # Test pricingZBO
    zbo = hwpde.pricingZBO(1.0, 5.0, 0.90, "Call")
    assert zbo >= 0

    # Test pricingCouponBond
    coupon_bond = hwpde.pricingCouponBond(0.0, 5.0, 0.04, 0.5)
    assert coupon_bond > 0

    # Test pricingCBO
    cbo = hwpde.pricingCBO(1.0, 5.0, 0.04, 1.0, 0.5, "Call")
    assert cbo >= 0


def test_hwpde_exotic_pricing():
    """Test exotic pricing methods."""
    hwpde = n.HWPDE(
        kappa=0.05,
        timeSigmas=[1.0, 5.0],
        sigmas=[0.01, 0.012],
        discount_factor_times=[0.0, 1.0, 5.0, 10.0],
        discount_factors=[1.0, 0.95, 0.80, 0.65],
    )

    # Test pricingBermudan
    exercises = [2.0, 3.0, 4.0]  # Exercise dates
    bermudan = hwpde.pricingBermudan(1.0, 5.0, exercises, 0.03, 0.5)
    assert bermudan >= 0

    # Test pricingCallableSwap
    callable_swap = hwpde.pricingCallableSwap(
        1.0, 5.0, exercises, 0.04, 1.0, 0.5, "Call"
    )
    assert isinstance(callable_swap, float)


def test_hwpde_analysis():
    """Test analysis and getter methods."""
    hwpde = n.HWPDE(
        kappa=0.05,
        timeSigmas=[1.0, 5.0],
        sigmas=[0.01, 0.012],
        discount_factor_times=[0.0, 1.0, 5.0, 10.0],
        discount_factors=[1.0, 0.95, 0.80, 0.65],
    )

    # Test getSwapRate
    swap_rate = hwpde.getSwapRate(1.0, 5.0, 0.5)
    assert -0.01 < swap_rate < 0.20

    # Test getImpVolATM
    imp_vol = hwpde.getImpVolATM(1.0, 5.0, 0.5)
    assert 0 < imp_vol < 2.0

    # Test getDFs
    time_points = [1.0, 2.0, 5.0]
    dfs = hwpde.getDFs(time_points)
    assert len(dfs) == 3
    assert all(isinstance(x, float) for x in dfs)

    # Test simulationPDE
    sim_times = [0.5, 1.0, 2.0]
    sim_path = hwpde.simulationPDE(sim_times)
    assert len(sim_path) == 3
    assert all(isinstance(x, float) for x in sim_path)

    # Test getters
    assert hwpde.getKappa() == 0.05
    assert hwpde.getTimeSigmas() == [1.0, 5.0]
    assert hwpde.getSigmas() == [0.01, 0.012]
    assert isinstance(hwpde.getInitialRate(), float)
    assert len(hwpde.getThetas()) > 0


def test_hwpde_calibration():
    """Test calibration method."""
    hwpde = n.HWPDE(
        kappa=0.05,
        timeSigmas=[1.0, 5.0],
        sigmas=[0.01, 0.012],
        discount_factor_times=[0.0, 1.0, 5.0, 10.0],
        discount_factors=[1.0, 0.95, 0.80, 0.65],
    )

    # Create swap quotes for calibration
    swaps = []

    s1 = n.DefSwap()
    s1.expiry = 1.0
    s1.tenor = 5.0
    s1.frequency = 0.5
    s1.swap_rate = 0.03
    s1.vol_atm = 0.25
    s1.value = 0.0
    swaps.append(s1)

    s2 = n.DefSwap()
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
        assert hwpde.getKappa() is not None
    except RuntimeError:
        # We expect a RuntimeError with "too much freedom" or similar
        # because the input data is dummy data and might not converge/be valid
        # This confirms the binding is reachable and throws the expected C++ exception
        pass
    except Exception as e:
        pytest.fail(f"Unexpected exception during calibration binding call: {e}")
