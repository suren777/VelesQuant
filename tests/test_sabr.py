import pytest
import velesquant._core as m
import math


def test_sabr_initialization():
    s = m.Sabr(1.0, 0.05, 0.5, 0.2, 0.3, -0.2)
    assert s.getMaturity() == 1.0
    assert s.getForward() == 0.05
    assert s.getBeta() == 0.5
    assert s.alpha == 0.2
    assert s.nu == 0.3
    assert s.rho == -0.2


def test_sabr_default_initialization():
    s = m.Sabr()
    # Check what defaults are if any (the default constructor was empty in C++, likely uninitialized garbage or default member init?)
    # Looking at sabr.h: sabr() {};
    # Members: maturity_, forward_ etc. are NOT initialized in default ctor.
    # This might be dangerous. But let's check if it crashes using bindings.
    pass


def test_sabr_atm_vol():
    # ATM: Strike = Forward = 0.05
    forward = 0.05
    s = m.Sabr(1.0, forward, 0.5, 0.2, 0.0, 0.0)
    # With nu=0, rho=0, it simplifies?
    iv = s.impliedVol(forward)
    assert iv > 0


def test_sabr_property_updates():
    s = m.Sabr(1.0, 0.05, 0.5, 0.2, 0.3, -0.2)
    original_iv = s.impliedVol(0.05)

    s.alpha = 0.4
    new_iv = s.impliedVol(0.05)
    assert new_iv != original_iv


def test_sabr_negative_strike():
    s = m.Sabr(1.0, 0.05)
    # SABR usually handles positive strikes.
    # Just checking it doesn't segfault.
    try:
        s.impliedVol(-0.01)
    except:
        pass
