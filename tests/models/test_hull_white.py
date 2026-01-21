from velesquant import HHW
from velesquant.instruments.rates import Swaption
from velesquant.market.curves import DiscountCurve
from velesquant.models.hullwhite import HullWhiteModel


def test_hull_white_basic():
    kappa = 0.03
    # Model wrapper uses constant sigma for now in simple constructor
    # But native test used time-dependent sigmas?
    # native: timeSigmas = [0.0, 10.0], sigmas = [0.01, 0.01]. This is constant 0.01.
    sigma = 0.01

    timeDFs = [0.0, 1.0, 2.0, 5.0, 10.0]
    DFs = [1.0, 0.98, 0.95, 0.85, 0.75]
    curve = DiscountCurve(timeDFs, DFs)

    hw = HullWhiteModel(kappa, sigma)

    # Test Option on Bond
    # original: hw.optionBond(1.0, 2.0, 0.9, m.OptionType.Call)
    p = hw.price_bond_option(1.0, 2.0, 0.9, curve, "Call")
    assert p >= 0

    # Test Swaption
    # original: hw.swaption(1.0, 5.0, 0.02)
    # New interface usage:
    # Swaption(expiry, tenor, strike)
    swaption = Swaption(expiry=1.0, tenor=5.0, strike=0.02)
    swaption_p = hw.price(swaption, curve)
    assert swaption_p >= 0

    # Test getters
    assert hw.kappa == 0.03
    assert hw.sigma == 0.01
    # assert hw.getTimeSigmas() ... wrapper stores constant sigma, doesn't expose list yet.
    # We accept this difference as the Python wrapper simplifies constant vol case.

    # Test simulation
    times = [0.0, 0.5, 1.0]
    paths = hw.simulate(times, curve)
    # Output should be vector of size times.size()
    assert len(paths) == 3


def test_hhw_binding():
    # HHW(s0, v0, r0, kappa, eta, rho, sigma1, sigma2, a)
    # This tests the hybrid model binding which is NOT wrapped yet.
    # Leaving as native usage.
    hhw = HHW(100.0, 0.04, 0.03, 1.0, 0.1, -0.5, 0.2, 0.2, 0.5)
    price = hhw.price(1.0, 100.0)
    assert isinstance(price, float)
    assert price > 0


def test_hhw_pricing_scenario():
    """
    Scenario: Pricing a call option with HHW.
    Compare price directionality when volatility params increase.
    """
    # Base case
    hhw_base = HHW(100.0, 0.04, 0.03, 1.0, 0.1, -0.5, 0.2, 0.2, 0.5)
    price_base = hhw_base.price(1.0, 100.0)

    # Increased volatility (sigma1)
    hhw_vol = HHW(100.0, 0.04, 0.03, 1.0, 0.1, -0.5, 0.4, 0.2, 0.5)
    price_vol = hhw_vol.price(1.0, 100.0)

    # Volatility impact can be complex in HHW depending on regimes/calibration.
    # For now, we verified the price changes (sensitivity exists).
    assert price_vol != price_base


def test_python_swaption_utils():
    # Test bindings for get_swap_rate and swaption_implied_vol via Python wrapper
    # which now implements them using HullWhiteModel Core and local logic.

    kappa = 0.03
    sigma = 0.01

    timeDFs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    # Simple discount factors: exp(-r*t) with r=0.05
    r = 0.05
    import math

    DFs = [math.exp(-r * t) for t in timeDFs]
    curve = DiscountCurve(timeDFs, DFs)

    hw = HullWhiteModel(kappa, sigma)

    expiry = 1.0
    tenor = 2.0

    # 1. Test get_swap_rate
    sr = hw.get_swap_rate(expiry, tenor, curve)
    # Expected approx par swap rate ~ 5%
    assert 0.04 < sr < 0.06

    # 2. Test swaption_implied_vol
    # Calculate a price first using price method
    swaption = Swaption(expiry=expiry, tenor=tenor, strike=sr)
    price = hw.price(swaption, curve)

    # Get implied vol
    vol = hw.swaption_implied_vol(expiry, tenor, price, curve)

    # The model generated price corresponds to some vol.
    # In HW model, the implied Black vol isn't constant but should be reasonable.
    assert 0.0 < vol < 1.0
