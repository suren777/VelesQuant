import pytest

from velesquant.models.swaption_sabr import SabrSwaptionModel


def test_swaption_init():
    s = SabrSwaptionModel(
        expiry=1.0,
        tenor=5.0,
        forward=0.03,
        annuity=1.0,
        beta=0.5,
        alpha=0.2,
        nu=0.3,
        rho=-0.4,
    )

    assert s.alpha == pytest.approx(0.2)
    assert s.nu == pytest.approx(0.3)
    assert s.rho == pytest.approx(-0.4)


def test_swaption_defaults():
    s = SabrSwaptionModel(expiry=1.0, tenor=5.0, forward=0.03, annuity=1.0)

    # Defaults: beta=0.85, alpha=0.5, nu=0.25, rho=-0.75
    assert s.alpha == pytest.approx(0.5)
    assert s.nu == pytest.approx(0.25)
    assert s.rho == pytest.approx(-0.75)


def test_swaption_pricing():
    s = SabrSwaptionModel(
        expiry=1.0,
        tenor=5.0,
        forward=0.03,
        annuity=1.0,
        beta=0.85,
        alpha=0.2,
        nu=0.3,
        rho=-0.4,
    )

    price = s.fair_value(0.03, "Call")
    assert price >= 0.0

    vol = s.implied_vol(0.03)
    assert vol > 0.0


def test_swap_fair_value():
    s = SabrSwaptionModel(expiry=1.0, tenor=5.0, forward=0.03, annuity=1.0)

    # Forward is 0.03
    # Strike 0.03 -> Swap Value should be 0 (Forward - Strike) * Annuity
    val_atm = s.swap_value(0.03)
    assert val_atm == pytest.approx(0.0)

    # Strike 0.02 -> (0.03 - 0.02) * 1.0 = 0.01
    val_itm = s.swap_value(0.02)
    assert val_itm == pytest.approx(0.01)
