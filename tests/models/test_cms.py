import pytest

from velesquant.models.cms import CMSModel


def test_cms_init():
    cms = CMSModel(
        expiry=1.0,
        tenor=5.0,
        forward=0.03,
        annuity=1.0,
        pay_cms=1.0,
        discount_cms=0.95,
        beta=0.5,
        alpha=0.2,
        nu=0.3,
        rho=-0.4,
    )

    # getForward returns the convexity-adjusted forward, so it should NOT equal forwardSR exactly
    assert cms.adjusted_forward != pytest.approx(0.03)
    assert cms.discount_cms == pytest.approx(0.95)
    assert cms.maturity == pytest.approx(1.0)
    assert cms.alpha == pytest.approx(0.2)


def test_cms_defaults():
    cms = CMSModel(
        expiry=1.0, tenor=5.0, forward=0.03, annuity=1.0, pay_cms=1.0, discount_cms=0.95
    )

    # Defaults: beta=0.85, alpha=0.5, nu=0.25, rho=-0.75
    assert cms.alpha == pytest.approx(0.5)


def test_cms_fair_value():
    cms = CMSModel(
        expiry=1.0,
        tenor=5.0,
        forward=0.03,
        annuity=1.0,
        pay_cms=1.0,
        discount_cms=0.95,
        beta=0.85,
        alpha=0.2,
        nu=0.25,
        rho=-0.5,
    )

    val = cms.fair_value(0.03, "Call")
    assert val >= 0.0

    vol = cms.implied_vol(0.03)
    assert vol > 0.0
