import math

import pytest

import velesquant.native as m


def test_utilities():
    assert m.cdf_normal(0.0) == pytest.approx(0.5)
    assert m.pdf_normal(0.0) == pytest.approx(1.0 / math.sqrt(2 * math.pi))

    # Test implied vol
    # Call BS for price, then back it out
    spot, strike, r, d, vol, T = 100.0, 100.0, 0.05, 0.02, 0.2, 1.0
    price = m.black_scholes_call(spot, strike, r, d, vol, T)
    _ = m.implied_vol(
        T, spot * math.exp((r - d) * T), strike, price
    )  # black formula uses forward
    # wait, BS takes spot. black_formula_call takes Forward.
    fwd = spot * math.exp((r - d) * T)
    price_black = m.black_formula_call(fwd, strike, vol, T)
    iv_black = m.implied_vol(T, fwd, strike, price_black)
    assert iv_black == pytest.approx(vol, rel=1e-5)


def test_equity_models_basic():
    spot, strike, r, d, vol, T = 100.0, 105.0, 0.03, 0.01, 0.25, 0.5
    price = m.black_scholes_call(spot, strike, r, d, vol, T)
    assert price > 0
    vega = m.black_scholes_call_vega(spot, strike, r, d, vol, T)
    assert vega > 0
