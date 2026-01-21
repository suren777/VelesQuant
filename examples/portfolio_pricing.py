import numpy as np

from velesquant.instruments.bonds import ZeroCouponBond
from velesquant.instruments.portfolio import Portfolio
from velesquant.instruments.rates import Swaption

from velesquant.market.container import Market
from velesquant.market.curves import DiscountCurve
from velesquant.models import HullWhiteModel, SabrModel


def build_market() -> Market:
    market = Market()

    usd_curve = DiscountCurve(
        times=[0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0],
        dfs=[1.0, 0.995, 0.99, 0.97, 0.94, 0.85, 0.70],
    )
    market.add("USD", usd_curve)

    return market


def load_market_data():
    return {
        "equity": {
            "spot": 100.0,
            "forward": 100.0,
            "maturity": 1.0,
            "strikes": np.array([90.0, 95.0, 100.0, 105.0, 110.0]),
            "vols": np.array([0.25, 0.22, 0.20, 0.22, 0.25]),
        },
        "swaptions": [
            {"expiry": 1.0, "tenor": 5.0, "vol_atm": 0.30, "rate": 0.03},
            {"expiry": 2.0, "tenor": 5.0, "vol_atm": 0.28, "rate": 0.035},
            {"expiry": 5.0, "tenor": 5.0, "vol_atm": 0.25, "rate": 0.04},
        ],
    }


def calibrate_sabr(market_data) -> SabrModel:
    print("\n--- Step 1: Calibrating SABR Volatility Model ---")

    eq = market_data["equity"]

    sabr = SabrModel(
        maturity=eq["maturity"],
        forward=eq["forward"],
        beta=0.5,
        alpha=0.3,
        nu=0.4,
        rho=-0.3,
    )

    sabr.calibrate(
        strikes=eq["strikes"],
        quotes=eq["vols"],
        calibration_target="Volatility",
    )

    print(f"  Forward: {eq['forward']}")
    print(f"  Calibrated alpha: {sabr.alpha:.4f}")
    print(f"  Calibrated nu:    {sabr.nu:.4f}")
    print(f"  Calibrated rho:   {sabr.rho:.4f}")

    return sabr


def calibrate_hull_white(market: Market, market_data) -> HullWhiteModel:
    print("\n--- Step 2: Creating Hull-White Interest Rate Model ---")

    curve = market.get("USD", DiscountCurve)
    _ = market_data["swaptions"]

    hw_model = HullWhiteModel(kappa=0.05, sigma=0.01)

    print(f"  Kappa: {hw_model.kappa:.4f}")
    print(f"  Sigma: {hw_model.sigma:.4f}")
    print(f"  Curve: {len(curve.times)} points")

    return hw_model


def price_portfolio(
    sabr: SabrModel, hw_model: HullWhiteModel, market: Market, market_data
):
    print("\n--- Step 3: Pricing Portfolio ---")

    portfolio_items = []

    atm_strike = market_data["equity"]["forward"]
    atm_vol = sabr.implied_vol(atm_strike)
    portfolio_items.append(("ATM Implied Vol", atm_vol))
    print(f"  ATM Implied Vol: {atm_vol:.4f}")

    otm_strike = 90.0
    otm_vol = sabr.implied_vol(otm_strike)
    portfolio_items.append(("90 Strike Implied Vol", otm_vol))
    print(f"  90 Strike Vol:   {otm_vol:.4f}")

    swaption = Swaption(expiry=1.0, tenor=5.0, strike=0.03)
    swaption_price = hw_model.price(swaption, market_data=market, curve_name="USD")
    portfolio_items.append(("1Y5Y Swaption", swaption_price))
    print(f"  1Y5Y Swaption:   {swaption_price:.6f}")

    pf = Portfolio()
    pf.add(ZeroCouponBond(maturity=2.0))
    pf.add(ZeroCouponBond(maturity=5.0))
    pf.add(Swaption(expiry=2.0, tenor=5.0, strike=0.035))

    portfolio_npv = hw_model.price(pf, market_data=market, curve_name="USD")
    portfolio_items.append(("Bond+Swaption Portfolio", portfolio_npv))
    print(f"  Bond+Swaption Portfolio: {portfolio_npv:.6f}")

    return portfolio_items


def main():
    print("=" * 60)
    print("VelesQuant Portfolio Pricing Example (Using Market Container)")
    print("=" * 60)

    market = build_market()
    print("\n  Market built with USD discount curve")

    market_data = load_market_data()

    sabr = calibrate_sabr(market_data)
    hw_model = calibrate_hull_white(market, market_data)

    portfolio_items = price_portfolio(sabr, hw_model, market, market_data)

    print("\n" + "=" * 60)
    print("Portfolio Summary")
    print("=" * 60)
    total_value = 0.0
    for name, price in portfolio_items:
        print(f"  {name:25s}: {price:12.6f}")
        total_value += price
    print("-" * 60)
    print(f"  {'Total':25s}: {total_value:12.6f}")
    print("=" * 60)

    print("\n--- Serialization Demo ---")
    import json

    market_json = json.dumps(market.to_dict(), indent=2)
    print(f"  Market serialized to JSON ({len(market_json)} bytes)")
    model_json = json.dumps(hw_model.to_dict(), indent=2)
    print(f"  HullWhiteModel: {model_json}")


if __name__ == "__main__":
    main()
