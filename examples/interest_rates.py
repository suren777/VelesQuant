""" """

import numpy as np

# Instruments
from velesquant.instruments.bonds import ZeroCouponBond
from velesquant.instruments.rates import Swaption

# Market Data & Containers
from velesquant.market.container import Market
from velesquant.market.curves import DiscountCurve

# Models
from velesquant.models import CMSModel, CMSSpreadModel, HullWhiteModel, HWPDEModel


def build_market() -> Market:
    """Build a Market container with a USD Discount Curve."""
    market = Market()

    # Market Construction
    times = np.linspace(0.0, 30.0, 31)
    rates = 0.03 * np.ones_like(times)
    dfs = np.exp(-rates * times)

    usd_curve = DiscountCurve(times=times.tolist(), dfs=dfs.tolist())
    market.add("USD", usd_curve)
    return market


def run_hull_white_examples(market: Market):
    print("\n" + "=" * 50)
    print(" Hull-White Model (Analytic & Tree)")
    print("=" * 50)

    hw = HullWhiteModel(kappa=0.1, sigma=0.01)
    print(f"Model Initialized: Kappa={hw.kappa}, Sigma={hw.sigma}")

    swaption = Swaption(expiry=5.0, tenor=10.0, strike=0.03, pay_frequency=0.5)
    price = hw.price(swaption, market_data=market, curve_name="USD")
    print(f"European Swaption (1Y into 10Y) Price: {price:.6f}")

    zcb = ZeroCouponBond(maturity=10.0)
    zcb_price = hw.price(zcb, market_data=market, curve_name="USD")
    print(f"Zero Coupon Bond (10Y) Price (Model):  {zcb_price:.6f}")


def run_pde_examples(market: Market):
    print("\n" + "=" * 50)
    print(" PDE Solver (Hull-White Framework)")
    print("=" * 50)

    curve = market.get("USD", DiscountCurve)
    times = curve.times
    dfs = curve.dfs

    hwpde = HWPDEModel(
        kappa=0.1,
        time_sigmas=[0.0, 30.0],
        sigmas=[0.01, 0.01],
        discount_factor_times=times,
        discount_factors=dfs,
    )

    exercises = [1.0, 2.0, 3.0, 4.0]
    expiry = 1.0
    tenor = 4.0
    strike = 0.03
    try:
        berm_price = hwpde.price_bermudan(
            expiry=expiry,
            tenor=tenor,
            exercises=exercises,
            strike=strike,
            pay_freq=0.5,
        )
        print(f"Bermudan Swaption Price: {berm_price:.6f}")
    except Exception as e:
        print(f"Bermudan Pricing skipped: {e}")


def run_cms_examples(market: Market):
    print("\n" + "=" * 50)
    print(" Constant Maturity Swap (CMS) Products")
    print("=" * 50)

    print("Pricing CMS Cap/Floor...")

    fwd_rate = 0.035
    annuity = 8.5
    expiry = 5.0
    tenor = 10.0

    cms_model = CMSModel(
        expiry=expiry,
        tenor=tenor,
        forward=fwd_rate,
        annuity=annuity,
        pay_cms=expiry,
        discount_cms=0.85,
        beta=0.5,
        alpha=0.2,
        nu=0.4,
        rho=-0.2,
    )

    call_price = cms_model.fair_value(strike=0.04, option_type="Call")
    put_price = cms_model.fair_value(strike=0.03, option_type="Put")

    print(f"CMS Call (K=4%): {call_price:.6f}")
    print(f"CMS Put (K=3%):  {put_price:.6f}")

    print("\nPricing CMS Spread Option (10Y - 2Y)...")
    try:
        spread_model = CMSSpreadModel(
            expiry1=5.0,
            tenor1=10.0,
            fwd1=0.035,
            annuity1=8.5,
            pay1=5.0,
            disc1=0.85,
            beta1=0.5,
            strikes1=np.array([0.035]),
            quotes1=np.array([0.20]),
            type1="Volatility",
            expiry2=5.0,
            tenor2=2.0,
            fwd2=0.025,
            annuity2=1.9,
            pay2=5.0,
            disc2=0.85,
            beta2=0.5,
            strikes2=np.array([0.025]),
            quotes2=np.array([0.25]),
            type2="Volatility",
            corr=0.6,
        )

        print("CMS Spread Model initialized successfully.")
    except Exception as e:
        print(f"CMS Spread setup skipped: {e}")


def main():
    print("VelesQuant Interest Rate Examples")
    market = build_market()

    run_hull_white_examples(market)
    run_pde_examples(market)
    run_cms_examples(market)


if __name__ == "__main__":
    main()
