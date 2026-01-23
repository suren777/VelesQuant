""" """

import numpy as np

# Models
from velesquant.models import (
    HestonModel,
    LocalVolModel,
    LogNormalBasketModel,
    SabrModel,
)


def run_sabr_example():
    print("\n" + "=" * 50)
    print(" SABR Stochastic Volatility Model")
    print("=" * 50)

    maturity = 1.0
    forward = 1.20

    strikes = [1.15, 1.20, 1.25]
    vols = [0.11, 0.10, 0.115]

    print("Calibrating SABR to market smile...")
    sabr = SabrModel(
        maturity=maturity,
        forward=forward,
        beta=0.5,
        alpha=0.1,
        nu=0.4,
        rho=-0.1,
    )

    sabr.calibrate(strikes=strikes, quotes=vols, calibration_target="Volatility")

    print(
        f"Calibrated Params: Alpha={sabr.alpha:.4f}, Nu={sabr.nu:.4f}, Rho={sabr.rho:.4f}"
    )

    from velesquant import black_formula_call

    strike_price = 1.22
    vol = sabr.implied_vol(strike_price)

    call_price = black_formula_call(forward, strike_price, vol, maturity)
    print(f"Call Option (K={strike_price}) Price: {call_price:.6f}")


def run_heston_example():
    print("\n" + "=" * 50)
    print(" Heston Stochastic Volatility Model")
    print("=" * 50)

    heston = HestonModel(
        spot=100.0,
        var0=0.04,
        kappa=2.0,
        theta=0.04,
        xi=0.3,
        rho=-0.7,
        seed=12345,
    )
    print("Heston Model initialized.")

    expiry = 1.0
    strike = 100.0
    fwd = 100.0 * np.exp(0.03 * expiry)

    price = heston.price_option(
        maturity=expiry, forward=fwd, strike=strike, option_type="call"
    )
    print(f"Heston Call (ATM) Price: {price:.6f}")


def run_local_vol_example():
    print("\n" + "=" * 50)
    print(" Local Volatility Model (PDE)")
    print("=" * 50)

    s1 = SabrModel(maturity=1.0, forward=100.0, alpha=0.2, beta=0.5, nu=0.4, rho=-0.5)
    s2 = SabrModel(maturity=2.0, forward=100.0, alpha=0.2, beta=0.5, nu=0.4, rho=-0.5)

    try:
        lv = LocalVolModel(sabr_models=[s1, s2], spot=100.0)
        print("Local Volatility Surface constructed.")

        price = lv.price_european(
            maturity=1.5, strike=100.0, option_type="call", n_steps=100
        )
        print(f"Local Vol Call Price (PDE): {price:.6f}")
    except Exception as e:
        print(f"Local Vol example skipped: {e}")


def run_basket_example():
    print("\n" + "=" * 50)
    print(" Basket Option (LogNormal Simulation)")
    print("=" * 50)

    spots = [100.0, 50.0]
    vols = [[0.2], [0.3]]
    corrs = [[1.0, 0.5], [0.5, 1.0]]
    forwards = [[100.0], [50.0]]

    try:
        basket = LogNormalBasketModel(
            spots=spots,
            strikes=[100.0, 50.0],
            maturities=[1.0],
            forwards=forwards,
            ivs=vols,
            correlation=corrs,
        )

        schedule = np.linspace(0.0, 1.0, 10)

        print("Basket Model initialized (Simulation ready).")
    except Exception as e:
        print(f"Basket example skipped: {e}")


def main():
    print("VelesQuant Equity & FX Examples")

    run_sabr_example()
    run_heston_example()
    run_local_vol_example()
    run_basket_example()


if __name__ == "__main__":
    main()
