# Supported Instruments in VelesQuant

This document outlines the financial instruments that the **VelesQuant** library can price, based on the current codebase analysis (Python wrappers, Native bindings, and Pricing Engines).

## Interest Rate Products

These instruments are primarily priced using the `HullWhiteModel` or `HWPDE` engines.

| Instrument | Description | Supported by |
| :--- | :--- | :--- |
| **Zero Coupon Bond** | A bond that pays a single face value at maturity with no coupons. | `HWPDE`, `HullWhiteModel`, `ZeroCouponBond` class |
| **Coupon Bond** | A bond that pays regular coupons and face value at maturity. | `HWPDE`, `CouponBond` class |
| **Interest Rate Swap** | A vanilla swap exchanging fixed for floating rate payments. | `HWPDE`, `Swaption` (underlying) |
| **Callable Swap** | A swap where one party has the right to terminate the contract on specified dates. | `HWPDE` |
| **DefSwap** | Defined Swap structure used for calibration. | `DefSwap` native class |

## Interest Rate Derivatives

Options on interest rate products, supported by various models including Hull-White, Short Rate models, and SABR.

| Instrument | Description | Supported by |
| :--- | :--- | :--- |
| **Swaption (European)** | Option to enter into a swap at a future date. | `HWPDE`, `HullWhiteModel`, `Swaption` class, `ShortRate1/2FPDE` |
| **Swaption (Bermudan)** | Swaption that can be exercised on multiple specified dates. | `HWPDE`, `ExerciseStyle.Bermudan` |
| **Bond Option** | Option to buy/sell a bond (Zero or Coupon). | `HWPDE` (Zero/Coupon), `HullWhiteModel`, `bond_option` |
| **CMS Cap/Floor** | Constant Maturity Swap rate option. | `CMS` class |
| **CMS Spread Option** | Option on the spread between two CMS rates. | `CmsSpread` class |
| **Quanto CMS** | CMS option with FX quanto adjustment. | `QuantoedCMS` class |
| **Quanto CMS Spread** | CMS Spread option with FX quanto adjustment. | `QuantoedCmsSpread` class |

## Equity & FX Derivatives

These generic option types are supported by volatility models like Local Volatility, Heston, and SABR.

| Instrument | Description | Supported by |
| :--- | :--- | :--- |
| **European Call/Put** | Standard option to buy/sell asset at strike. | `LocalVol`, `Heston`, `SchobelZhu`, `BlackScholes` helpers |
| **Basket Option** | Option on a basket of underlying assets. | `LogNormalBasket` (Simulation) |

## Portfolio Management

- **Portfolio**: A container class (`src/velesquant/instruments/portfolio.py`) allows grouping multiple instruments for batch pricing and risk management.

## Pricing Engines & Models

The following models are available to price the above instruments:
- **Hull-White (1F)**: Analytics, Trees, PDE
- **Short Rate (1F/2F)**: PDE
- **SABR**: Analytic, PDE
- **Heston**: Semi-analytic, Simulation
- **Local Volatility**: PDE
- **Schobel-Zhu**: Stochastic Volatility
