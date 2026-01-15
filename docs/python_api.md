# VelesQuant Python API Reference

This document provides a detailed overview of the core functions and classes exposed to Python from the `VelesQuant` C++ library. For practical examples and usage scenarios, please refer to the linked test files.

## Utility Functions
*Use cases: [`tests/core/test_general.py`](../tests/core/test_general.py)*

### `cdf_normal(x: float) -> float`
Calculates the cumulative distribution function for a standard normal distribution.

### `pdf_normal(x: float) -> float`
Calculates the probability density function for a standard normal distribution.

### `implied_vol(maturity: float, forward: float, strike: float, price: float) -> float`
Calculates the implied volatility using the Black formula for a given option price.

---

## Equity & Volatility Models

### Black-Scholes Functions
*Use cases: [`tests/core/test_general.py`](../tests/core/test_general.py)*
*   **`black_scholes_call(spot, strike, r, d, vol, expiry)`**: Standard Black-Scholes call price.
*   **`black_scholes_call_vega(spot, strike, r, d, vol, expiry)`**: Vega of a Black-Scholes call.

### `Sabr` Class
*Use cases: [`tests/models/test_sabr.py`](../tests/models/test_sabr.py)*
A class representing the SABR volatility model.
*   **Constructor**: `Sabr(maturity, forward, beta=0.85, alpha=0.5, nu=0.25, rho=-0.75)`
*   **Properties**: `alpha`, `nu`, `rho`, `maturity`, `forward`, `beta` (all getter/setter).
*   **Methods**:
    *   `impliedVol(strike)`: Returns ATM implied volatility.
    *   `premiumBlackScholes(strike, callORput="call")`: Standard BS premium using SABR vol.

### `Heston` Class
*Use cases: [`tests/models/test_heston.py`](../tests/models/test_heston.py)*
A class representing the Heston stochastic volatility model.
*   **Constructor**: `Heston(spot, var0, kappa, theta, xi, rho, seed=42)`
*   **Methods**:
    *   `hestonPrice(maturity, forward, strike, optType="call")`: Pricing via characteristic function.
    *   `simulationHeston(times, forwards)`: Returns a Monte Carlo simulation path.

### `SchobelZhu` Class
*Use cases: [`tests/models/test_schobzhu.py`](../tests/models/test_schobzhu.py)*
A class representing the Schobel-Zhu stochastic volatility model.
*   **Constructor**: `SchobelZhu(spot, var0, kappa, theta, xi, rho)`
*   **Methods**:
    *   `SchobelPrice(maturity, forward, strike)`: Analytical pricing.
    *   `simulationSchobelZhu(times, forwards)`: Monte Carlo simulation.
    *   `calibrator(maturities, forwards, strikes, marketQuotes)`: Calibrate model parameters to market data.

### `LocalVol` Class
*Use cases: [`tests/models/test_local_vol.py`](../tests/models/test_local_vol.py)*
A class for Local Volatility modeling calibrated to SABR slices.
*   **Constructor**: `LocalVol(sabr_models: List[Sabr])`
*   **Properties**: `spot` (getter/setter).
*   **Methods**:
    *   `callPDE(maturity, strike, N=100)`: Price a call option using a Finite Difference solver.
    *   `putPDE(maturity, strike, N=100)`: Price a put option using a Finite Difference solver.
    *   `density(maturity, Nt)`: Returns the transition density for a given maturity.

### SABR PDE Solvers (`AfSabr`, `AntonovSabr`)
*Use cases: [`tests/pde/test_sabr_pde.py`](../tests/pde/test_sabr_pde.py)*
Classes for solving the SABR PDE directly.
*   **Constructors**:
    *   `AfSabr(alpha, beta, nu, rho, shift, maturity, F, sizeX, sizeT, nd)`
    *   `AntonovSabr(alpha, beta, nu, rho, shift, maturity, F, sizeX, sizeT, nd)`
*   **Methods**:
    *   `getDensity()`: Returns the probability density vector.
    *   `getFgrid()`: Returns the forward grid coordinates.

---

## Interest Rate & Hybrid Models

### `HullWhite` Class
*Use cases: [`tests/models/test_hull_white.py`](../tests/models/test_hull_white.py)*
A class representing the Hull-White short-rate model.
*   **Constructor**: `HullWhite(kappa, timeSigmas, sigmas, timeDFs, DFs)`
*   **Methods**:
    *   `optionBond(expiry, maturity, strike, callORput)`: Price a European option on a Zero-Coupon Bond.
    *   `swaption(expiry, tenor, strike, payFrequency=0.5)`: Price a European swaption.
    *   `simulationHW(times)`: Returns a simulated short-rate path.

### `HHW` Class (Hybrid Hull-White)
*Use cases: [`tests/models/test_hull_white.py`](../tests/models/test_hull_white.py)*
A hybrid model combining Hull-White interest rates with stochastic equity volatility.
*   **Constructor**: `HHW(s0, v0, r0, kappa, eta, rho, sigma1, sigma2, a)`
*   **Methods**:
    *   `HHWPrice(maturity, strike)`: Price a European option under the hybrid model.

### `Termstructure` Class
*Use cases: [`tests/termstructure/test_termstructure.py`](../tests/termstructure/test_termstructure.py)*
A wrapper around QuantLib term structures.
*   **Constructor**: `Termstructure(days, rates, calendar, daycount)`
*   **Methods**:
    *   `discount(date_serial)`: Returns the discount factor for a given date serial.
    *   `rate(date_serial, tenor_days)`: Returns the forward rate.

---

## Tree Models

### `CTree` Class
*Use cases: [`tests/models/test_tree.py`](../tests/models/test_tree.py)*
Binomial and Trinomial tree implementations for option pricing.
*   **Constructor**: `CTree(S, T, F, IV)` where `S` is spot, `T` is times vector, `F` is forwards, `IV` is implied vols.
*   **Enums**:
    *   `OptionType` (`Call`, `Put`)
    *   `TreeType` (`Recombining`, `NonRecombining`)
    *   `ExerciseStyle` (`American`, `European`, `Bermudan`)
*   **Methods**:
    *   `calculateBinomial(strike, Maturity, Nnodes, style, pay, tree)`: Binomial tree pricing.
    *   `calculateTrinomial(strike, Maturity, Nnodes, style, pay, tree)`: Trinomial tree pricing.

---

## Monte Carlo

### `SkewMC` Class
*Use cases: [`tests/models/test_skew_mc.py`](../tests/models/test_skew_mc.py)*
Monte Carlo engine using skewed volatility surfaces via SABR models.
*   **Constructor**: `SkewMC(sabrModels: List[Sabr])`
*   **Methods**:
    *   `simulation(times, spot, kappa)`: Returns simulated spot values at given times.

---

## Fixed Income PDE Solvers

### `HWPDE` Class
*Use cases: [`tests/pde/test_fi_pde.py`](../tests/pde/test_fi_pde.py)*
Hull-White PDE solver for European and Bermudan swaptions, callable swaps.
*   **Constructors**:
    *   `HWPDE(R0, kappa, timeSigmas, sigmas, timeThetas, thetas)`
    *   `HWPDE(kappa, timeSigmas, sigmas, timeDF, DF)` (calibration-based)
*   **Methods**:
    *   `pricingSwaption(Expiry, Tenor, Strike, PayFrequency)`: European swaption.
    *   `pricingBermudan(Expiry, Tenor, Exercises, Strike, PayFrequency)`: Bermudan swaption.
    *   `pricingCallableSwap(Expiry, Tenor, Exercises, Coupon, Strike, PayFrequency, type)`: Callable swap.

### `ShortRate1FPDE` Class
*Use cases: [`tests/pde/test_fi_pde.py`](../tests/pde/test_fi_pde.py)*
Generic 1-factor short-rate model PDE solver.
*   **Constructor**: `ShortRate1FPDE(R0, kappa, alpha, beta, gamma, timeSigmas, sigmas, timeThetas, thetas)`
*   **Methods**:
    *   `pricingSwaption(Expiry, Tenor, Strike, PayFrequency=0.5)`: European swaption.

### `ShortRate2FPDE` Class
*Use cases: [`tests/pde/test_fi_pde.py`](../tests/pde/test_fi_pde.py)*
2-factor short-rate model (G2++) PDE solver using ADI methods.
*   **Constructor**: `ShortRate2FPDE(kappa1, kappa2, lambda, timeSigma1s, sigma1s, timeSigma2s, sigma2s, timeAlphas, alphas)`
*   **Methods**:
    *   `pricingSwaption(Expiry, Tenor, Strike, PayFrequency=0.5)`: European swaption.
