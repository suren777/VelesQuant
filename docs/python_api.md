# VelesQuant Python API Reference

This document provides a detailed overview of the core functions and classes exposed to Python from the `VelesQuant` C++ library (`velesquant.native`). For higher-level Pythonic usage, consider using the wrappers in `velesquant.models`.

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
*   **Properties**: `alpha`, `nu`, `rho`, `maturity`, `forward`, `beta` (getter/setter).
*   **Methods**:
    *   `implied_vol(strike)`: Returns ATM implied volatility.
    *   `local_vol(spot)`: Returns local volatility.
    *   `normal_vol(strike)`: Returns normal volatility.
    *   `premium_black_scholes(strike, call_or_put)`: Standard BS premium using SABR vol.
    *   `premium_bachelier(strike, call_or_put)`: Bachelier premium.
    *   `simulate(times, forwards)`: Monte Carlo simulation.
    *   `calibrate(strikes, quotes, target)`: Calibrate model parameters.

### `Heston` Class
*Use cases: [`tests/models/test_heston.py`](../tests/models/test_heston.py)*
A class representing the Heston stochastic volatility model.
*   **Constructor**: `Heston(spot, var0, kappa, theta, xi, rho, seed=42)`
*   **Methods**:
    *   `price(maturity, forward, strike, opt_type)`: Pricing via characteristic function.
    *   `simulate(times, forwards)`: Returns a Monte Carlo simulation path.
    *   `calibrate(maturities, forwards, strikes, quotes, target)`: Calibrate model parameters.

### `SchobelZhu` Class
*Use cases: [`tests/models/test_schobzhu.py`](../tests/models/test_schobzhu.py)*
A class representing the Schobel-Zhu stochastic volatility model.
*   **Constructor**: `SchobelZhu(spot, var0, kappa, theta, xi, rho)`
*   **Methods**:
    *   `price(maturity, forward, strike)`: Analytical pricing.
    *   `simulate(times, forwards)`: Monte Carlo simulation.
    *   `calibrate(maturities, forwards, strikes, market_quotes)`: Calibrate model parameters.

### `LocalVol` Class
*Use cases: [`tests/models/test_local_vol.py`](../tests/models/test_local_vol.py)*
A class for Local Volatility modeling calibrated to SABR slices.
*   **Constructor**: `LocalVol(sabr_models: List[Sabr])`
*   **Properties**: `spot` (getter/setter).
*   **Methods**:
    *   `call_pde(maturity, strike, N=100)`: Price a call option using a Finite Difference solver.
    *   `put_pde(maturity, strike, N=100)`: Price a put option using a Finite Difference solver.
    *   `density(maturity, Nt)`: Returns the transition density for a given maturity.

### `LogNormalBasket` Class
*Use cases: [`tests/models/test_basket.py`](../tests/models/test_basket.py)*
A class representing a multi-asset Log-Normal basket model.
*   **Constructor**: `LogNormalBasket(Spot, Strike, Maturities, Forwards, IV, correlation)`
*   **Methods**:
    *   `simulate(schedule)`: Monte Carlo simulation of the basket assets.
    *   `simulate_with_rebalancing(schedule)`: Simulation with rebalancing logic.
    *   `get_n_assets()`: Returns number of assets.

### SABR PDE Solvers (`AfSabr`, `AntonovSabr`)
*Use cases: [`tests/pde/test_sabr_pde.py`](../tests/pde/test_sabr_pde.py)*
Classes for solving the SABR PDE directly.
*   **Constructors**:
    *   `AfSabr(alpha, beta, nu, rho, shift, maturity, F, sizeX, sizeT, nd)`
    *   `AntonovSabr(alpha, beta, nu, rho, maturity, F, sizeX, sizeT, nd)`
*   **Methods**:
    *   `getDensity()`: Returns the probability density vector.
    *   `getFgrid()`: Returns the forward grid coordinates.

### `SabrPDE` Class
*Use cases: [`tests/pde/test_sabr_pde.py`](../tests/pde/test_sabr_pde.py)*
Base or generic SABR PDE solver class.
*   **Methods**:
    *   `get_density()`: Returns density.
    *   `get_f_grid()`: Returns grid.
    *   `get_alpha()`, `get_beta()`, `get_nu()`, `get_rho()`: Getters.

---

## Interest Rate & Hybrid Models

### `HullWhiteModel` & `HullWhiteAnalyticEngine`
*Use cases: [`tests/models/test_hull_white.py`](../tests/models/test_hull_white.py)*
Classes representing the Hull-White 1-Factor short-rate model.
*   **HullWhiteModel Constructor**: `HullWhiteModel(kappa, time_sigmas, sigmas, discount_factor_times, discount_factors)`
*   **HullWhiteModel Methods**:
    *   `simulate(times)`: Returns a simulated short-rate path (list).
    *   `get_discount_factor(time)`: Get ZCB price.
    *   `get_kappa()`, `get_sigmas()`, `get_time_sigmas()`: Getters.
*   **HullWhiteAnalyticEngine Constructor**: `HullWhiteAnalyticEngine(model)`
*   **HullWhiteAnalyticEngine Methods**:
    *   `option_bond(expiry, maturity, strike, type)`: Price a European option on a Zero-Coupon Bond.
    *   `swaption(expiry, tenor, strike, pay_frequency)`: Price a European swaption.

### `HHW` Class (Heston Hull-White)
*Use cases: [`tests/models/test_hull_white.py`](../tests/models/test_hull_white.py)*
A hybrid model combining Hull-White interest rates with Heston stochastic equity volatility.
*   **Constructor**: `HHW(s0, v0, r0, kappa, eta, rho, sigma1, sigma2, a)`
*   **Methods**:
    *   `HHWPrice(maturity, strike)`: Price a European option under the hybrid model.

### `CMS` Class
*Use cases: [`tests/models/test_cms.py`](../tests/models/test_cms.py)*
Constant Maturity Swap pricing.
*   **Constructor**: `CMS(expiry_sr, tenor_sr, forward_sr, annuity_sr, pay_cms, discount_cms, beta, alpha, nu, rho)`
*   **Methods**:
    *   `fair_value(strike, type)`: Returns fair value of CMS option.
    *   `get_implied_vol(strike)`: Returns implied volatility.

### `CmsSpread` Class
*Use cases: [`tests/models/test_cms_spread.py`](../tests/models/test_cms_spread.py)*
CMS Spread option pricing using copulas.
*   **Constructor**: `CmsSpread(expiry1, tenor1, fwd1, ..., expiry2, tenor2, fwd2, ..., corr)`
*   **Methods**:
    *   `spread_option(K, a, b)`: Price a spread option (payoff $max(a S_1 - b S_2 - K, 0)$).
    *   `simulate()`: Simulate spread paths.

### `Swaption` Class
*Use cases: [`tests/models/test_swaption.py`](../tests/models/test_swaption.py)*
Native Swaption pricing class.
*   **Constructor**: `Swaption(expiry, tenor, forward, annuity, beta, alpha, nu, rho)`
*   **Methods**:
    *   `fair_value(strike, call_put)`: Price swaption using SABR.
    *   `swap_fair_value(strike)`: Price underlying swap.
    *   `get_implied_vol(strike)`: Get implied volatility.

### `Termstructure` Class
*Use cases: [`tests/termstructure/test_termstructure.py`](../tests/termstructure/test_termstructure.py)*
A wrapper around QuantLib term structures.
*   **Constructor**: `Termstructure(days, rates, calendar, daycount)`
*   **Methods**:
    *   `discount(date_serial)`: Returns the discount factor for a given date serial.
    *   `rate(date_serial, tenor_days)`: Returns the forward rate.

### `QuantoedCMS` Class
*Use cases: [`tests/models/test_quantoed_cms.py`](../tests/models/test_quantoed_cms.py)*
Quantoed Constant Maturity Swap pricing.
*   **Methods**:
    *   `fair_value(strike, type)`: Returns fair value.
    *   `get_forward()`: Returns the forward rate.
    *   `simulate(times)`: Returns simulation path.

### `QuantoedCmsSpread` Class
*Use cases: [`tests/models/test_quantoed_cms_spread.py`](../tests/models/test_quantoed_cms_spread.py)*
Quantoed CMS Spread option pricing.
*   **Methods**:
    *   `simulate()`: Simulate spread paths.

---

## Tree Models

### `CTree` Class
*Use cases: [`tests/models/test_tree.py`](../tests/models/test_tree.py)*
Binomial and Trinomial tree implementations for option pricing.
*   **Constructor**: `CTree(S, T, F, IV)`
*   **Enum Arguments**: `OptionType.{Call, Put}`, `TreeType.{Recombining, NonRecombining}`, `ExerciseStyle.{American, European, Bermudan}`.
*   **Methods**:
    *   `calculate_binomial(strike, maturity, n_nodes, style, call_put, tree_type)`: Binomial tree pricing.
    *   `calculate_trinomial(strike, maturity, n_nodes, style, call_put, tree_type)`: Trinomial tree pricing.

---

## Monte Carlo

### `SkewMC` Class
*Use cases: [`tests/models/test_skew_mc.py`](../tests/models/test_skew_mc.py)*
Monte Carlo engine using skewed volatility surfaces via SABR models.
*   **Constructor**: `SkewMC(sabrModels: List[Sabr])`
*   **Methods**:
    *   `simulate(times, spot, kappa)`: Returns simulated spot values at given times.

---

## Fixed Income PDE Solvers

### `HWPDE` Class
*Use cases: [`tests/pde/test_fi_pde.py`](../tests/pde/test_fi_pde.py)*
Hull-White PDE solver for European and Bermudan swaptions, callable swaps.
*   **Constructors**:
    *   `HWPDE(R0, kappa, timeSigmas, sigmas, timeThetas, thetas)`
    *   `HWPDE(kappa, timeSigmas, sigmas, discount_factor_times, discount_factors)`
*   **Methods**:
    *   `price_swaption(expiry, tenor, strike, pay_frequency)`: European swaption.
    *   `price_bermudan(expiry, tenor, exercises, strike, pay_frequency)`: Bermudan swaption.
    *   `price_callable_swap(...)`: Callable swap pricing.
    *   `price_zero_bond(maturity)`: Zero bond.
    *   `price_zero_bond_option(...)`: Option on zero bond.
    *   `simulate(times)`: Simulation grid.
    *   `calibrate(...)`: Calibrate to swaps.

### `ShortRate1FPDE` Class
*Use cases: [`tests/pde/test_fi_pde.py`](../tests/pde/test_fi_pde.py)*
Generic 1-factor short-rate model PDE solver.
*   **Constructor**: `ShortRate1FPDE(initial_rate, kappa, alpha, beta, gamma, time_sigmas, sigmas, time_thetas, thetas)`
*   **Methods**:
    *   `pricingBermudan(...)`: Bermudan pricing.
    *   `price_swaption(...)`: European swaption.

### `ShortRate2FPDE` Class
*Use cases: [`tests/pde/test_fi_pde.py`](../tests/pde/test_fi_pde.py)*
2-factor short-rate model (G2++) PDE solver using ADI methods.
*   **Constructor**: `ShortRate2FPDE(kappa1, kappa2, lambda, time_sigma1s, sigma1s, time_sigma2s, sigma2s, time_alphas, alphas)`
*   **Methods**:
    *   `price_swaption(expiry, tenor, strike, pay_frequency)`: European swaption.
