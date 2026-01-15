# VelesQuant Python API

VelesQuant is a high-performance library for quantitative finance, providing models for volatility, interest rates, and hybrid derivatives.

## Installation

```bash
# Using poetry (recommended)
poetry install

# Or using pip
pip install .
```

## Quick Start

```python
import numpy as np
import velesquant.native as n

# --- SABR Volatility Model ---
sabr = n.Sabr(maturity=1.0, forward=100.0, beta=0.5)

# Calibrate to market data
strikes = np.array([[90.0, 100.0, 110.0]])
quotes = np.array([[0.25, 0.20, 0.22]])  # Implied volatilities
result = sabr.calibrate(strikes=strikes, quotes=quotes, quote_type=n.CalibrationTarget.Volatility)

print(f"Calibrated SABR: alpha={sabr.alpha:.4f}, nu={sabr.nu:.4f}, rho={sabr.rho:.4f}")

# Get implied volatility for a strike
iv = sabr.impliedVol(105.0)

# --- Heston Stochastic Volatility ---
heston = n.Heston(spot=100.0, var0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7, seed=42)

# Price a vanilla option
price = heston.hestonPrice(maturity=1.0, forward=100.0, strike=100.0, optType="call")

# --- Hull-White Interest Rate Model ---
hw = n.HullWhite(
    kappa=0.1,
    timeSigmas=[1.0, 5.0],
    sigmas=[0.01, 0.012],
    timeDFs=[0.0, 1.0, 5.0, 10.0],
    DFs=[1.0, 0.95, 0.80, 0.65]
)

# Price a swaption
swaption_price = hw.swaption(Expiry=1.0, Tenor=5.0, Strike=0.03)

# Simulate short rate paths
paths = hw.simulation(times=[0.25, 0.5, 0.75, 1.0])
```

## Available Models

| Model | Class | Description |
|-------|-------|-------------|
| SABR | `Sabr` | Stochastic Alpha Beta Rho volatility model |
| Heston | `Heston` | Stochastic variance model |
| Local Vol | `LocalVol` | Local volatility from SABR slices |
| Hull-White | `HullWhite` | 1-Factor short rate model |
| ShortRate2F | `ShortRate2FPDE` | 2-Factor G2++ short rate model |
| CMS | `CMS` | Constant Maturity Swap pricing |
| Swaption | `Swaption` | Interest rate swaption pricing |

## Calibration API

All models support calibration via numpy arrays:

```python
# DefSwap structure for interest rate calibration
swap = n.DefSwap()
swap.Expiry = 1.0
swap.Tenor = 5.0
swap.SwapRate = 0.03
swap.VolATM = 0.20

swaps = [swap]
hw.calibrate(swapQuotes=swaps, target=n.CalibrationTarget.Volatility)
```

## Running Tests

```bash
poetry run pytest tests/
```

## Benchmarks

```bash
poetry run python benchmarks/bench_calibration.py
```

## License

Proprietary - VelesQuant
