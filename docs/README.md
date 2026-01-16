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

VelesQuant provides high-level Python wrappers in `velesquant.models` for ease of use.

```python
import numpy as np
from velesquant.models import SabrModel, HestonModel, HullWhiteModel
from velesquant.market.curves import DiscountCurve

# --- SABR Volatility Model ---
sabr = SabrModel(maturity=1.0, forward=100.0, beta=0.5)

# Calibrate to market data
strikes = [90.0, 100.0, 110.0]
quotes = [0.25, 0.20, 0.22]  # Implied volatilities
sabr.calibrate(strikes=strikes, quotes=quotes, calibration_target="Volatility")

print(f"Calibrated SABR: alpha={sabr.alpha:.4f}, nu={sabr.nu:.4f}, rho={sabr.rho:.4f}")

# Get implied volatility for a strike
iv = sabr.implied_vol(105.0)

# --- Heston Stochastic Volatility ---
heston = HestonModel(spot=100.0, var0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7, seed=42)

# Price a vanilla option
price = heston.price_option(maturity=1.0, forward=100.0, strike=100.0, option_type="call")

# --- Hull-White Interest Rate Model ---
# 1. Define a Discount Curve
curve = DiscountCurve(
    times=[0.0, 1.0, 5.0, 10.0],
    dfs=[1.0, 0.95, 0.80, 0.65]
)

# 2. Initialize Model
hw = HullWhiteModel(kappa=0.1, sigma=0.01)

# 3. Price a Swaption
from velesquant.instruments.rates import Swaption
swaption = Swaption(expiry=1.0, tenor=5.0, strike=0.03, pay_frequency=0.5)
price = hw.price(swaption, curve)

# 4. Simulate Short Rate (Vectorized)
times = np.array([0.25, 0.5, 0.75, 1.0])
paths = hw.simulate(times, curve)
```

## Available Models

| Model | Class | Description |
|-------|-------|-------------|
| SABR | [`SabrModel`](models/sabr.md) | Stochastic Alpha Beta Rho volatility model |
| Heston | [`HestonModel`](models/heston.md) | Stochastic variance model |
| Local Vol | [`LocalVolModel`](models/localvol.md) | Local volatility from SABR slices |
| Hull-White | [`HullWhiteModel`](models/hullwhite.md) | 1-Factor short rate model |
| ShortRate2F | `ShortRate2FPDEModel` | 2-Factor G2++ short rate model |
| CMS | [`CMSModel`](models/cms.md) | Constant Maturity Swap pricing |
| CMS Spread | [`CMSSpreadModel`](models/cms_spread.md) | CMS Spread Option pricing (Copula) |
| Basket | [`LogNormalBasketModel`](models/basket.md) | Multi-asset Log-Normal Basket |
| Quantoed CMS | [`QuantoedCMSModel`](models/quantoed.md) | Quantoed Constant Maturity Swap |
| Quantoed Spread | [`QuantoedCMSSpreadModel`](models/quantoed.md) | Quantoed CMS Spread Option |
| Heston Hull-White | [`HybridHWModel`](models/hybrid_hw.md) | Hybrid Equity-Interest Rate Model |

## Calibration Example

All models support calibration via numpy arrays or lists:

```python
# Heston Calibration
maturities = [1.0, 1.0, 2.0]
strikes = [100.0, 110.0, 100.0]
forwards = [100.0, 100.0, 100.0]
quotes = [5.0, 2.0, 8.0] # Prices

heston.calibrate(
    maturities=maturities, 
    forwards=forwards, 
    strikes=strikes, 
    quotes=quotes, 
    calibration_target="Price"
)
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
