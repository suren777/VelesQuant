# VelesQuant

**VelesQuant** is a high-performance financial valuation library written in C++ and exposed to Python using **pybind11**. It provides tools for pricing complex financial derivatives using advanced models and numerical techniques.

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
| SABR | [`SabrModel`](docs/models/sabr.md) | Stochastic Alpha Beta Rho volatility model |
| Heston | [`HestonModel`](docs/models/heston.md) | Stochastic variance model |
| Local Vol | [`LocalVolModel`](docs/models/localvol.md) | Local volatility from SABR slices |
| Hull-White | [`HullWhiteModel`](docs/models/hullwhite.md) | 1-Factor short rate model |
| ShortRate2F | `ShortRate2FPDEModel` | 2-Factor G2++ short rate model |
| CMS | [`CMSModel`](docs/models/cms.md) | Constant Maturity Swap pricing |
| CMS Spread | [`CMSSpreadModel`](docs/models/cms_spread.md) | CMS Spread Option pricing (Copula) |
| Basket | [`LogNormalBasketModel`](docs/models/basket.md) | Multi-asset Log-Normal Basket |
| Quantoed CMS | [`QuantoedCMSModel`](docs/models/quantoed.md) | Quantoed Constant Maturity Swap |
| Quantoed Spread | [`QuantoedCMSSpreadModel`](docs/models/quantoed.md) | Quantoed CMS Spread Option |
| Heston Hull-White | [`HybridHWModel`](docs/models/hybrid_hw.md) | Hybrid Equity-Interest Rate Model |


## Supported Instruments

For a comprehensive list of supported financial instruments and their pricing methods, please refer to [Supported Instruments](docs/supported_instruments.md).

## Testing


### Python Tests
```bash
poetry run pytest
```

### C++ Tests
```bash
./build/tests_cpp/velesquant_tests
```

## Contributing

We welcome contributions! Please follow these steps:

1.  **Bug Reports**: Open an issue describing the bug and providing a reproducible example.
2.  **Feature Requests**: Use the issue tracker to propose new features.
3.  **Pull Requests**: 
    - Fork the repository.
    - Create a feature branch.
    - Ensure all C++ and Python tests pass.
    - Submit a PR with a clear description of changes.

## License
[MIT License](LICENSE)
