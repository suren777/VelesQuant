# Hull-White 1-Factor Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Hull-White 1-Factor Short Rate Model |
| **Model ID** | `VQ-IR-001` |
| **Model Tier / Materiality** | `Tier 1 (Critical)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price and hedge interest rate derivatives including swaptions, bond options, callable/cancellable swaps, and Bermudan swaptions.
* **Methodology:** Mean-reverting short-rate model with time-dependent volatility, solved via both analytical (Jamshidian decomposition) and PDE methods.
* **Performance:** Calibrates to ATM swaption volatilities with typical RMSE < 1bp.
* **Primary Limitation:** Single-factor structure cannot capture decorrelation between different tenor rates.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Valuation and risk management of vanilla and exotic interest rate derivatives.
* **Intended Use:** 
  - Swaption pricing and hedging
  - Bond option valuation
  - Callable/Bermudan derivative pricing
  - Portfolio-level NPV computation
* **User Group:** Quantitative analysts, traders, risk managers.

### 2.2 Scope and Coverage

* **Target Population:** Interest rate derivatives across major currencies (USD, EUR, GBP).
* **Exclusions:**
  - Multi-factor dynamics (use 2-factor models for spread options)
  - Negative rate modeling without shift extension
  - Credit-contingent derivatives

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The Hull-White model specifies short rate dynamics under the risk-neutral measure:

$$
dr_t = (\theta(t) - \kappa \cdot r_t) \, dt + \sigma(t) \, dW_t
$$

Where:
- $r_t$: Instantaneous short rate
- $\kappa$: Mean reversion speed (constant)
- $\theta(t)$: Time-dependent drift calibrated to fit the initial term structure
- $\sigma(t)$: Piecewise-constant volatility term structure
- $W_t$: Standard Brownian motion

**Key Properties:**
- Gaussian distribution of rates (can go negative)
- Affine term structure: $P(t,T) = A(t,T) \cdot e^{-B(t,T) \cdot r_t}$
- Analytical bond option pricing via Jamshidian decomposition

**Justification:** Chosen for:
1. Analytical tractability for European products
2. Efficient calibration to swaption market
3. Extensibility to Bermudans via PDE solver

### 3.2 Alternatives Considered

| Alternative | Reason for Rejection |
|-------------|---------------------|
| Black-Karasinski | Log-normal rates add complexity without proportional benefit for EUR/GBP markets |
| G2++ (2-Factor) | Higher computational cost; not needed for single-tenor hedging |
| SABR-LMM | Overkill for vanilla swaptions; use for exotics requiring skew |

### 3.3 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Gaussian short rate | Analytical tractability | **Medium:** Mispricing at low/negative rates | Monitor rate distribution | **Amber** |
| A.2 | Constant mean reversion | Calibration stability | **Low:** Minor hedging error | Annual recalibration | **Green** |
| A.3 | Perfect correlation across tenors | Single-factor limitation | **High:** CMS spread mispricing | Benchmark vs 2-factor | **Red** |
| A.4 | Continuous trading | Standard derivative pricing | **Low:** Negligible for liquid markets | — | **Green** |

---

## 4. Data

### 4.1 Data Sources

* **Discount Curve:** Zero-coupon bond prices or discount factors at specified tenors
* **Swaption Volatilities:** ATM normal or lognormal vols for calibration
* **Source Systems:** Market data feeds, Bloomberg, Reuters

### 4.2 Processing & Cleaning

* **Interpolation:** Log-linear interpolation of discount factors
* **Missing Data:** Require complete curve; no imputation

### 4.3 Data Quality & Proxies

* **Quality Checks:** Arbitrage-free curve construction verified
* **Proxies:** None required for standard calibration

---

## 5. Model Development

### 5.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Python Wrapper                           │
│  hullwhite.py: HullWhiteModel                               │
│  ├── price(instrument, market_data)                        │
│  ├── simulate(times, curve)                                │
│  └── price_bond_option(expiry, maturity, strike, curve)    │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                    C++ Native Module                        │
├─────────────────────────────────────────────────────────────┤
│  hw.cpp: Analytical Methods                                 │
│  ├── swaption() - Jamshidian decomposition                 │
│  ├── optionBond() - Black formula on bonds                 │
│  ├── calibrator() - Levenberg-Marquardt optimization       │
│  └── calibratorBstrp() - Bootstrapping approach            │
├─────────────────────────────────────────────────────────────┤
│  hw_pde.cpp: PDE Methods (HWPDE class)                     │
│  ├── pricingSwaption() - European swaptions                │
│  ├── pricingBermudan() - Bermudan exercise                 │
│  ├── pricingCallableSwap() - Callable structures           │
│  └── termStructureCalibrator() - θ(t) fitting              │
└─────────────────────────────────────────────────────────────┘
```

### 5.2 Parameter Estimation

**Parameters:**
- $\kappa$: Mean reversion speed, calibrated from swaption tenor structure
- $\sigma(t)$: Piecewise-constant volatilities at swaption expiries

**Calibration Methods:**

1. **Levenberg-Marquardt (`calibrator`):**
   - Joint optimization of all $\sigma_i$ and $\kappa$
   - Objective: Minimize squared error vs market swaption prices/vols

2. **Bootstrap (`calibratorBstrp`):**
   - Sequential fitting: first $\kappa$ from vol ratios, then $\sigma_i$ per expiry
   - Faster but less robust for scarce data

**Calibration Target Options:**
- `Price`: Match swaption prices directly
- `Volatility`: Match implied volatilities

---

## 6. Testing and Performance

### 6.1 Statistical Performance

* **In-Sample:** Typical calibration error < 0.5bp on ATM vols
* **Out-of-Sample:** Validated against QuantLib reference implementation
* **Stability:** PSI not applicable (deterministic model)

### 6.2 Supported Products

| Product | Method | Validated |
|---------|--------|-----------|
| European Swaption | Analytical | ✅ |
| Zero-Coupon Bond Option | Analytical | ✅ |
| Coupon Bond Option | PDE | ✅ |
| Bermudan Swaption | PDE | ✅ |
| Callable Swap | PDE | ✅ |

### 6.3 Benchmarking

* Compared against QuantLib `HullWhite` implementation
* Monte Carlo cross-validation for complex products

---

## 7. Limitations and Post-Model Adjustments

### 7.1 Known Limitations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Negative rates possible | Underpricing floors in neg-rate environment | Use shifted model variant |
| Single factor | Cannot price CMS spreads accurately | Use 2-factor or copula overlay |
| Constant κ | May misfit long-dated vols | Periodic recalibration |

### 7.2 Post-Model Adjustments (Overlays)

* **Shift Parameter:** Shift curve by fixed amount for negative rate environments (not currently implemented in Python wrapper)

---

## 8. Ongoing Monitoring and Governance

### 8.1 Governance Framework

> This model operates under the VelesQuant Model Risk Management Framework.
> * **1st Line (Developer):** Monitor calibration quality daily
> * **2nd Line (Validation):** Independent quarterly review
> * **3rd Line (Audit):** Annual compliance verification

### 8.2 Performance Triggers

| Metric | Green | Amber | Red |
| --- | --- | --- | --- |
| **Calibration Error** | < 1bp | 1-5bp | > 5bp |
| **Pricing vs Market** | < 2bp | 2-10bp | > 10bp |
| **Recalibration Frequency** | Daily | Weekly | Monthly |

---

## 9. Appendices

### A. API Reference

```python
from velesquant.models import HullWhiteModel
from velesquant.market.curves import DiscountCurve

# Initialize model
model = HullWhiteModel(kappa=0.03, sigma=0.01)

# Price a swaption
from velesquant.instruments.rates import Swaption
swaption = Swaption(expiry=1.0, tenor=5.0, strike=0.02, pay_frequency=0.5)
price = model.price(swaption, curve)

# Price a bond option
price = model.price_bond_option(
    expiry=1.0, 
    maturity=5.0, 
    strike=0.95, 
    curve=curve,
    option_type="Call"
)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [hullwhite.py](../../src/velesquant/models/hullwhite.py) |
| C++ Analytical | [hw.cpp](../../src/models/hw.cpp) |
| C++ PDE Solver | [hw_pde.cpp](../../src/pde_solvers/hw_pde.cpp) |

### C. References

1. Hull, J. & White, A. (1990). "Pricing Interest-Rate-Derivative Securities"
2. Brigo, D. & Mercurio, F. (2006). "Interest Rate Models - Theory and Practice"
3. PRA SS1/23: Model Risk Management Principles
4. Federal Reserve SR 11-7: Guidance on Model Risk Management

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
