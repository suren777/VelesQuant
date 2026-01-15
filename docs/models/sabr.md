# SABR Stochastic Volatility Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | SABR Stochastic Alpha Beta Rho Model |
| **Model ID** | `VQ-VOL-001` |
| **Model Tier / Materiality** | `Tier 1 (Critical)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Model volatility smile dynamics for pricing and hedging interest rate and equity options.
* **Methodology:** Stochastic volatility model with CEV-type forward dynamics and correlated volatility process.
* **Performance:** Calibrates to full smile with typical error < 0.5 vol points.
* **Primary Limitation:** Requires careful handling of low strikes ( β < 1) to avoid numerical instability.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Capture volatility smile/skew observed in option markets.
* **Intended Use:**
  - Swaption smile modeling
  - Caplet/floorlet pricing
  - Local volatility surface construction
  - Monte Carlo simulation of forward rates
* **User Group:** Quantitative analysts, derivatives traders.

### 2.2 Scope and Coverage

* **Target Population:** European-style options on rates or equities.
* **Exclusions:**
  - American/Bermudan options (use PDE extensions)
  - Path-dependent exotics (use simulation mode)

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The SABR model specifies joint dynamics under the forward measure:

$$
dF_t = \alpha_t F_t^\beta dW_t^F
$$

$$
d\alpha_t = \nu \alpha_t dW_t^\alpha
$$

$$
dW_t^F \cdot dW_t^\alpha = \rho \, dt
$$

Where:
- $F_t$: Forward rate/price
- $\alpha_t$: Stochastic volatility (initial value $\alpha_0 = \alpha$)
- $\beta$: CEV exponent (typically fixed, 0 ≤ β ≤ 1)
- $\nu$: Volatility of volatility
- $\rho$: Correlation between forward and volatility (typically negative)

**Approximate Implied Volatility Formula (Hagan et al.):**

For strike $K$ and forward $F$:

$$
\sigma_{BS}(K) = \frac{\alpha}{(FK)^{(1-\beta)/2}} \cdot \frac{z}{x(z)} \cdot \left(1 + \epsilon T\right)
$$

where $z = \frac{\nu}{\alpha}(FK)^{(1-\beta)/2}\log(F/K)$ and correction terms depend on $\beta$, $\rho$, $\nu$.

**Justification:** Industry-standard model for smile dynamics, analytically tractable, and widely accepted by regulators.

### 3.2 Alternatives Considered

| Alternative | Reason for Rejection |
|-------------|---------------------|
| Local Volatility | Cannot reproduce smile dynamics over time |
| Heston | More parameters, less transparent for IR applications |
| SVI | Purely parametric, no dynamic interpretation |

### 3.3 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Log-normal volatility | Standard SABR formulation | **Low:** Minor bias in extreme scenarios | N/A | **Green** |
| A.2 | Constant β | Fixed during calibration | **Medium:** Smile shape constraint | Validate β choice | **Amber** |
| A.3 | No jumps | Diffusion-only dynamics | **Medium:** Fat tails underestimated | Compare to market | **Amber** |
| A.4 | ATM approximation accuracy | Hagan formula is asymptotic | **Low:** Small error near ATM | Monitor OTM fit | **Green** |

---

## 4. Data

### 4.1 Data Sources

* **Forward Rate:** From yield curve or forward swap rate
* **Option Quotes:** Market volatilities or prices at multiple strikes
* **Source Systems:** Bloomberg, Reuters, internal pricing systems

### 4.2 Processing & Cleaning

* **Strike Grid:** Typically 5-11 strikes around ATM
* **Outlier Treatment:** Exclude quotes with bid-ask > threshold

### 4.3 Data Quality & Proxies

* **Quality Checks:** Arbitrage-free smile verification
* **Proxies:** None required for calibration

---

## 5. Model Development

### 5.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Python Wrapper                           │
│  sabr.py: SabrModel                                        │
│  ├── implied_vol(strike)                                   │
│  └── calibrate(strikes, quotes, target)                    │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                    C++ Native Module                        │
│  sabr.cpp: Sabr class                                      │
│  ├── impliedVol() - Hagan approximation                    │
│  ├── normalVol() - Bachelier volatility                    │
│  ├── premiumBlackScholes() / premiumBachelier()            │
│  ├── localVol() - Dupire extraction                        │
│  ├── calibrator() - Levenberg-Marquardt                    │
│  └── simulation() - Quantile-based Monte Carlo             │
└─────────────────────────────────────────────────────────────┘
```

### 5.2 Parameter Estimation

**Parameters to Calibrate:**
- $\alpha$: Initial volatility level
- $\nu$: Vol-of-vol
- $\rho$: Correlation (constrained to [-1, 1])

**Fixed Parameter:**
- $\beta$: User-specified (typically 0.5 for rates, 1.0 for equities)

**Calibration Methods:**

1. **Full Calibration (`calibrator`):**
   - Optimizes α, ν, ρ jointly via Levenberg-Marquardt
   - Objective: Minimize squared error vs market vols/prices

2. **ATM-Constrained (`calibratorWithInitialATM`):**
   - Fits α analytically to match ATM vol
   - Optimizes ν, ρ only

**Calibration Target Options:**
- `Volatility`: Match implied volatilities
- `Price`: Match option premiums

---

## 6. Testing and Performance

### 6.1 Statistical Performance

* **In-Sample:** Typical fit error < 0.3 vol points
* **Out-of-Sample:** Validated against QuantLib SABR implementation
* **Stability:** Robust for β ∈ [0, 1], ρ ∈ (-1, 1)

### 6.2 Supported Features

| Feature | Method | Validated |
|---------|--------|-----------|
| Black-Scholes Implied Vol | Hagan Formula | ✅ |
| Normal (Bachelier) Vol | Separate formula | ✅ |
| European Option Pricing | BS / Bachelier | ✅ |
| Local Volatility Extraction | Dupire | ✅ |
| Monte Carlo Simulation | Quantile Table | ✅ |

### 6.3 Benchmarking

* Compared against QuantLib `SabrSmileSection`
* Validated local vol extraction against Dupire formula

---

## 7. Limitations and Post-Model Adjustments

### 7.1 Known Limitations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Low-strike instability | Negative density possible for β < 1 | Use shifted SABR |
| Asymptotic approximation | Less accurate far OTM | Use full integration |
| Static smile | No term structure dynamics | Use multiple expiry slices |

### 7.2 Post-Model Adjustments

* **Shift Parameter:** Add constant shift to handle negative rates
* **Backbone Interpolation:** Smooth parameters across expiries

---

## 8. Ongoing Monitoring and Governance

### 8.1 Governance Framework

> This model operates under the VelesQuant Model Risk Management Framework.
> * **1st Line (Developer):** Daily calibration quality monitoring
> * **2nd Line (Validation):** Quarterly independent review
> * **3rd Line (Audit):** Annual compliance verification

### 8.2 Performance Triggers

| Metric | Green | Amber | Red |
| --- | --- | --- | --- |
| **Smile Fit Error** | < 0.5 vol | 0.5-2 vol | > 2 vol |
| **ρ Bounds** | \|ρ\| < 0.95 | \|ρ\| > 0.95 | ρ at boundary |
| **Arbitrage** | None | Wing issues | Butterfly violations |

---

## 9. Appendices

### A. API Reference

```python
from velesquant.models import SabrModel

# Initialize model
model = SabrModel(
    maturity=1.0,
    forward=0.03,
    beta=0.5,
    alpha=0.2,
    nu=0.5,
    rho=-0.3
)

# Get implied volatility
vol = model.implied_vol(strike=0.025)

# Calibrate to market
model = model.calibrate(
    strikes=[0.02, 0.025, 0.03, 0.035, 0.04],
    quotes=[0.25, 0.22, 0.20, 0.22, 0.26],
    calibration_target="Volatility"
)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [sabr.py](../../src/velesquant/models/sabr.py) |
| C++ Implementation | [sabr.cpp](../../src/volatility/sabr.cpp) |

### C. References

1. Hagan, P. et al. (2002). "Managing Smile Risk"
2. Oblój, J. (2008). "Fine-tune your smile: Correction to Hagan et al."
3. PRA SS1/23: Model Risk Management Principles
4. Federal Reserve SR 11-7: Guidance on Model Risk Management

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
