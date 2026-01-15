# Heston Stochastic Volatility Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Heston Stochastic Volatility Model |
| **Model ID** | `VQ-VOL-002` |
| **Model Tier / Materiality** | `Tier 1 (Critical)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price and hedge equity and FX options with stochastic volatility, capturing smile dynamics.
* **Methodology:** Mean-reverting variance process with semi-analytical pricing via characteristic functions.
* **Performance:** Calibrates to full volatility surface; accurate pricing across strikes and maturities.
* **Primary Limitation:** CIR variance process can become negative in discrete simulation (Milstein/Euler).

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Capture equity volatility smile and term structure for accurate derivative pricing.
* **Intended Use:**
  - European option pricing (calls/puts)
  - Exotic pricing via Monte Carlo (barriers, cliquets)
  - Volatility surface calibration
* **User Group:** Equity derivatives traders, quantitative analysts.

### 2.2 Scope and Coverage

* **Target Population:** Equity and FX options.
* **Exclusions:**
  - Interest rate derivatives (use Hull-White)
  - Credit-contingent products

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The Heston model specifies joint dynamics:

$$
dS_t = S_t \sqrt{V_t} \, dW_t^S
$$

$$
dV_t = \kappa (\theta - V_t) \, dt + \xi \sqrt{V_t} \, dW_t^V
$$

$$
dW_t^S \cdot dW_t^V = \rho \, dt
$$

Where:
- $S_t$: Spot price
- $V_t$: Instantaneous variance (with initial value $V_0$)
- $\kappa$: Mean reversion speed
- $\theta$: Long-run variance
- $\xi$: Volatility of volatility
- $\rho$: Spot-vol correlation (typically negative for equities)

**Feller Condition:** $2\kappa\theta > \xi^2$ ensures variance stays positive.

**Characteristic Function:**
The model admits a closed-form characteristic function, enabling semi-analytical European option pricing via Fourier inversion (Gauss-Laguerre quadrature).

**Justification:** Industry-standard for equity derivatives; balances analytical tractability with realistic smile dynamics.

### 3.2 Alternatives Considered

| Alternative | Reason for Rejection |
|-------------|---------------------|
| SABR | Less natural for equity (no mean-reverting vol) |
| Local Volatility | Static; doesn't capture forward smile dynamics |
| Bates (Heston + Jumps) | Added complexity; use when jumps observed |

### 3.3 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Mean-reverting variance | Empirical observation | **Low:** Core model feature | Calibration stability | **Green** |
| A.2 | Feller condition | Ensures positive variance | **High:** Negative variance in simulation | Check $2\kappa\theta > \xi^2$ | **Amber** |
| A.3 | Continuous paths | No jumps in spot/vol | **Medium:** Fat tails underestimated | Compare to market | **Amber** |
| A.4 | Constant parameters | No term structure in params | **Low:** Use separate calibrations | Per-maturity fit | **Green** |

---

## 4. Data

### 4.1 Data Sources

* **Spot Price:** Real-time market feed
* **Forward Curve:** Dividend-adjusted forwards
* **Option Quotes:** Implied vols or prices across strike/maturity grid
* **Source Systems:** Bloomberg, Reuters, exchange feeds

### 4.2 Processing & Cleaning

* **Strike Grid:** Typically 25-delta, 10-delta, ATM strikes per maturity
* **Quote Type:** Can calibrate to prices or implied vols

### 4.3 Data Quality & Proxies

* **Quality Checks:** Arbitrage-free surface verification
* **Proxies:** None required

---

## 5. Model Development

### 5.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Python Wrapper                           │
│  heston.py: HestonModel                                    │
│  ├── price_option(maturity, forward, strike, type)         │
│  ├── simulate(times, forwards)                             │
│  └── calibrate(maturities, forwards, strikes, quotes)      │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                    C++ Native Module                        │
│  s_vol.cpp: sVol class                                     │
│  ├── hestonPriceCF() - Characteristic function pricing     │
│  ├── hestonPrice() - Alternative integration               │
│  ├── simulationHeston() - Milstein scheme MC               │
│  ├── simulationHestonDNT() - Double-no-touch barriers      │
│  ├── simulationHestonCliq() - Cliquet simulation           │
│  └── calibrator() - Levenberg-Marquardt optimization       │
└─────────────────────────────────────────────────────────────┘
```

### 5.2 Parameter Estimation

**Parameters:**
- $V_0$: Initial variance
- $\kappa$: Mean reversion speed
- $\theta$: Long-run variance
- $\xi$: Vol-of-vol
- $\rho$: Spot-vol correlation

**Calibration Method:**
- Levenberg-Marquardt optimization
- Objective: Minimize squared pricing error vs market
- Weights: ATM-focused (exponential decay from spot)
- Feller condition enforced via penalty function

**Calibration Variants:**
- `calibrator`: Standard 5-parameter calibration
- `IVcalibrator`: Match implied volatilities
- `FXcalibrator`: 4-parameter variant with fixed κ

---

## 6. Testing and Performance

### 6.1 Statistical Performance

* **In-Sample:** Typical fit error < 1% of option price
* **Out-of-Sample:** Validated against QuantLib `HestonProcess`
* **Simulation:** Milstein scheme with variance floor for stability

### 6.2 Supported Products

| Product | Method | Validated |
|---------|--------|-----------|
| European Call/Put | Characteristic Function | ✅ |
| Down-and-Out Barrier | Monte Carlo | ✅ |
| Cliquet / Accumulator | Monte Carlo | ✅ |
| Maximum Lookback | Monte Carlo | ✅ |

### 6.3 Benchmarking

* Compared against QuantLib `AnalyticHestonEngine`
* Monte Carlo convergence verified

---

## 7. Limitations and Post-Model Adjustments

### 7.1 Known Limitations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Feller violation | Negative variance in simulation | Absorbing boundary / reflection |
| No jumps | Underestimates crash risk | Use Bates extension |
| Single maturity calibration | No term structure | Separate calibration per slice |

### 7.2 Post-Model Adjustments

* **Variance Floor:** Apply `abs(V)` in simulation to prevent negative values
* **Weighted Calibration:** ATM-focused weights for stability

---

## 8. Ongoing Monitoring and Governance

### 8.1 Governance Framework

> This model operates under the VelesQuant Model Risk Management Framework.
> * **1st Line (Developer):** Daily calibration monitoring
> * **2nd Line (Validation):** Quarterly independent review
> * **3rd Line (Audit):** Annual compliance verification

### 8.2 Performance Triggers

| Metric | Green | Amber | Red |
| --- | --- | --- | --- |
| **Calibration Error** | < 2% | 2-5% | > 5% |
| **Feller Ratio** | > 1.1 | 0.8-1.1 | < 0.8 |
| **ρ Stability** | Stable | Jumpy | Unrealistic |

---

## 9. Appendices

### A. API Reference

```python
from velesquant.models import HestonModel

# Initialize model
model = HestonModel(
    spot=100.0,
    var0=0.04,
    kappa=2.0,
    theta=0.04,
    xi=0.5,
    rho=-0.7
)

# Price a European call
price = model.price_option(
    maturity=1.0,
    forward=100.0,
    strike=100.0,
    option_type="call"
)

# Simulate paths
paths = model.simulate(
    times=[0.25, 0.5, 0.75, 1.0],
    forwards=[100.5, 101.0, 101.5, 102.0]
)

# Calibrate to market
model = model.calibrate(
    maturities=[0.25, 0.5, 1.0],
    forwards=[100.0, 100.0, 100.0],
    strikes=[90.0, 100.0, 110.0],
    quotes=[12.5, 8.0, 5.5],
    calibration_target="Price"
)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [heston.py](../../src/velesquant/models/heston.py) |
| C++ Implementation | [s_vol.cpp](../../src/volatility/s_vol.cpp) |
| C++ Facade | [heston_facade.cpp](../../src/models/heston_facade.cpp) |

### C. References

1. Heston, S. (1993). "A Closed-Form Solution for Options with Stochastic Volatility"
2. Gatheral, J. (2006). "The Volatility Surface"
3. PRA SS1/23: Model Risk Management Principles
4. Federal Reserve SR 11-7: Guidance on Model Risk Management

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
