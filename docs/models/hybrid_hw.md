# Heston Hull-White (Hybrid) Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Heston Hull-White Hybrid Model |
| **Model ID** | `VQ-HYB-001` |
| **Model Tier / Materiality** | `Tier 2 (Significant)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price long-dated equity or FX derivatives where interest rate stochasticity is significant.
* **Methodology:** Heston stochastic volatility combined with Hull-White 1-Factor interest rates.
* **Performance:** Monte Carlo simulation or semi-analytical pricing (approximate).
* **Primary Limitation:** Computationally intensive due to multi-factor simulation.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Interest rate risk cannot be ignored for long-dated equity options.
* **Intended Use:**
  - Long-dated equity options (LEAPS)
  - Equity-Interest Rate Hybrids
* **User Group:** Structured Equity desk.

### 2.2 Scope and Coverage

* **Target Population:** Equity/FX Hybrids > 5 years maturity.
* **Exclusions:**
  - Short-dated vanilla options (use Heston)

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The model couples the Heston and Hull-White dynamics:

1.  **Equity:** $dS_t = r_t S_t dt + \sqrt{V_t} S_t dW_S$
2.  **Variance:** $dV_t = \kappa (\theta - V_t) dt + \xi \sqrt{V_t} dW_V$
3.  **Rate:** $dr_t = (\theta_r(t) - a r_t) dt + \sigma_r dW_r$

With correlations $\rho_{SV}, \rho_{Sr}, \rho_{Vr}$.

### 3.2 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Constant Correlations | Simplification | **Medium:** Decorrelation risk | Monitor cross-asset corr | **Amber** |
| A.2 | Independence of Vol/Rate | Often assumed $\rho_{Vr}=0$ | **Low:** Second order effect | N/A | **Green** |

---

## 4. Model Development

### 4.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  hybrid_hw.py: HybridHWModel                               │
│  ├── fair_value(maturity, strike)                          │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  native.HHW                                                │
│  ├── price() - Pricing engine                              │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Appendices

### A. API Reference

```python
from velesquant.models import HybridHWModel

model = HybridHWModel(
    s0=100.0,      # Spot
    v0=0.04,       # Initial Variance
    r0=0.03,       # Initial Short Rate
    kappa=2.0,     # Heston Mean Reversion
    eta=0.04,      # Heston Long Run Var (Theta)
    rho=-0.7,      # Spot-Vol Correlation
    sigma1=0.3,    # Heston Vol-of-Vol (Xi)
    sigma2=0.01,   # Hull-White Volatility
    a=0.03         # Hull-White Mean Reversion
)

# Price a European Call option under the Hybrid model
price = model.fair_value(maturity=5.0, strike=100.0)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [hybrid_hw.py](../../src/velesquant/models/hybrid_hw.py) |
| C++ Implementation | `HHW` in native module |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
