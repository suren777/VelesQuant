# SABR Swaption Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | SABR-based Swaption Pricing Model |
| **Model ID** | `VQ-IR-003` |
| **Model Tier / Materiality** | `Tier 1 (Critical)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price European swaptions with SABR smile dynamics.
* **Methodology:** SABR volatility model applied to forward swap rates.
* **Performance:** Accurate swaption pricing across strike range.
* **Primary Limitation:** European exercise only (use HW PDE for Bermudans).

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Swaption smile requires SABR parameterization.
* **Intended Use:**
  - Payer/receiver swaption pricing
  - Swaption vol surface calibration
  - Delta/vega hedging
* **User Group:** Interest rate derivatives traders.

### 2.2 Scope and Coverage

* **Target Population:** European swaptions.
* **Exclusions:**
  - Bermudan/callable (use Hull-White PDE)
  - CMS-linked swaptions (use CMS model)

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The model applies SABR dynamics to the forward swap rate:

$$
dS_t = \alpha_t S_t^\beta dW_t^S
$$

With swaption value:
$$
V = A \cdot \text{Black}(S, K, \sigma_{SABR}(K), T)
$$

Where $A$ is the annuity (PV01) of the underlying swap.

---

## 4. Model Development

### 4.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  swaption_sabr.py: SabrSwaptionModel                       │
│  ├── fair_value(strike, option_type)                       │
│  ├── implied_vol(strike)                                   │
│  └── swap_value(strike)                                    │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  swaption.cpp: Swaption class                              │
│  ├── swaptionFairValue() - Black pricing with SABR vol     │
│  ├── swapFairValue() - Underlying swap value               │
│  └── getImpliedVol() - SABR implied vol                    │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Appendices

### A. API Reference

```python
from velesquant.models import SabrSwaptionModel

model = SabrSwaptionModel(
    expiry=1.0,
    tenor=5.0,
    forward=0.025,
    annuity=4.5,
    beta=0.85,
    alpha=0.5,
    nu=0.25,
    rho=-0.75
)

# Price a payer swaption
price = model.fair_value(strike=0.03, option_type="Call")

# Get implied volatility
vol = model.implied_vol(strike=0.03)

# Underlying swap value
swap = model.swap_value(strike=0.03)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [swaption_sabr.py](../../src/velesquant/models/swaption_sabr.py) |
| C++ Implementation | [swaption.cpp](../../src/models/swaption.cpp) |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
