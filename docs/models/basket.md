# Log-Normal Basket Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Log-Normal Basket Model |
| **Model ID** | `VQ-EQ-003` |
| **Model Tier / Materiality** | `Tier 2 (Significant)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price options on a basket of assets (e.g., equity basket).
* **Methodology:** Multivariate Log-Normal dynamics with constant correlation.
* **Performance:** Standard closed-form (Gentle approximation) or Monte Carlo.
* **Primary Limitation:** Assumes constant correlation and volatility structure.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Valuation of multi-asset structures where correlation is key.
* **Intended Use:**
  - Basket options
  - Worst-of / Best-of options (via simulation)
  - Portfolio risk modeling
* **User Group:** Equity derivatives traders, structurers.

### 2.2 Scope and Coverage

* **Target Population:** Equity baskets, FX baskets.
* **Exclusions:**
  - Local volatility baskets (use SLV)
  - Stochastic correlation models

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The model assumes each asset $S_i$ follows geometric Brownian motion:

$$
dS_i(t) = (r - q_i) S_i(t) dt + \sigma_i S_i(t) dW_i(t)
$$

With correlation structure:
$$
d W_i(t) d W_j(t) = \rho_{ij} dt
$$

Simulation is performed using Cholesky decomposition of the correlation matrix to generate correlated random numbers.

### 3.2 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Constant Correlation | Tractability | **High:** Mispricing in stress | Monitor skew sensitivity | **Amber** |
| A.2 | Log-normal dynamics | Standard approximation | **Medium:** Fat tails missing | Compare to market skew | **Amber** |

---

## 4. Model Development

### 4.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  basket.py: LogNormalBasketModel                           │
│  ├── simulate(schedule)                                    │
│  └── simulate_with_rebalancing(schedule)                   │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  native.LogNormalBasket                                    │
│  ├── simulate() - Monte Carlo engine                       │
│  └── get_n_assets()                                        │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Appendices

### A. API Reference

```python
from velesquant.models import LogNormalBasketModel

# Define basket parameters (2 assets)
spots = [100.0, 100.0]
strikes = [100.0, 100.0]
maturities = [1.0, 1.0]

# Forwards and IVs are lists of lists (term structure per asset)
forwards = [[100.0, 100.0], [100.0, 100.0]]
ivs = [[0.2, 0.2], [0.2, 0.2]]

# Correlation Matrix
corr = [
    [1.0, 0.5],
    [0.5, 1.0]
]

model = LogNormalBasketModel(
    spots=spots,
    strikes=strikes,
    maturities=maturities,
    forwards=forwards,
    ivs=ivs,
    correlation=corr
)

# Run simulation
paths = model.simulate(schedule=[0.5, 1.0])
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [basket.py](../../src/velesquant/models/basket.py) |
| C++ Implementation | `LogNormalBasket` in native module |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
