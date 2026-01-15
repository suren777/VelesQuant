# Local Volatility Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Local Volatility Model (Dupire) |
| **Model ID** | `VQ-VOL-003` |
| **Model Tier / Materiality** | `Tier 2 (Significant)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Construct a state-dependent volatility surface for PDE-based pricing.
* **Methodology:** Dupire local volatility extracted from SABR smile slices.
* **Performance:** Exact fit to vanilla prices by construction.
* **Primary Limitation:** Forward smile dynamics not captured (static model).

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Price path-dependent exotics consistent with vanilla market.
* **Intended Use:**
  - PDE pricing of barriers, Asians
  - Risk-neutral density extraction
  - Calibration target for hybrid models
* **User Group:** Exotic derivatives traders, quantitative analysts.

### 2.2 Scope and Coverage

* **Target Population:** Equity/FX options where smile consistency is critical.
* **Exclusions:**
  - Products sensitive to forward smile (use SABR/Heston)

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

**Dupire Local Volatility:**

Given European call prices $C(K,T)$, local volatility is:

$$
\sigma_{loc}^2(K,T) = \frac{\frac{\partial C}{\partial T}}{\frac{1}{2}K^2 \frac{\partial^2 C}{\partial K^2}}
$$

**Implementation:**
- Input: SABR models per maturity slice
- Extract local vol via finite difference on SABR call prices
- PDE solver for option pricing

**Justification:** Standard approach for consistent exotic pricing.

### 3.2 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Diffusion dynamics | No jumps | **Medium:** Fat tails missed | N/A | **Amber** |
| A.2 | Static surface | No re-calibration | **High:** Forward smile wrong | Compare to stochastic vol | **Red** |

---

## 4. Model Development

### 4.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  localvol.py: LocalVolModel                                │
│  ├── call_pde(maturity, strike, num_steps)                 │
│  ├── put_pde(maturity, strike, num_steps)                  │
│  └── density(maturity, num_steps)                          │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  l_vol.cpp / local_vol_facade.cpp                          │
│  ├── callPDE() / putPDE() - PDE pricing                    │
│  └── density() - Risk-neutral density                      │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Appendices

### A. API Reference

```python
from velesquant.models import LocalVolModel, SabrModel

# Build SABR slices
sabr_1y = SabrModel(maturity=1.0, forward=100, beta=0.5, alpha=0.2, nu=0.4, rho=-0.3)
sabr_2y = SabrModel(maturity=2.0, forward=102, beta=0.5, alpha=0.22, nu=0.38, rho=-0.35)

# Construct local vol model
model = LocalVolModel(sabr_models=[sabr_1y, sabr_2y], spot=100.0)

# Price via PDE
call_price = model.call_pde(maturity=1.5, strike=105, num_steps=100)
put_price = model.put_pde(maturity=1.5, strike=95, num_steps=100)

# Get density
density = model.density(maturity=1.0, num_steps=50)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [localvol.py](../../src/velesquant/models/localvol.py) |
| C++ Implementation | [l_vol.cpp](../../src/volatility/l_vol.cpp) |
| C++ Facade | [local_vol_facade.cpp](../../src/models/local_vol_facade.cpp) |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
