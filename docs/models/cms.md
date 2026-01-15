# CMS (Constant Maturity Swap) Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | CMS Convexity-Adjusted Pricing Model |
| **Model ID** | `VQ-IR-002` |
| **Model Tier / Materiality** | `Tier 2 (Significant)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price CMS caplets, floorlets, and CMS-linked structures with convexity adjustment.
* **Methodology:** SABR-based smile with convexity-adjusted forward rate.
* **Performance:** Accurate CMS pricing across strike range.
* **Primary Limitation:** Convexity adjustment is approximation-based.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** CMS rates exhibit convexity bias requiring adjustment.
* **Intended Use:**
  - CMS caplet/floorlet pricing
  - CMS swap valuation
  - CMS spread products (via extension)
* **User Group:** Interest rate derivatives traders, structurers.

### 2.2 Scope and Coverage

* **Target Population:** CMS-linked derivatives.
* **Exclusions:**
  - Path-dependent CMS exotics (use simulation)
  - CMS spread options (use CMS Spread model)

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

**CMS Convexity Adjustment:**

The CMS rate is a non-linear function of the underlying swap rate. Under the CMS measure, the forward CMS rate differs from the forward swap rate:

$$
E^{CMS}[S_T] = S_0 + \text{Convexity Adjustment}
$$

The adjustment depends on:
- Swap rate volatility (ATM vol from SABR)
- Duration and convexity of the underlying swap
- Payment timing mismatch

**Implementation:**
1. Calibrate SABR to swaption smile
2. Compute convexity-adjusted forward
3. Price CMS options using SABR with adjusted forward

**Justification:** Standard market practice for CMS; balances accuracy with computational efficiency.

### 3.2 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | SABR dynamics | Industry standard | **Low:** Well-tested | N/A | **Green** |
| A.2 | Convexity approximation | Quick pricing | **Medium:** Bias for long tenors | Compare to replication | **Amber** |
| A.3 | Static adjustment | No dynamic hedging | **Low:** Conservative | N/A | **Green** |

---

## 4. Model Development

### 4.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  cms.py: CMSModel                                          │
│  ├── adjusted_forward (property)                           │
│  ├── fair_value(strike, option_type)                       │
│  └── implied_vol(strike)                                   │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  cms.cpp: cms class                                        │
│  ├── fairValue() - CMS option pricing                      │
│  ├── convAdj() / convAdjAlt() - Convexity adjustment       │
│  └── Internal SABR calibration                             │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Appendices

### A. API Reference

```python
from velesquant.models import CMSModel

model = CMSModel(
    expiry=1.0,
    tenor=10.0,
    forward=0.03,
    annuity=8.5,
    pay_cms=1.0,
    discount_cms=0.97,
    beta=0.85,
    alpha=0.5,
    nu=0.25,
    rho=-0.75
)

# Get convexity-adjusted forward
adj_fwd = model.adjusted_forward

# Price a CMS caplet
price = model.fair_value(strike=0.035, option_type="Call")
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [cms.py](../../src/velesquant/models/cms.py) |
| C++ Implementation | [cms.cpp](../../src/models/cms.cpp) |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
