# CMS Spread Model

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | CMS Spread Option Pricing Model |
| **Model ID** | `VQ-IR-004` |
| **Model Tier / Materiality** | `Tier 2 (Significant)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price options on the spread between two constant maturity swap rates (e.g., CMS10 - CMS2).
* **Methodology:** Copula-based approach combining two marginal CMS distributions.
* **Performance:** Efficient integration of joint density.
* **Primary Limitation:** Gaussian copula assumption.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Valuation of steepeners/flatteners.
* **Intended Use:**
  - CMS Spread options
  - Range accruals on spreads
* **User Group:** Interest rate derivatives traders.

### 2.2 Scope and Coverage

* **Target Population:** Spread options on major indices (EURIBOR, SOFR).
* **Exclusions:**
  - Path-dependent spread exotics (use multi-factor HWPDE)

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The model constructs the joint distribution of two CMS rates $S_1, S_2$:

1.  **Marginals:** Derived from SABR-convexity-adjusted distributions for individual CMS rates.
2.  **Dependence:** Modeled via a Gaussian Copula with correlation $\rho$.

The option price is the expectation of the payoff under this joint measure:
$$
V = E[(\alpha S_1 - \beta S_2 - K)^+]
$$

This is solved via numerical integration over the joint density.

### 3.2 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Gaussian Copula | Standard market model | **Medium:** Tail dependence missed | Compare to historical spread | **Amber** |
| A.2 | Constant Correlation | Simplification | **Medium:** Spread volatility risk | Monitor correlation sensitivity | **Amber** |

---

## 4. Model Development

### 4.1 Implementation Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  cms_spread.py: CMSSpreadModel                             │
│  ├── spread_option(K, a, b)                                │
│  └── simulate()                                            │
└─────────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  native.CmsSpread                                          │
│  ├── spread_option() - Integration engine                  │
│  └── simulate()                                            │
└─────────────────────────────────────────────────────────────┘
```

---

## 5. Appendices

### A. API Reference

```python
from velesquant.models import CMSSpreadModel
from velesquant.native import CalibrationTarget

# Define two legs (CMS10y and CMS2y parameters)
# ... [Parameter setup omitted for brevity, see test cases] ...

model = CMSSpreadModel(
    # Leg 1 (e.g. 10y)
    expiry1=1.0, tenor1=10.0, fwd1=0.03, annuity1=8.5, pay1=1.0, disc1=0.97, beta1=0.5,
    strikes1=[0.02, 0.04], quotes1=[0.20, 0.25], type1=CalibrationTarget.Volatility,
    # Leg 2 (e.g. 2y)
    expiry2=1.0, tenor2=2.0, fwd2=0.02, annuity2=1.9, pay2=1.0, disc2=0.97, beta2=0.5,
    strikes2=[0.01, 0.03], quotes2=[0.25, 0.30], type2=CalibrationTarget.Volatility,
    # Correlation
    corr=0.6
)

# Price spread option: Payoff = (1.0 * CMS10 - 1.0 * CMS2 - 0.005)+
price = model.spread_option(K=0.005, a=1.0, b=1.0)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrapper | [cms_spread.py](../../src/velesquant/models/cms_spread.py) |
| C++ Implementation | `CmsSpread` in native module |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
