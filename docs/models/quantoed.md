# Quantoed CMS & Spread Models

**Regulatory Compliance:** PRA SS1/23 (UK) & SR 11-7 (US)

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Quantoed CMS Pricing Model |
| **Model ID** | `VQ-IR-005` |
| **Model Tier / Materiality** | `Tier 2 (Significant)` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price CMS derivatives settled in a currency different from the underlying swap.
* **Methodology:** Convexity adjustment augmented with Quanto adjustment (FX correlation).
* **Performance:** Closed-form adjustment.
* **Primary Limitation:** Assumes log-normal FX-Rate covariance.

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Valuation of multi-currency interest rate products.
* **Intended Use:**
  - Quanto CMS Caplets
  - Quanto CMS Spreads
* **User Group:** Cross-currency desk.

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

The drift of the underlying swap rate is adjusted for the covariance with the FX rate:

$$
\mu_{quanto} = \mu_{domestic} - \rho_{FX, Rate} \sigma_{FX} \sigma_{Rate}
$$

This adjustment is applied on top of the standard CMS convexity adjustment.

### 3.2 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | Constant FX Correlation | Standard approximation | **Medium:** FX/Rate decorrelation risk | Monitor cross-gamma | **Amber** |

---

## 4. Model Development

### 4.1 Implementation Architecture

Wraps `native.QuantoedCMS` and `native.QuantoedCmsSpread`.

### A. API Reference

```python
from velesquant.models import QuantoedCMSModel

# ... (Initialization similar to CMSModel with added FX params) ...
# FX Volatility and Correlation
sigma_fx = 0.1
corr_fx_rate = 0.3

# (Hypothetical usage - derived from test_quantoed_cms.py)
# model = QuantoedCMSModel(..., sigma_fx, corr_fx_rate)
# price = model.fair_value(strike, option_type)
```

### B. Source Files

| Component | File |
|-----------|------|
| Python Wrappers | [quantoed_cms.py](../../src/velesquant/models/quantoed_cms.py) |
| C++ Implementation | `QuantoedCMS` in native module |

---

### Final Check

- [x] **Replicability:** Model can be reproduced from this document
- [ ] **Version Control:** Update doc version when code changes
- [x] **Independence:** Documentation is objective and critical
