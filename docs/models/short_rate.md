# Short Rate Models (1F & 2F)

**Regulatory Compliance:** [Relevant Regulation if known, e.g., SR 11-7]

---

## Document Control & Metadata

| Field | Details |
| --- | --- |
| **Model Name** | Normalized Short Rate Models (1-Factor & 2-Factor/G2++) |
| **Model ID** | `VQ-IR-002` |
| **Model Tier / Materiality** | `Tier 2` |
| **Model Owner** | `[TBD]` |
| **Model Developer(s)** | VelesQuant Team |
| **Version** | `v1.0` |
| **Last Approval Date** | `[TBD]` |
| **Review Frequency** | Annual |

---

## 1. Executive Summary

* **Purpose:** Price and hedge interest rate derivatives using generalized short rate dynamics. Ideally suited for products requiring time-dependent volatility or multi-factor correlation (e.g., swaptions, options on bonds).
* **Methodology:** 
    *   **1-Factor:** Generalized PDE solver for models of the form $dr = (\theta - \kappa r)dt + \sigma r^\gamma dW$ (includes Black-Karasinski, Extended Vasicek).
    *   **2-Factor (G2++):** Two-factor additive Gaussian model solved via ADI (Alternating Direction Implicit) PDE scheme.
* **Performance:** 
    *   **1F:** Efficient for Bermudan swaptions and standard American products.
    *   **2F:** Captures decorrelation and humped volatility structures; computationally more intensive.
* **Primary Limitation:** Grid-based PDE methods scale poorly beyond 2-3 factors (Curse of Dimensionality).

---

## 2. Model Scope and Usage

### 2.1 Business Purpose

* **Problem Statement:** Valuation of interest rate products where the assumption of a constant or single-factor term structure is insufficient, or where specific dynamics (like log-normality in Black-Karasinski) are required.
* **Intended Use:** 
    *   Valuation of European and Bermudan Swaptions.
    *   Pricing Zero Coupon and Coupon Bond Options.
    *   Risk management (Greeks) via finite difference grid.
* **User Group:** Quantitative Analysts, Rate Traders, Risk Management.

### 2.2 Scope and Coverage

* **Target Population:** Vanilla and light-exotic interest rate derivatives.
* **Exclusions:** 
    *   Path-dependent options requiring full Monte Carlo (unless handled by PDE state augmentation).
    *   High-dimensional baskets (multi-asset).

---

## 3. Mathematical Theory and Methodology

### 3.1 Theoretical Framework

#### One-Factor Generalized Model (`ShortRate1FPDE`)
The model dynamics are defined by the generalized stochastic differential equation (SDE):

$$ dr_t = [\theta(t) - \kappa r_t] dt + \sigma(t) (\alpha + \beta r_t)^\gamma dW_t $$

The associated PDE for a contingent claim $V(t, r)$ is:

$$ \frac{\partial V}{\partial t} + \frac{1}{2}\sigma(t)^2 (\alpha + \beta r)^{2\gamma} \frac{\partial^2 V}{\partial r^2} + (\theta(t) - \kappa r) \frac{\partial V}{\partial r} - rV = 0 $$

Key parameters:
- $\kappa$: Mean reversion speed.
- $\sigma(t)$: Time-dependent volatility.
- $\alpha, \beta, \gamma$: Control the distribution shape (e.g., $\gamma=0$ for Normal/Vasicek, $\gamma=1$ for Log-Normal/Black-Karasinski with shift).

#### Two-Factor G2++ Model (`ShortRate2FPDE`)
The short rate is given by $r_t = x_t + y_t + \varphi(t)$, where $x_t$ and $y_t$ are correlated Gaussian factors:

$$ dx_t = -\kappa_1 x_t dt + \sigma_1 dW_1 $$
$$ dy_t = -\kappa_2 y_t dt + \sigma_2 dW_2 $$
$$ dW_1 dW_2 = \rho dt $$

This model allows for a richer correlation structure between rates of different maturities compared to 1-factor models.

### 3.2 Alternatives Considered

| Alternative | Reason for Rejection |
|-------------|---------------------|
| **Hull-White Analytic (`HWModel`)** | Faster for European options but less flexible for Bermudan/American features or non-normal dynamics. |
| **Monte Carlo Simulation** | Necessary for path-dependent options but slower and less precise for Greeks compared to PDE grid methods for low-dimensional problems. |

### 3.3 Key Assumptions Register

| ID | Assumption | Justification | Impact of Failure | Monitoring | Rating |
| --- | --- | --- | --- | --- | --- |
| A.1 | **Continuous Diffusion** | Core assumption of SDE interactions. | Large for jump-diffusion markets. | Monitor market jumps. | **Amber** |
| A.2 | **Discretization Bias** | Finite Difference / ADI schemes approximate continuous derivatives. | Convergence error in prices. | Grid convergence tests. | **Green** |
| A.3 | **Parameter Stationarity** | $\kappa, \rho$ are often assumed constant piecewise. | Misfit to long-dated term structures. | Regular calibration. | **Green** |

---

## 4. Data

### 4.1 Data Sources

*   **Yield Curves:** Discount factors for the reporting currency ($DF(t)$).
*   **Volatility Surface:** Swaption volatilities (Normal or Log-Normal) for calibration.
*   **Instrument Data:** Term sheets for swaps, bonds, and options.

### 4.2 Processing & Cleaning

*   **Yield Curves:** Construction via bootstrapping or optimization (e.g., Monotone Convex Splines) before input.
*   **Interpolation:** Linear interpolation on log-discount factors or rates within the PDE time grid.

### 4.3 Data Quality & Proxies

*   **Quality Checks:** Check for non-negative forward rates (unless mandated by market) and arbitrage-free surfaces.

---

## 5. Model Development

### 5.1 Implementation Architecture

The implementation uses a split architecture:
1.  **C++ Core (`velesquant::native`):**
    *   `ShortRate1FPDE`: Crank-Nicolson or Implicit finite difference solver.
    *   `ShortRate2FPDE`: ADI (Douglas Scheme) solver for 2D spatial grids.
2.  **Python Wrappers (`velesquant.models.pde_solvers`):**
    *   `ShortRate1FPDEModel`: User-friendly interface for 1F model.
    *   `ShortRate2FPDEModel`: User-friendly interface for G2++.

### 5.2 Parameter Estimation

**Calibration Process:**
1.  **1-Factor:** Calibrate mean reversion ($\kappa$) and volatility term structure ($\sigma(t)$) to a strip of co-terminal or diagonal swaptions.
    *   Objective: Minimize $\sum (ModelPrice - MarketPrice)^2$.
2.  **2-Factor:** Calibrate $\kappa_1, \kappa_2, \rho, \sigma_1, \sigma_2$ to reference swaptions.
    *   $\varphi(t)$ is automatically fitted to match the initial discount curve exactly.

---

## 6. Testing and Performance

### 6.1 Statistical Performance

*   **Calibration:** Levenberg-Marquardt optimizer used for stable convergence.
*   **Accuracy:** PDE schemes generally converge to $O(\Delta t^2 + \Delta x^2)$ accuracy.

### 6.2 Supported Products

| Product | Method | Validated |
|---------|--------|-----------|
| **Zero Coupon Bond** | PDE | ✅ |
| **Coupon Bond Option** | PDE | ✅ |
| **European Swaption** | PDE | ✅ |
| **Bermudan Swaption** | PDE | ✅ |

### 6.3 Benchmarking

*   **Benchmarks:** Compared against Analytical Hull-White (for 1F, $\gamma=0$) and QuantLib G2++ implementation.

---

## 7. Limitations and Post-Model Adjustments

### 7.1 Known Limitations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| **Dimensionality** | 2F is the practical limit for grid PDE. | Use Monte Carlo for >2 factors. |
| **Negative Rates** | Log-normal 1F (Black-Karasinski) fails with negative rates. | Use Shifted Log-normal or Normal (Hull-White) dynamics. |
| **Grid Boundaries** | Artificial boundaries can reflect probability mass. | Ensure grid boundaries are sufficiently wide ($>\pm 5$ std deviations). |

### 7.2 Post-Model Adjustments (Overlays)

*   None explicitly programmed in the model wrapper. Users apply shifts to input curves if needed for specific negative rate handling outside the model logic.

---

## 8. Ongoing Monitoring and Governance

### 8.1 Governance Framework

> Operates under VelesQuant general MRM (Model Risk Management) policy.
> *   **Daily:** Calibration checks by Trading/Quant team.
> *   **Periodic:** Validation of grid convergence configurations.

### 8.2 Performance Triggers

| Metric | Green | Amber | Red |
| --- | --- | --- | --- |
| **Grid Convergence** | < 1bp change | 1-5bp change | > 5bp change (doubling steps) |
| **Calibration Error** | < 1% Vol | 1-5% Vol | > 5% Vol |

---

## 9. Appendices

### A. API Reference

#### 1-Factor Usage
```python
from velesquant.models.pde_solvers import ShortRate1FPDEModel

model = ShortRate1FPDEModel(
    initial_rate=0.03, kappa=0.01, sigma=0.01,
    alpha=0.0, beta=1.0, gamma=0.0, # Hull-White / Normal Dynamics
    time_sigmas=[1.0], sigmas=[0.01]
)
price = model.price_swaption(expiry=1.0, tenor=5.0, strike=0.03)
```

#### 2-Factor Usage
```python
from velesquant.models.pde_solvers import ShortRate2FPDEModel

model = ShortRate2FPDEModel(
    kappa1=0.03, kappa2=0.01, lam=-0.5,
    time_sigma1s=[1.0], sigma1s=[0.01],
    time_sigma2s=[1.0], sigma2s=[0.005],
    time_alphas=[1.0], alphas=[0.0]
)
price = model.price_swaption(expiry=1.0, tenor=5.0, strike=0.03)
```

### B. Source Files

| Component | File |
|-----------|------|
| **Python Wrapper** | `src/velesquant/models/pde_solvers.py` |
| **C++ 1F Header** | `include/velesquant/pde_solvers/short_rate_1f_pde.h` |
| **C++ 2F Header** | `include/velesquant/pde_solvers/short_rate_2f_pde.h` |

### C. References

1.  Brigo, D. & Mercurio, F. (2006). "Interest Rate Models - Theory and Practice"
2.  Numerical Recipes in C++ (for PDE solver schemes)
3.  Hull, J. & White, A. (1990). "Pricing Interest-Rate-Derivative Securities"

---

### Final Check

- [x] **Replicability:** Model can be reproduced from SDE/PDE descriptions.
- [x] **Version Control:** Linked to specific source files.
- [x] **Independence:** Validated against benchmarks.
