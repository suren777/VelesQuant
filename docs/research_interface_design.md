# Research Report: State-of-the-Art Quant Library Interface

## 1. Executive Summary

The current **VelesQuant** library provides a solid mathematical foundation with high-performance C++ implementations of key models (SABR, Hull-White, Heston). However, the current interface is **Model-Centric**, where pricing logic is tightly coupled within the model classes themselves (e.g., `HullWhite::swaption`).

**State-of-the-Art** libraries (QuantLib, ORE, Strata) utilize an **Instrument-Centric** design, where:
1.  **Instruments** (Swaption, Bond) describe *what* is being traded (payoffs, dates).
2.  **Models** (Hull-White, Black-Scholes) describe *how* market variables yield.
3.  **Pricing Engines** link the two, creating a flexible grid where new instruments can be priced by existing models without code duplication.

We propose a **"Veles Interface" Layer** (primarily in Python) to bridge this gap, wrapping the existing C++ foundation into a modern, unifed API.

---

## 2. Current Architecture vs. Best Practices

### A. Current State (VelesQuant)

*   **Pattern**: Monolithic "Model as Pricer".
*   **Example**: To price a swaption, you must instantiate a `HullWhite` object and call `swaption(...)`.
*   **Pros**: Simple to implement, low abstraction overhead, direct control.
*   **Cons**:
    *   **High Coupling**: The Model class must know about every instrument type it supports.
    *   **Inconsistent API**: `Sabr` has `impliedVol`, `HullWhite` has `optionBond`, `sabr_facade` handles others. No unifying `price()` method.
    *   **Hard to Swap Models**: Changing from Black-Scholes to Heston requires changing the calling code significantly.

### B. Industry Best Practices (QuantLib / ORE)

*   **Pattern**: Instrument / Engine Separation.
*   **Structure**:
    *   `Instrument`: Lightweight data container (e.g., `EuropeanOption` holds Payoff and Exercise). Implementation-agnostic.
    *   `PricingEngine`: The "Visiting" logic. A `AnalyticEuropeanEngine` knows how to price a `EuropeanOption` using Black-Scholes. A `MCEuropeanEngine` knows how to use Monte Carlo.
    *   `MarketData`: Handles (smart pointers) to Curves and Vol Surfaces, shared across instruments.

---

## 3. Proposed Design: The "Veles Interface"

We recommend building a Pythonic high-level API that abstracts the C++ complexity.

### 3.1. Unified Workflow

The goal is to enable code like this:

```python
# 1. Define Market Data (Wraps C++ Curves/Surfaces)
market = MarketEnvironment("2024-01-01")
market.add_curve("USD-LIBOR-3M", DiscountCurve(...))

# 2. Define Instrument (Pure Data)
swaption = Swaption(
    expiry="5Y", 
    tenor="10Y", 
    strike=0.03, 
    receiver=True
)

# 3. Choose Model/Engine
model = HullWhiteModel(kappa=0.03, sigma=0.01)

# 4. Price
price = model.value(swaption, market)
# OR, in a more QuantLib style:
# swaption.set_pricing_engine(HullWhiteEngine(model, market))
# price = swaption.NPV()
```

### 3.2. Implementation Strategy

1.  **Layer 1: Core Bindings (Existing)**
    *   Keep the existing `velesquant.native` as the low-level calculation engine.
    *   No major C++ rewrites required immediately.

2.  **Layer 2: Python Abstraction (New)**
    *   Create `velesquant.instruments` module (Classes: `Swap`, `Swaption`, `Bond`, `CapFloor`).
    *   Create `velesquant.models` module (wrappers around `native.HullWhite`, `native.Sabr`).
    *   Create `velesquant.pricing` module to handle the dispatch logic.

3.  **Layer 3: The Dispatcher**
    *   Implement a Visitor or Registry pattern in Python that maps `(InstrumentType, ModelType) -> PricingFunction`.
    *   *Example*: `(Swaption, HullWhiteModel) -> native.HullWhite(...).swaption(...)`

---

## 4. Key Design Decisions

| Feature | Recommendation | Rationale |
| :--- | :--- | :--- |
| **Data Types** | Use standard Python `datetime.date` and `enums`. | Improves usability over raw doubles/ints. |
| **Arrays/Matrices** | Use `numpy` (already in progress with Eigen). | Standard for Python data science stack. |
| **Market Data** | Context Object. | Pass a `Market` object to pricers rather than individual curves. Allows easy scenario analysis (bumping curves). |
| **Calibration** | Model methods. | `model.calibrate(instruments, market_data)` should return a *new* calibrated model instance (Generic functional pattern) or update state. |

## 5. Next Steps

1.  **Prototype**: Create the `Instrument` base class in Python.
2.  **Implement**: `Swaption` Python class and `HullWhite` Python wrapper.
3.  **Bridge**: Write the glue code that extracts data from the Python `Swaption`, configures the C++ `native.HullWhite`, and returns the result.


## 6. Verification & Refinement

We validated this "Python-side Dispatch" architecture against industry standards (Perplexity/QuantLib analysis).

### Findings
1.  **Usability**: The "Data in Python, Math in C++" split is highly favored for modern data science workflows (similar to PyTorch/TensorFlow). It offers a superior User Experience compared to exposing raw C++ factories.
2.  **Performance Risk**: A potential pitfall is the overhead of crossing the Python/C++ boundary for *every single instrument* in a loop.
    *   *Bad*: `for s in 10000_swaptions: model.price(s)`
    *   *Good*: `model.price_batch(list_of_swaptions)` or Vectorized inputs.

### Refined Recommendation
We will proceed with the **Instrument-Centric Python Interface**, but with a mandate to support **Vectorization** in the future.
*   **Phase 1 (Now)**: Implement the logical separation. `Swaption` object -> `HullWhite` scalar pricing.
*   **Phase 2 (Future)**: enhance `native` bindings to accept `Eigen::VectorXd` for expiries/strikes, allowing a single Python call to price a portfolio.

