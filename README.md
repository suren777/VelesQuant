# VelesQuant

**VelesQuant** is a high-performance financial valuation library written in C++ and exposed to Python using **pybind11**. It provides tools for pricing complex financial derivatives using advanced models and numerical techniques.

## Features

- **Models**:
  - **Local Volatility**: SABR, Heston, Schobel-Zhu.
  - **Fixed Income**: Hull-White (1F/2F).
  - **Exotics**: CMS, Swaptions, Quantos.
- **Solvers**:
  - Finite Difference Methods (PDE solvers).
  - Monte Carlo Simulation.
  - Skew-adjusted simulations.

## Documentation

For detailed documentation, please refer to:
*   [**Python API Reference**](docs/python_api.md)
*   [**Pricing Model Documentation**](docs/models/README.md) (Regulatory compliant)

For practical examples, see the [tests](tests/) directory.

## Prerequisites

To build `velesquant`, you need the following system dependencies:

- **C++ Compiler** (supporting C++17)
- **CMake** (>= 3.15)
- **QuantLib** & **Boost**
- **Poetry** (for Python dependency management)

### macOS Installation
```bash
brew install quantlib boost cmake libomp
brew install poetry
```

## Installation & Build

### Python Package (Recommended)
This project uses **Poetry** and **scikit-build-core**.

```bash
# Install dependencies and the project
poetry install

# Alternatively, for development (editable mode)
poetry run pip install -e .
```

### C++ Standalone (For Core Development)
```bash
mkdir build && cd build
cmake ..
make -j4
```

## Usage

```python
import velesquant._core as v_core

# Initialize a SABR model
# sabr(maturity, forward, beta, alpha, nu, rho)
model = v_core.Sabr(1.0, 0.05, 0.5, 0.2, 0.3, -0.2)

# Calculate implied volatility
strike = 0.05
iv = model.impliedVol(strike)
print(f"SABR ATM Implied Vol: {iv:.5f}")
```

## Testing

### Python Tests
```bash
poetry run pytest
```

### C++ Tests
```bash
./build/tests_cpp/velesquant_tests
```

## Contributing

We welcome contributions! Please follow these steps:

1.  **Bug Reports**: Open an issue describing the bug and providing a reproducible example.
2.  **Feature Requests**: Use the issue tracker to propose new features.
3.  **Pull Requests**: 
    - Fork the repository.
    - Create a feature branch.
    - Ensure all C++ and Python tests pass.
    - Submit a PR with a clear description of changes.

## License
[Insert License Here]
