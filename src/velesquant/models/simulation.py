from typing import List, Union
import numpy as np
from velesquant.models.localvol import LocalVolModel


def simulate_multi_asset(
    models: List[LocalVolModel],
    correlation_matrix: Union[List[List[float]], np.ndarray],
    times: Union[List[float], np.ndarray],
    n_paths: int,
    seed: int = 42,
) -> np.ndarray:
    """
    Simulate multiple correlated assets using Local Volatility models.

    Parameters:
    - models: List of LocalVolModel instances.
    - correlation_matrix: N x N correlation matrix.
    - times: Simulation time points.
    - n_paths: Number of Monte Carlo paths.
    - seed: Random seed.

    Returns:
    - np.ndarray: Shape (n_paths, len(times), n_assets) containing spot paths.
                  Or (n_assets * len(times), n_paths) to match C++ matrix layout if needed.
                  C++ facade returned Matrix(n*N_assets, NP). Let's return (n_paths, n_assets, n_steps) for Pythonic usage.
    """
    n_assets = len(models)

    # Validation
    corr = np.array(correlation_matrix)
    if corr.shape != (n_assets, n_assets):
        raise ValueError(
            f"Correlation matrix shape {corr.shape} must match number of models {n_assets}"
        )

    times_arr = np.array(times, dtype=np.float64)
    n_steps = len(times_arr)

    # Cholesky decomposition for correlation
    # corr = L * L.T
    try:
        L = np.linalg.cholesky(corr)
    except np.linalg.LinAlgError:
        raise ValueError("Correlation matrix is not positive definite")

    rng = np.random.default_rng(seed)

    # Generate uncorrelated random numbers: (n_steps, n_paths, n_assets)
    # We need n_steps random numbers per path per asset.
    # The C++ code generates randoms step-by-step or all at once?
    # C++: for each path of size n, generates n rands per asset.

    # Shape: (n_assets, n_paths, n_steps)
    uncorrelated_rands = rng.standard_normal((n_assets, n_paths, n_steps))

    # Apply correlation
    # We want correlated rands C such that Cov(C) = corr.
    # If Z is uncorrelated (n_assets, ...), then C = L @ Z
    # reshape for matmul: (n_assets, n_paths * n_steps)
    uncorrelated_reshaped = uncorrelated_rands.reshape(n_assets, -1)
    correlated_reshaped = L @ uncorrelated_reshaped
    correlated_rands = correlated_reshaped.reshape(n_assets, n_paths, n_steps)

    # Run simulations
    # Result container: (n_paths, n_assets, n_steps)
    results = np.zeros((n_paths, n_assets, n_steps))

    for i, model in enumerate(models):
        # We need to simulate path by path or vectorised?
        # native model.simulate takes (times, rands) for ONE path?
        # Let's check lVol::simulation signature.
        # std::vector<double> simulation(std::vector<double> times, std::vector<double> rands) const;
        # It returns vector of spots for ONE path.

        # This loop is slow in Python. But lVol doesn't expose vectorised simulation across paths yet?
        # Simulation facade did `modelLV->simulation(times, rands1)` inside a loop over NP.
        # So we emulate that layer in Python.

        # Optimization: We could push this loop to C++ later if needed (expose a batch simulation).
        # For now, strictly refactoring logic to Python.

        asset_rands = correlated_rands[i]  # (n_paths, n_steps)

        for p in range(n_paths):
            path_rands = asset_rands[p, :]
            # Convert to list for binding if needed, or binding handles numpy array
            path_spots = model.simulate(times_arr, path_rands)
            results[p, i, :] = path_spots

    return results
