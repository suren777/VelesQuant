import time

import numpy as np
import velesquant as n


def benchmark():
    print("Starting HWPDE Benchmark...")
    # Setup model
    kappa = 0.05
    time_sigmas = [1.0, 5.0, 10.0, 20.0]
    sigmas = [0.01, 0.012, 0.015, 0.010]
    time_dfs = [0.0, 1.0, 5.0, 10.0, 20.0, 30.0]
    # Simple exponential discount factors for setup
    dfs = [np.exp(-0.03 * t) for t in time_dfs]

    hw_model = n.HullWhiteModel(kappa, time_sigmas, sigmas, time_dfs, dfs)

    # Use standard grid parameters
    grid_points = 512
    time_step = 0.01  # 100 steps per year

    print(f"Grid: {grid_points} points, dt: {time_step}")
    pde = n.HWPDE(hw_model, grid_points, time_step)

    # Warmup
    pde.price_swap(0.0, 5.0, 0.03, 0.5)

    iterations = 50
    start = time.time()
    for i in range(iterations):
        # Price a 10Y Swap (more steps)
        pde.price_swap(0.0, 10.0, 0.03, 0.5)

    end = time.time()
    avg_time = (end - start) / iterations
    print(f"Time per simple swap pricing (10Y): {avg_time:.4f}s")
    print(f"Total time for {iterations} iterations: {end-start:.4f}s")

    # Benchmark Calibration (more intensive)
    print("\nBenchmarking Calibration (10 instruments)...")
    swaps = []
    for t in range(1, 11):
        s = n.DefSwap()
        s.expiry = float(t)
        s.tenor = 5.0  # 5Y into 1Y, 2Y... 10Y
        s.frequency = 0.5
        s.swap_rate = 0.03 + 0.001 * t
        s.vol_atm = 0.20
        s.value = 0.0
        swaps.append(s)

    start_cal = time.time()
    try:
        pde.calibrate(time_dfs, dfs, swaps)
    except RuntimeError as e:
        print(f"Calibration finished (or failed with expected error): {e}")

    end_cal = time.time()
    print(f"Calibration time: {end_cal - start_cal:.4f}s")


if __name__ == "__main__":
    benchmark()
