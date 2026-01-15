"""
Benchmark script for measuring Python binding overhead vs C++ performance.

This script measures the time taken to calibrate SABR and Heston models
over multiple iterations to assess binding overhead.
"""

import time

import numpy as np

import velesquant.native as n


def benchmark_sabr_calibration(iterations: int = 100):
    """Benchmark SABR calibration performance."""
    print(f"\n=== SABR Calibration Benchmark ({iterations} iterations) ===")

    # Setup test data
    strikes = np.array([[90.0, 95.0, 100.0, 105.0, 110.0]])
    quotes = np.array([[0.25, 0.22, 0.20, 0.22, 0.25]])

    # Warm-up run
    sabr = n.Sabr(maturity=1.0, forward=100.0, beta=0.5, alpha=0.5, nu=0.5, rho=-0.5)
    sabr.calibrate(
        strikes=strikes, quotes=quotes, quote_type=n.CalibrationTarget.Volatility
    )

    # Timed runs
    start_time = time.perf_counter()
    for _ in range(iterations):
        sabr = n.Sabr(
            maturity=1.0, forward=100.0, beta=0.5, alpha=0.5, nu=0.5, rho=-0.5
        )
        sabr.calibrate(
            strikes=strikes, quotes=quotes, quote_type=n.CalibrationTarget.Volatility
        )
    end_time = time.perf_counter()

    total_time = end_time - start_time
    avg_time_ms = (total_time / iterations) * 1000

    print(f"Total time: {total_time:.4f} seconds")
    print(f"Average per calibration: {avg_time_ms:.4f} ms")
    print(f"Calibrations per second: {iterations / total_time:.2f}")

    return avg_time_ms


def benchmark_heston_calibration(iterations: int = 50):
    """Benchmark Heston calibration performance."""
    print(f"\n=== Heston Calibration Benchmark ({iterations} iterations) ===")

    # Setup test data (volatility surface)
    maturities = np.array([[0.25], [0.5], [1.0]])
    forwards = np.array([[100.0], [100.0], [100.0]])
    strikes = np.array([[90.0, 100.0, 110.0]])
    quotes = np.array(
        [
            [0.25, 0.20, 0.22],
            [0.24, 0.19, 0.21],
            [0.23, 0.18, 0.20],
        ]
    )

    # Warm-up run
    heston = n.Heston(
        spot=100.0, var0=0.04, kappa=1.0, theta=0.04, xi=0.1, rho=-0.5, seed=42
    )
    heston.calibrate(maturities, forwards, strikes, quotes, "Volatility")

    # Timed runs
    start_time = time.perf_counter()
    for _ in range(iterations):
        heston = n.Heston(
            spot=100.0, var0=0.04, kappa=1.0, theta=0.04, xi=0.1, rho=-0.5, seed=42
        )
        heston.calibrate(maturities, forwards, strikes, quotes, "Volatility")
    end_time = time.perf_counter()

    total_time = end_time - start_time
    avg_time_ms = (total_time / iterations) * 1000

    print(f"Total time: {total_time:.4f} seconds")
    print(f"Average per calibration: {avg_time_ms:.4f} ms")
    print(f"Calibrations per second: {iterations / total_time:.2f}")

    return avg_time_ms


def benchmark_hull_white_simulation(iterations: int = 1000):
    """Benchmark Hull-White simulation performance."""
    print(f"\n=== Hull-White Simulation Benchmark ({iterations} iterations) ===")

    # Setup
    hw = n.HullWhite(
        kappa=0.1,
        timeSigmas=[1.0, 5.0, 10.0],
        sigmas=[0.01, 0.012, 0.015],
        timeDFs=[0.0, 1.0, 5.0, 10.0, 30.0],
        DFs=[1.0, 0.95, 0.80, 0.65, 0.30],
    )

    times = [0.25 * i for i in range(1, 121)]  # 30 years quarterly

    # Warm-up
    _ = hw.simulation(times)

    # Timed runs
    start_time = time.perf_counter()
    for _ in range(iterations):
        hw.simulation(times)
    end_time = time.perf_counter()

    total_time = end_time - start_time
    avg_time_ms = (total_time / iterations) * 1000

    print(f"Total time: {total_time:.4f} seconds")
    print(f"Average per simulation: {avg_time_ms:.4f} ms")
    print(f"Simulations per second: {iterations / total_time:.2f}")

    return avg_time_ms


if __name__ == "__main__":
    print("=" * 60)
    print("VelesQuant Python Binding Performance Benchmark")
    print("=" * 60)

    sabr_time = benchmark_sabr_calibration(iterations=100)
    heston_time = benchmark_heston_calibration(iterations=50)
    hw_time = benchmark_hull_white_simulation(iterations=1000)

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"SABR Calibration:       {sabr_time:.4f} ms/call")
    print(f"Heston Calibration:     {heston_time:.4f} ms/call")
    print(f"Hull-White Simulation:  {hw_time:.4f} ms/call")
    print("=" * 60)
