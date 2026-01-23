
import numpy as np

from velesquant import LocalVol

from ..instruments.base import Instrument
from .base import MarketDataInput, Model
from .sabr import SabrModel


class LocalVolModel(Model):
    """
    Local Volatility Model Wrapper.

    Constructs a local volatility surface from a list of SABR models
    (one per maturity slice) OR from direct vectors of SABR parameters.
    """

    def __init__(
        self,
        sabr_models: list[SabrModel] | None = None,
        spot: float = 100.0,
        maturities: list[float] | np.ndarray | None = None,
        forwards: list[float] | np.ndarray | None = None,
        betas: list[float] | np.ndarray | None = None,
        alphas: list[float] | np.ndarray | None = None,
        nus: list[float] | np.ndarray | None = None,
        rhos: list[float] | np.ndarray | None = None,
    ):
        """
        Initialize LocalVolModel.

        Args:
            sabr_models: List of SabrModel instances (one per maturity).
            spot: Initial spot price.
            maturities: List of maturities (if sabr_models not provided).
            forwards: List of forwards.
            betas: List of betas.
            alphas: List of alphas.
            nus: List of nus.
            rhos: List of rhos.
        """
        self._spot = spot
        self._sabr_models = sabr_models

        if sabr_models is not None:
            # Extract native SABR objects for C++ construction
            native_sabrs = [s._cpp_model for s in sabr_models]
            self._cpp_model = LocalVol(native_sabrs)
            self._cpp_model.spot = spot
        elif maturities is not None:
            # Vector based construction
            self.maturities = self._to_flat_list(maturities)
            self.forwards = self._to_flat_list(forwards) if forwards is not None else []
            self.betas = self._to_flat_list(betas) if betas is not None else []
            self.alphas = self._to_flat_list(alphas) if alphas is not None else []
            self.nus = self._to_flat_list(nus) if nus is not None else []
            self.rhos = self._to_flat_list(rhos) if rhos is not None else []

            # Basic validation of lengths?
            # Pybind11 will throw if vector sizes mismatch in lVol constructor potentially?
            # Or lVol constructor assumes them.

            self._cpp_model = LocalVol(
                self.maturities,
                self.forwards,
                self.betas,
                self.alphas,
                self.nus,
                self.rhos,
                self._spot,
            )
        else:
            # Default empty constructor?
            self._cpp_model = LocalVol()
            self._cpp_model.spot = spot

    @staticmethod
    def _to_flat_list(
        data: list[float] | np.ndarray | list[list[float]],
    ) -> list[float]:
        if isinstance(data, np.ndarray):
            return data.flatten().tolist()
        if isinstance(data, list):
            # Check if it's a list of lists (matrix)
            if data and isinstance(data[0], list):
                return [item for sublist in data for item in sublist]
            return data
        return list(data)

    @property
    def spot(self) -> float:
        """Get current spot price."""
        return self._cpp_model.spot

    @spot.setter
    def spot(self, value: float):
        """Set spot price."""
        self._spot = value
        self._cpp_model.spot = value

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "LocalVolModel":
        """Calibration not yet implemented for LocalVol."""
        raise NotImplementedError("LocalVol calibration not yet implemented")

    def call_pde(self, maturity: float, strike: float, num_steps: int = 50) -> float:
        """Helper for call option price."""
        return self._cpp_model.call_pde(maturity, strike, num_steps)

    def put_pde(self, maturity: float, strike: float, num_steps: int = 50) -> float:
        """Helper for put option price."""
        return self._cpp_model.put_pde(maturity, strike, num_steps)

    def price_european(
        self,
        maturity: float,
        strike: float,
        option_type: str = "Call",
        n_steps: int = 100,
    ) -> float:
        """
        Price a European option using PDE.
        Matches facade wrapper style.
        """
        if option_type.lower() == "call":
            return self._cpp_model.call_pde(maturity, strike, n_steps)
        elif option_type.lower() == "put":
            return self._cpp_model.put_pde(maturity, strike, n_steps)
        else:
            raise ValueError(f"Unknown option type: {option_type}")

    def price_barrier(
        self,
        maturity: float,
        upper_barrier: float,
        lower_barrier: float,
        n_steps: int = 100,
    ) -> float:
        """
        Price a Double-No-Touch (DNT) option using PDE.
        """
        return self._cpp_model.dnt_pde(maturity, upper_barrier, lower_barrier, n_steps)

    def density(self, maturity: float, num_steps: int = 50) -> list[float]:
        """Compute risk-neutral density."""
        return self._cpp_model.density(maturity, num_steps)

    def get_density(self, maturity: float, n_points: int = 100) -> list[float]:
        """Alias for density to match facade wrapper."""
        return self.density(maturity, n_points)

    def export_local_vol_surface(
        self, times: list[float] | np.ndarray
    ) -> list[list[float]]:
        """
        Export the local volatility surface for the given time points.
        Returns a matrix (list of lists) of local volatilities.
        """
        flat_times = self._to_flat_list(times)
        return self._cpp_model.export_lv(flat_times)

    def simulate(
        self,
        times: list[float] | np.ndarray,
        rands: list[float] | np.ndarray | None = None,
    ) -> list[float]:
        """
        Simulate a path using the local volatility model.
        """
        flat_times = self._to_flat_list(times)
        if rands is not None:
            flat_rands = self._to_flat_list(rands)
            return self._cpp_model.simulate(flat_times, flat_rands)
        else:
            return self._cpp_model.simulate_paths(flat_times)

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "LocalVolModel",
            "spot": self._spot,
            "sabr_models": [
                s.to_dict() if hasattr(s, "to_dict") else {}
                for s in (self._sabr_models or [])
            ],
        }
