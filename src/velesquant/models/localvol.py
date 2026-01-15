from typing import List

from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model
from .sabr import SabrModel


class LocalVolModel(Model):
    """
    Local Volatility Model Wrapper.

    Constructs a local volatility surface from a list of SABR models
    (one per maturity slice).
    """

    def __init__(self, sabr_models: List[SabrModel], spot: float = 100.0):
        """
        Initialize LocalVolModel from SABR model slices.

        Args:
            sabr_models: List of SabrModel instances (one per maturity).
            spot: Initial spot price.
        """
        self._spot = spot
        self._sabr_models = sabr_models

        # Extract native SABR objects for C++ construction
        native_sabrs = [s._cpp_model for s in sabr_models]
        self._cpp_model = native.LocalVol(native_sabrs)
        self._cpp_model.spot = spot

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
        """
        Price a European call option using PDE method.

        Args:
            maturity: Time to expiry in years.
            strike: Strike price.
            num_steps: Number of time steps for PDE solver.

        Returns:
            Call option price.
        """
        return self._cpp_model.callPDE(maturity, strike, num_steps)

    def put_pde(self, maturity: float, strike: float, num_steps: int = 50) -> float:
        """
        Price a European put option using PDE method.

        Args:
            maturity: Time to expiry in years.
            strike: Strike price.
            num_steps: Number of time steps for PDE solver.

        Returns:
            Put option price.
        """
        return self._cpp_model.putPDE(maturity, strike, num_steps)

    def density(self, maturity: float, num_steps: int = 50) -> List[float]:
        """
        Compute risk-neutral density at maturity.

        Args:
            maturity: Time to maturity.
            num_steps: Number of steps for density computation.

        Returns:
            List of density values.
        """
        return self._cpp_model.density(maturity, num_steps)

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "LocalVolModel",
            "spot": self._spot,
            "sabr_models": [
                s.to_dict() if hasattr(s, "to_dict") else {} for s in self._sabr_models
            ],
        }
