from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model
from .sabr import SabrModel


class SkewMCModel(Model):
    """
    Skew-adjusted Monte Carlo Simulation Model.

    Uses a set of SABR models (one per maturity) to simulate asset paths
    consistent with the market volatility skew.
    """

    def __init__(self, sabr_models: list[SabrModel]):
        """
        Initialize Skew MC model.

        Args:
            sabr_models: A list of SabrModel instances, one for each maturity slice.
        """
        self.sabr_models = sabr_models

        # Extract native Sabr objects for the C++ constructor
        native_sabrs = [m._cpp_model for m in sabr_models]
        self._cpp_model = native.SkewMC(native_sabrs)

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "SkewMCModel":
        """Calibration is Typically handled by calibrating the underlying SABR models."""
        return self

    def simulate(self, times: list[float], spot: float, kappa: float) -> list[float]:
        """
        Run Monte Carlo simulation.
        """
        return self._cpp_model.simulate(times, spot, kappa)

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError(
            "General price method not implemented for SkewMCModel"
        )

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "SkewMCModel",
            "sabr_models": [m.to_dict() for m in self.sabr_models],
        }
