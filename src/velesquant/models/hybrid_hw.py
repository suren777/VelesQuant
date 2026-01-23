from velesquant import HHW

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class HybridHWModel(Model):
    """
    Hybrid Hull-White Model.

    Combines an equity process with stochastic volatility and a
    stochastic interest rate process (Hull-White).
    """

    def __init__(
        self,
        s0: float,
        v0: float,
        r0: float,
        kappa: float,
        eta: float,
        rho: float,
        sigma1: float,
        sigma2: float,
        a: float,
    ):
        """
        Initialize Hybrid Hull-White model.
        """
        self._s0 = s0
        self._v0 = v0
        self._r0 = r0
        self._kappa = kappa
        self._eta = eta
        self._rho = rho
        self._sigma1 = sigma1
        self._sigma2 = sigma2
        self._a = a

        self._cpp_model = HHW(s0, v0, r0, kappa, eta, rho, sigma1, sigma2, a)

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "HybridHWModel":
        """Calibration not yet implemented."""
        raise NotImplementedError("HybridHWModel calibration not yet implemented")

    def fair_value(self, maturity: float, strike: float) -> float:
        """Price a European option."""
        # HHWPrice is a method in HHW class
        return self._cpp_model.price(maturity, strike)

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError(
            "General price method not implemented for HybridHWModel"
        )

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "HybridHWModel",
            "s0": self._s0,
            "v0": self._v0,
            "r0": self._r0,
            "kappa": self._kappa,
            "eta": self._eta,
            "rho": self._rho,
            "sigma1": self._sigma1,
            "sigma2": self._sigma2,
            "a": self._a,
        }
