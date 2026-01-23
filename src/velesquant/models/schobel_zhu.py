from velesquant import CalibrationTarget, SchobelZhu

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class SchobelZhuModel(Model):
    """
    Schöbel-Zhu Stochastic Volatility Model.

    Similar to Heston but models mean-reversion in volatility (not variance).
    """

    def __init__(
        self,
        spot: float,
        var0: float,
        kappa: float,
        theta: float,
        xi: float,
        rho: float,
    ):
        """
        Initialize Schöbel-Zhu model.
        """
        self._spot = spot
        self._var0 = var0
        self._kappa = kappa
        self._theta = theta
        self._xi = xi
        self._rho = rho

        self._cpp_model = SchobelZhu(spot, var0, kappa, theta, xi, rho)

    @property
    def spot(self) -> float:
        return self._spot

    @spot.setter
    def spot(self, value: float):
        self._spot = value
        # spot is not exposed in C++ model, we just keep it in the wrapper

    @property
    def var0(self) -> float:
        return self._var0

    @var0.setter
    def var0(self, value: float):
        self._var0 = value
        self._cpp_model.var0 = value

    @property
    def kappa(self) -> float:
        return self._kappa

    @kappa.setter
    def kappa(self, value: float):
        self._kappa = value
        self._cpp_model.kappa = value

    @property
    def theta(self) -> float:
        return self._theta

    @theta.setter
    def theta(self, value: float):
        self._theta = value
        self._cpp_model.theta = value

    @property
    def xi(self) -> float:
        return self._xi

    @xi.setter
    def xi(self, value: float):
        self._xi = value
        self._cpp_model.xi = value

    @property
    def rho(self) -> float:
        return self._rho

    @rho.setter
    def rho(self, value: float):
        self._rho = value
        self._cpp_model.rho = value

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "SchobelZhuModel":
        """Calibration via Instrument list not yet implemented."""
        raise NotImplementedError(
            "SchobelZhuModel calibration via Instrument list not yet implemented"
        )

    def calibrate_to_quotes(
        self,
        maturities: list[float],
        forwards: list[float],
        strikes: list[float],
        quotes: list[float],
        target: str = "Price",
    ) -> "SchobelZhuModel":
        """
        Calibrate Schöbel-Zhu model to raw market quotes.
        """
        ctarget = (
            CalibrationTarget.Price
            if target.lower() == "price"
            else CalibrationTarget.Volatility
        )
        self._cpp_model.calibrate(maturities, forwards, strikes, quotes, ctarget)

        # Update attributes from calibrated C++ model
        self._var0 = self._cpp_model.var0
        self._kappa = self._cpp_model.kappa
        self._theta = self._cpp_model.theta
        self._xi = self._cpp_model.xi
        self._rho = self._cpp_model.rho
        return self

    def price_european(self, maturity: float, forward: float, strike: float) -> float:
        """Price a European option."""
        return self._cpp_model.price(maturity, forward, strike)

    def simulate(self, times: list[float], forwards: list[float]) -> list[float]:
        """Run Monte Carlo simulation."""
        return self._cpp_model.simulate(times, forwards)

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError(
            "General price method not implemented for SchobelZhuModel"
        )

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "SchobelZhuModel",
            "spot": self._spot,
            "var0": self._var0,
            "kappa": self._kappa,
            "theta": self._theta,
            "xi": self._xi,
            "rho": self._rho,
        }
