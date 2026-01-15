import numpy as np
from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class CMSSpreadModel(Model):
    """
    CMS Spread Pricing Model.

    Prices options on the spread between two CMS rates.
    Uses SABR dynamics for each underlying swap rate.
    """

    def __init__(
        self,
        expiry1: float,
        tenor1: float,
        fwd1: float,
        annuity1: float,
        pay1: float,
        disc1: float,
        beta1: float,
        strikes1: np.ndarray,
        quotes1: np.ndarray,
        type1: str,
        expiry2: float,
        tenor2: float,
        fwd2: float,
        annuity2: float,
        pay2: float,
        disc2: float,
        beta2: float,
        strikes2: np.ndarray,
        quotes2: np.ndarray,
        type2: str,
        corr: float,
    ):
        """
        Initialize CMS Spread model.
        """
        self.params = {
            "expiry1": expiry1,
            "tenor1": tenor1,
            "fwd1": fwd1,
            "annuity1": annuity1,
            "pay1": pay1,
            "disc1": disc1,
            "beta1": beta1,
            "strikes1": strikes1,
            "quotes1": quotes1,
            "type1": type1,
            "expiry2": expiry2,
            "tenor2": tenor2,
            "fwd2": fwd2,
            "annuity2": annuity2,
            "pay2": pay2,
            "disc2": disc2,
            "beta2": beta2,
            "strikes2": strikes2,
            "quotes2": quotes2,
            "type2": type2,
            "corr": corr,
        }

        # Convert strings to native enums
        ctype1 = (
            native.CalibrationTarget.Price
            if type1.lower() == "price"
            else native.CalibrationTarget.Volatility
        )
        ctype2 = (
            native.CalibrationTarget.Price
            if type2.lower() == "price"
            else native.CalibrationTarget.Volatility
        )

        self._cpp_model = native.CmsSpread(
            expiry1,
            tenor1,
            fwd1,
            annuity1,
            pay1,
            disc1,
            beta1,
            strikes1,
            quotes1,
            ctype1,
            expiry2,
            tenor2,
            fwd2,
            annuity2,
            pay2,
            disc2,
            beta2,
            strikes2,
            quotes2,
            ctype2,
            corr,
        )

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "CMSSpreadModel":
        """Calibration not yet implemented at high level."""
        raise NotImplementedError("CMSSpreadModel calibration not yet implemented")

    def spread_option(self, strike: float, a: float = 1.0, b: float = -1.0) -> float:
        """
        Price a spread option on a*S1 + b*S2.

        Args:
            strike: Strike of the spread option.
            a: Weight for first CMS rate.
            b: Weight for second CMS rate.
        """
        return self._cpp_model.spread_option(strike, a, b)

    def simulate(self) -> list[float]:
        """Run Monte Carlo simulation for CMS rates."""
        return self._cpp_model.simulate()

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError(
            "General price method not implemented for CMSSpreadModel"
        )

    def to_dict(self) -> dict:
        """Serialize model state."""
        d = {"type": "CMSSpreadModel"}
        d.update(self.params)
        for key in ["strikes1", "quotes1", "strikes2", "quotes2"]:
            if isinstance(d[key], np.ndarray):
                d[key] = d[key].tolist()
        return d
