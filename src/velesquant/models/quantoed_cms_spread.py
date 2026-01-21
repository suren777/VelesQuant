import numpy as np
from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class QuantoedCMSSpreadModel(Model):
    """
    Quantoed CMS Spread Pricing Model.

    Prices options on the spread between two CMS rates where one or both
    are quantoed into a domestic currency.
    """

    def __init__(
        self,
        expiry1: float,
        tenor1: float,
        fwd1: float,
        annuity1: float,
        pay1: float,
        disc1: float,
        cor_fx1: float,
        atm_vol_fx1: float,
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
        cor_fx2: float,
        atm_vol_fx2: float,
        beta2: float,
        strikes2: np.ndarray,
        quotes2: np.ndarray,
        type2: str,
        corr: float,
    ):
        """
        Initialize Quantoed CMS Spread model.
        """
        self.params = {
            "expiry1": expiry1,
            "tenor1": tenor1,
            "fwd1": fwd1,
            "annuity1": annuity1,
            "pay1": pay1,
            "disc1": disc1,
            "cor_fx1": cor_fx1,
            "atm_vol_fx1": atm_vol_fx1,
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
            "cor_fx2": cor_fx2,
            "atm_vol_fx2": atm_vol_fx2,
            "beta2": beta2,
            "strikes2": strikes2,
            "quotes2": quotes2,
            "type2": type2,
            "corr": corr,
        }

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

        self._cpp_model = native.QuantoedCmsSpread(
            expiry1,
            tenor1,
            fwd1,
            annuity1,
            pay1,
            disc1,
            cor_fx1,
            atm_vol_fx1,
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
            cor_fx2,
            atm_vol_fx2,
            beta2,
            strikes2,
            quotes2,
            ctype2,
            corr,
        )

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "QuantoedCMSSpreadModel":
        """Calibration not yet implemented at high level."""
        raise NotImplementedError(
            "QuantoedCMSSpreadModel calibration not yet implemented"
        )

    def simulate(self, cr1: float, cr2: float) -> list[float]:
        """
        Run Monte Carlo simulation.

        Args:
            cr1: Correlation parameter 1.
            cr2: Correlation parameter 2.
        """
        return self._cpp_model.simulate(cr1, cr2)

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError(
            "General price method not implemented for QuantoedCMSSpreadModel"
        )

    def to_dict(self) -> dict:
        """Serialize model state."""
        d = {"type": "QuantoedCMSSpreadModel"}
        d.update(self.params)
        for key in ["strikes1", "quotes1", "strikes2", "quotes2"]:
            if isinstance(d[key], np.ndarray):
                d[key] = d[key].tolist()
        return d
