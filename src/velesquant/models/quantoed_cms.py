import numpy as np
from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class QuantoedCMSModel(Model):
    """
    Quantoed CMS Pricing Model.

    Prices CMS instruments in a domestic currency while the underlying rate
    is in a foreign currency, applying a quanto adjustment based on FX correlation.
    """

    def __init__(
        self,
        expiry: float,
        tenor: float,
        fwd: float,
        annuity: float,
        pay: float,
        disc: float,
        cor_fx: float,
        atm_vol_fx: float,
        beta: float,
        strikes: np.ndarray,
        quotes: np.ndarray,
        calibration_type: str = "Price",
    ):
        """
        Initialize Quantoed CMS model.
        """
        self.params = {
            "expiry": expiry,
            "tenor": tenor,
            "fwd": fwd,
            "annuity": annuity,
            "pay": pay,
            "disc": disc,
            "cor_fx": cor_fx,
            "atm_vol_fx": atm_vol_fx,
            "beta": beta,
            "strikes": strikes,
            "quotes": quotes,
            "calibration_type": calibration_type,
        }

        ctype = (
            native.CalibrationTarget.Price
            if calibration_type.lower() == "price"
            else native.CalibrationTarget.Volatility
        )

        self._cpp_model = native.QuantoedCMS(
            expiry,
            tenor,
            fwd,
            annuity,
            pay,
            disc,
            cor_fx,
            atm_vol_fx,
            beta,
            strikes,
            quotes,
            ctype,
        )

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "QuantoedCMSModel":
        """Calibration not yet implemented at high level."""
        raise NotImplementedError("QuantoedCMSModel calibration not yet implemented")

    def fair_value(self, strike: float, option_type: str = "Call") -> float:
        """
        Price a quantoed CMS caplet/floorlet.
        """
        otype = native.OptionType.Call
        if option_type.lower() == "put":
            otype = native.OptionType.Put
        # Note: In bindings, fairValue return value and args are fixed
        return self._cpp_model.fair_value(strike, otype)

    def get_forward(self) -> float:
        """Get the quanto-adjusted forward rate."""
        return self._cpp_model.get_forward()

    def simulate(self, corr_rn: float) -> list[float]:
        """
        Run Monte Carlo simulation.
        """
        return self._cpp_model.simulate(corr_rn)

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError(
            "General price method not implemented for QuantoedCMSModel"
        )

    def to_dict(self) -> dict:
        """Serialize model state."""
        d = {"type": "QuantoedCMSModel"}
        d.update(self.params)
        if isinstance(d["strikes"], np.ndarray):
            d["strikes"] = d["strikes"].tolist()
        if isinstance(d["quotes"], np.ndarray):
            d["quotes"] = d["quotes"].tolist()
        return d
