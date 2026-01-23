from velesquant import OptionType, Swaption

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class SabrSwaptionModel(Model):
    """
    SABR-based Swaption Pricing Model.

    Prices swaptions using SABR volatility dynamics.
    This wraps the native.Swaption class which handles
    convexity-adjusted forward rates and SABR smile.
    """

    def __init__(
        self,
        expiry: float,
        tenor: float,
        forward: float,
        annuity: float,
        beta: float = 0.85,
        alpha: float = 0.5,
        nu: float = 0.25,
        rho: float = -0.75,
    ):
        """
        Initialize SABR Swaption model.

        Args:
            expiry: Option expiry in years.
            tenor: Underlying swap tenor in years.
            forward: Forward swap rate.
            annuity: Annuity (PV01) of the swap.
            beta: SABR beta parameter (CEV exponent).
            alpha: SABR alpha parameter (initial vol).
            nu: SABR nu parameter (vol of vol).
            rho: SABR rho parameter (correlation).
        """
        self.expiry = expiry
        self.tenor = tenor
        self.forward = forward
        self.annuity = annuity
        self.beta = beta
        self.alpha = alpha
        self.nu = nu
        self.rho = rho

        self._cpp_model = Swaption(
            expiry, tenor, forward, annuity, beta, alpha, nu, rho
        )

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "SabrSwaptionModel":
        """Calibration not yet implemented."""
        raise NotImplementedError("SabrSwaptionModel calibration not yet implemented")

    def fair_value(self, strike: float, option_type: str = "Call") -> float:
        """
        Price the swaption at a given strike.

        Args:
            strike: Fixed rate strike.
            option_type: "Call" (payer) or "Put" (receiver).

        Returns:
            Swaption fair value.
        """
        otype = OptionType.Call
        if option_type.lower() == "put":
            otype = OptionType.Put
        return self._cpp_model.fair_value(strike, otype)

    def implied_vol(self, strike: float) -> float:
        """
        Get SABR implied volatility at strike.

        Args:
            strike: Strike rate.

        Returns:
            Implied Black volatility.
        """
        return self._cpp_model.get_implied_vol(strike)

    def swap_value(self, strike: float) -> float:
        """
        Value of the underlying swap (forward - strike) * annuity.

        Args:
            strike: Fixed rate strike.

        Returns:
            Swap fair value.
        """
        return self._cpp_model.swap_fair_value(strike)

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "SabrSwaptionModel",
            "expiry": self.expiry,
            "tenor": self.tenor,
            "forward": self.forward,
            "annuity": self.annuity,
            "beta": self.beta,
            "alpha": self.alpha,
            "nu": self.nu,
            "rho": self.rho,
        }
