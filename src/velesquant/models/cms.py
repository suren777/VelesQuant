from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class CMSModel(Model):
    """
    Constant Maturity Swap (CMS) Pricing Model.

    Prices CMS instruments with convexity adjustment using SABR dynamics.
    """

    def __init__(
        self,
        expiry: float,
        tenor: float,
        forward: float,
        annuity: float,
        pay_cms: float,
        discount_cms: float,
        beta: float = 0.85,
        alpha: float = 0.5,
        nu: float = 0.25,
        rho: float = -0.75,
    ):
        """
        Initialize CMS model.

        Args:
            expiry: CMS rate fixing expiry in years.
            tenor: Underlying swap rate tenor in years.
            forward: Forward swap rate.
            annuity: Annuity (PV01) of the swap rate.
            pay_cms: Payment time for CMS cashflow.
            discount_cms: Discount factor to CMS payment date.
            beta: SABR beta parameter.
            alpha: SABR alpha parameter.
            nu: SABR nu parameter.
            rho: SABR rho parameter.
        """
        self.expiry = expiry
        self.tenor = tenor
        self.forward = forward
        self.annuity = annuity
        self.pay_cms = pay_cms
        self.discount_cms = discount_cms
        self.beta = beta
        self.alpha = alpha
        self.nu = nu
        self.rho = rho

        self._cpp_model = native.CMS(
            expiry, tenor, forward, annuity, pay_cms, discount_cms, beta, alpha, nu, rho
        )

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "CMSModel":
        """Calibration not yet implemented."""
        raise NotImplementedError("CMSModel calibration not yet implemented")

    @property
    def adjusted_forward(self) -> float:
        """Get convexity-adjusted forward CMS rate."""
        return self._cpp_model.getForward()

    @property
    def maturity(self) -> float:
        """Get CMS maturity."""
        return self._cpp_model.getMaturity()

    def fair_value(self, strike: float, option_type: str = "Call") -> float:
        """
        Price a CMS caplet/floorlet.

        Args:
            strike: Strike rate.
            option_type: "Call" (caplet) or "Put" (floorlet).

        Returns:
            Option fair value.
        """
        otype = native.OptionType.Call
        if option_type.lower() == "put":
            otype = native.OptionType.Put
        return self._cpp_model.fairValue(strike, otype)

    def implied_vol(self, strike: float) -> float:
        """
        Get SABR implied volatility at strike.

        Args:
            strike: Strike rate.

        Returns:
            Implied Black volatility.
        """
        return self._cpp_model.getImpliedVol(strike)

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "CMSModel",
            "expiry": self.expiry,
            "tenor": self.tenor,
            "forward": self.forward,
            "annuity": self.annuity,
            "pay_cms": self.pay_cms,
            "discount_cms": self.discount_cms,
            "beta": self.beta,
            "alpha": self.alpha,
            "nu": self.nu,
            "rho": self.rho,
        }
