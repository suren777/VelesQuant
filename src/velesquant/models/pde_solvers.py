from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class HWPDEModel(Model):
    """
    Hull-White PDE Pricing Model.
    """

    def __init__(
        self,
        kappa: float,
        time_sigmas: list[float],
        sigmas: list[float],
        discount_factor_times: list[float] = None,
        discount_factors: list[float] = None,
        initial_rate: float = None,
        time_thetas: list[float] = None,
        thetas: list[float] = None,
        discount_curve: list[tuple[float, float]] = None,
    ):
        """
        Initialize HWPDE model.
        Supports two constructor variants:
        1. (kappa, timeSigmas, sigmas, discount_factor_times, discount_factors)
        2. (initial_rate, kappa, timeSigmas, sigmas, timeThetas, thetas)

        Optional `discount_curve` can be passed as a list of (time, discount_factor) tuples
        instead of separate lists.
        """
        self.params = locals()
        self.params.pop("self")

        # Handle discount_curve convenience argument
        if discount_curve is not None:
            if discount_factor_times is not None or discount_factors is not None:
                raise ValueError(
                    "Cannot provide both discount_curve and discount_factor_times/discount_factors"
                )
            discount_factor_times = [pt[0] for pt in discount_curve]
            discount_factors = [pt[1] for pt in discount_curve]

        if initial_rate is not None:
            self._cpp_model = native.HWPDE(
                initial_rate, kappa, time_sigmas, sigmas, time_thetas, thetas
            )
        else:
            self._cpp_model = native.HWPDE(
                kappa, time_sigmas, sigmas, discount_factor_times, discount_factors
            )

    def price_swaption(
        self, expiry: float, tenor: float, strike: float, pay_freq: float = 0.5
    ) -> float:
        """Price a swaption."""
        return self._cpp_model.price_swaption(expiry, tenor, strike, pay_freq)

    def price_bermudan(
        self,
        expiry: float,
        tenor: float,
        exercises: list[float],
        strike: float,
        pay_freq: float = 0.5,
    ) -> float:
        """Price a Bermudan swaption."""
        return self._cpp_model.price_bermudan(
            expiry, tenor, exercises, strike, pay_freq
        )

    def price_callable_swap(
        self,
        expiry: float,
        tenor: float,
        exercises: list[float],
        coupon: float,
        strike: float,
        pay_freq: float = 0.5,
        option_type: str = "Call",
    ) -> float:
        """Price a Callable Swap."""
        native_type = (
            native.OptionType.Call
            if option_type.capitalize() == "Call"
            else native.OptionType.Put
        )
        return self._cpp_model.price_callable_swap(
            expiry, tenor, exercises, coupon, strike, pay_freq, native_type
        )

    def price_zbo(
        self, expiry: float, maturity: float, strike: float, option_type: str = "Call"
    ) -> float:
        """Price a Zero Coupon Bond Option."""
        native_type = (
            native.OptionType.Call
            if option_type.capitalize() == "Call"
            else native.OptionType.Put
        )
        return self._cpp_model.price_zero_bond_option(
            expiry, maturity, strike, native_type
        )

    def price_cbo(
        self,
        expiry: float,
        tenor: float,
        coupon: float,
        strike: float,
        pay_freq: float = 0.5,
        option_type: str = "Call",
    ) -> float:
        """Price a Coupon Bond Option."""
        native_type = (
            native.OptionType.Call
            if option_type.capitalize() == "Call"
            else native.OptionType.Put
        )
        return self._cpp_model.price_coupon_bond_option(
            expiry, tenor, coupon, strike, pay_freq, native_type
        )

    def calibrate(
        self, time_dfs: list[float], dfs: list[float], swap_quotes: list[dict]
    ) -> "HWPDEModel":
        """
        Calibrate model.
        swap_quotes is a list of native.DefSwap-like objects.
        """
        # Convert dicts or simple objects to native.DefSwap
        native_swaps = []
        for q in swap_quotes:
            ds = native.DefSwap()
            ds.expiry = q.get("expiry", 0)
            ds.tenor = q.get("tenor", 0)
            ds.frequency = q.get("frequency", 0.5)
            ds.swap_rate = q.get("rate", 0)
            ds.vol_atm = q.get("vol", 0)
            native_swaps.append(ds)

        self._cpp_model.calibrate(time_dfs, dfs, native_swaps)
        return self

    def simulate(self, time_points: list[float]) -> list[float]:
        """Run Monte Carlo simulation using the PDE grid."""
        return self._cpp_model.simulate(time_points)

    def to_dict(self) -> dict:
        """Serialize model state."""
        d = {"type": "HWPDEModel"}
        d.update(self.params)
        return d


class ShortRate1FPDEModel(Model):
    """
    General 1-Factor Short Rate PDE Model.
    """

    def __init__(
        self,
        initial_rate: float,
        kappa: float,
        alpha: float,
        beta: float,
        gamma: float,
        time_sigmas: list[float],
        sigmas: list[float],
        time_thetas: list[float],
        thetas: list[float],
    ):
        self.params = locals()
        self.params.pop("self")
        self._cpp_model = native.ShortRate1FPDE(
            initial_rate,
            kappa,
            alpha,
            beta,
            gamma,
            time_sigmas,
            sigmas,
            time_thetas,
            thetas,
        )

    def price_swaption(
        self, expiry: float, tenor: float, strike: float, pay_freq: float = 0.5
    ) -> float:
        return self._cpp_model.price_swaption(expiry, tenor, strike, pay_freq)

    def price_zbo(
        self, expiry: float, maturity: float, strike: float, option_type: str = "Call"
    ) -> float:
        """Price a Zero Coupon Bond Option."""
        native_type = (
            native.OptionType.Call
            if option_type.capitalize() == "Call"
            else native.OptionType.Put
        )
        return self._cpp_model.price_zero_bond_option(
            expiry, maturity, strike, native_type
        )

    def price_cbo(
        self,
        expiry: float,
        tenor: float,
        coupon: float,
        strike: float,
        pay_freq: float = 0.5,
        option_type: str = "Call",
    ) -> float:
        """Price a Coupon Bond Option."""
        native_type = (
            native.OptionType.Call
            if option_type.capitalize() == "Call"
            else native.OptionType.Put
        )
        return self._cpp_model.price_coupon_bond_option(
            expiry, tenor, coupon, strike, pay_freq, native_type
        )

    def calibrate(self, instruments, market_data) -> "ShortRate1FPDEModel":
        raise NotImplementedError

    def to_dict(self) -> dict:
        d = {"type": "ShortRate1FPDEModel"}
        d.update(self.params)
        return d


class SabrPDEModel(Model):
    """
    SABR PDE Pricing Model (AfSabr or AntonovSabr).
    """

    def __init__(
        self,
        variant: str,
        alpha: float,
        beta: float,
        nu: float,
        rho: float,
        maturity: float,
        f: float,
        size_x: int = 100,
        size_t: int = 100,
        nd: float = 0.5,
        shift: float = 0.0,
    ):
        self.params = locals()
        self.params.pop("self")
        if variant.lower() == "af":
            self._cpp_model = native.AfSabr(
                alpha, beta, nu, rho, shift, maturity, f, size_x, size_t, nd
            )
        else:
            self._cpp_model = native.AntonovSabr(
                alpha, beta, nu, rho, maturity, f, size_x, size_t, nd
            )

    def get_density(self) -> list[float]:
        """Get risk-neutral density from PDE."""
        return self._cpp_model.get_density()

    def get_f_grid(self) -> list[float]:
        """Get forward price grid."""
        return self._cpp_model.get_f_grid()

    def calibrate(self, instruments, market_data) -> "SabrPDEModel":
        raise NotImplementedError

    def to_dict(self) -> dict:
        d = {"type": "SabrPDEModel"}
        d.update(self.params)
        return d
