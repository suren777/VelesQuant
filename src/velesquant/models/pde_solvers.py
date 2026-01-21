from velesquant import native

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


def _to_native_swaps(swap_quotes: list[dict]) -> list[native.DefSwap]:
    native_swaps = []
    for q in swap_quotes:
        ds = native.DefSwap()
        ds.expiry = q.get("expiry", 0.0)
        ds.tenor = q.get("tenor", 0.0)
        ds.frequency = q.get("frequency", 0.5)
        ds.swap_rate = q.get("rate", 0.0)
        ds.vol_atm = q.get("vol", 0.0)
        native_swaps.append(ds)
    return native_swaps


class HWPDEModel(Model):
    """
    Hull-White PDE Pricing Model.
    """

    def __init__(
        self,
        kappa: float,
        time_sigmas: list[float],
        sigmas: list[float],
        discount_factor_times: list[float] | None = None,
        discount_factors: list[float] | None = None,
        initial_rate: float | None = None,
        time_thetas: list[float] | None = None,
        thetas: list[float] | None = None,
        discount_curve: list[tuple[float, float]] | None = None,
        grid_points: int = 512,
        time_step: float = 0.001,
    ):
        """
        Initialize HWPDE model.
        Supports two constructor variants:
        1. (kappa, time_sigmas, sigmas, discount_factor_times, discount_factors)
        2. (initial_rate, kappa, time_sigmas, sigmas, time_thetas, thetas)

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

        if initial_rate is not None or time_thetas is not None or thetas is not None:
            # Maintain backward compatibility if possible, or warn/error
            # But since bindings don't support it, we must error or assume params map differently.
            # However, looking at ShortRate1F, it DOES support R0/Thetas. HullWhiteModel binding clearly takes DFs.
            # We will raise error for unsupported args to avoid runtime crash.
            raise ValueError(
                "HWPDEModel currently only supports initialization with discount_factor_times and discount_factors (or discount_curve). "
                "initial_rate/thetas are not supported by the underlying native model."
            )
        else:
            # HullWhiteModel(kappa, timeSigmas, sigmas, discount_factor_times, discount_factors)
            self._model_instance = native.HullWhiteModel(
                kappa, time_sigmas, sigmas, discount_factor_times, discount_factors
            )
            self._cpp_model = native.HWPDE(self._model_instance, grid_points, time_step)

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

    def price_zero_bond(self, maturity: float) -> float:
        """Price a Zero Coupon Bond."""
        return self._cpp_model.price_zero_bond(maturity)

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
        self,
        time_dfs: list[float],
        dfs: list[float],
        swap_quotes: list[dict],
        optimizer_params: dict[str, float] | None = None,
    ) -> "HWPDEModel":
        """
        Calibrate model.
        swap_quotes is a list of dicts with keys: expiry, tenor, frequency, rate, vol.
        """
        native_swaps = _to_native_swaps(swap_quotes)
        self._cpp_model.calibrate(time_dfs, dfs, native_swaps, optimizer_params or {})
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
        time_thetas: list[float] | None = None,
        thetas: list[float] | None = None,
        grid_points: int = 256,
        time_step: float = 0.0025,
    ):
        """
        Initialize ShortRate1FPDEModel.

        Parameters
        ----------
        initial_rate : float
            Initial short rate r0.
        kappa : float
            Mean reversion speed.
        alpha, beta, gamma : float
            Parameters for local volatility/drift specification (e.g. for extended models).
        time_sigmas : list[float]
            Times for time-dependent volatility.
        sigmas : list[float]
            Volatility values corresponding to time_sigmas.
        time_thetas, thetas : list[float], optional
            Time-dependent drift parameters (ignored by current ShortRate1FModel constructor but kept for interface consistency).
        grid_points : int, optional
            Number of spatial grid points for the PDE solver (default 256).
        time_step : float, optional
            Time step for the PDE solver (default 0.0025).
        """
        self.params = locals()
        self.params.pop("self")
        self._model_instance = native.ShortRate1FModel(
            initial_rate,
            kappa,
            alpha,
            beta,
            gamma,
            time_sigmas,
            sigmas,
        )
        self._cpp_model = native.ShortRate1FPDE(
            self._model_instance, grid_points, time_step
        )

    def price_swaption(
        self, expiry: float, tenor: float, strike: float, pay_freq: float = 0.5
    ) -> float:
        return self._cpp_model.price_swaption(expiry, tenor, strike, pay_freq)

    def price_zero_bond(self, maturity: float) -> float:
        """Price a Zero Coupon Bond."""
        return self._cpp_model.price_zero_bond(maturity)

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
        self,
        time_dfs: list[float],
        dfs: list[float],
        swap_quotes: list[dict],
        optimizer_params: dict[str, float] | None = None,
    ) -> "ShortRate1FPDEModel":
        """
        Calibrate model parameters.

        Parameters
        ----------
        time_dfs : list[float]
            Times for discount factors.
        dfs : list[float]
            Discount factors corresponding to time_dfs.
        swap_quotes : list[dict]
            List of dictionaries representing swap calibration instruments.
            Each dict should usually contain: 'expiry', 'tenor', 'rate', 'vol'.
        optimizer_params : dict[str, float], optional
            Parameters for the optimizer (e.g., 'ftol', 'xtol', 'maxfev').

        Returns
        -------
        ShortRate1FPDEModel
            Self, for method chaining.
        """
        native_swaps = _to_native_swaps(swap_quotes)
        self._cpp_model.calibrate(time_dfs, dfs, native_swaps, optimizer_params or {})
        return self

    def to_dict(self) -> dict:
        d = {"type": "ShortRate1FPDEModel"}
        d.update(self.params)
        return d


class ShortRate2FPDEModel(Model):
    """
    Two-Factor Short Rate PDE Model (G2++).
    """

    def __init__(
        self,
        kappa1: float,
        kappa2: float,
        lam: float,
        time_sigma1s: list[float],
        sigma1s: list[float],
        time_sigma2s: list[float],
        sigma2s: list[float],
        time_alphas: list[float],
        alphas: list[float],
        time_step: float = 0.020833333333,
    ):
        """
        Initialize ShortRate2FPDEModel (G2++).

        Parameters
        ----------
        kappa1, kappa2 : float
            Mean reversion speeds for the two factors.
        lam : float
            Correlation (lambda) between the two factors usually denoted as rho, but named lam here.
        time_sigma1s, sigma1s : list[float]
            Time-dependent volatility structure for factor 1.
        time_sigma2s, sigma2s : list[float]
            Time-dependent volatility structure for factor 2.
        time_alphas, alphas : list[float]
            Time-dependent shift/drift parameters.
        time_step : float, optional
            Time step for the 2D PDE solver (default ~1/48 or 0.020833).
        """
        self.params = locals()
        self.params.pop("self")
        self._model_instance = native.ShortRate2FModel(
            kappa1,
            kappa2,
            lam,
            time_sigma1s,
            sigma1s,
            time_sigma2s,
            sigma2s,
            time_alphas,
            alphas,
        )
        self._cpp_model = native.ShortRate2FPDE(self._model_instance, time_step)

    def price_swaption(
        self, expiry: float, tenor: float, strike: float, pay_freq: float = 0.5
    ) -> float:
        return self._cpp_model.price_swaption(expiry, tenor, strike, pay_freq)

    def price_zero_bond(self, maturity: float) -> float:
        return self._cpp_model.price_zero_bond(maturity)

    def calibrate(
        self,
        time_dfs: list[float],
        dfs: list[float],
        swap_quotes: list[dict],
        optimizer_params: dict[str, float] | None = None,
    ) -> "ShortRate2FPDEModel":
        """Calibrate model parameters."""
        native_swaps = _to_native_swaps(swap_quotes)
        self._cpp_model.calibrate(time_dfs, dfs, native_swaps, optimizer_params or {})
        return self

    def to_dict(self) -> dict:
        d = {"type": "ShortRate2FPDEModel"}
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
        """
        Initialize SABR PDE Model.

        Parameters
        ----------
        variant : str
            "Af" for AfSabr or "Antonov" for AntonovSabr (Approximations/PDE schemes).
        alpha : float
            Initial volatility.
        beta : float
            CEV exponent (0 <= beta <= 1).
        nu : float
            Vol-of-vol.
        rho : float
            Correlation between spot and volatility.
        maturity : float
            Time to expiry.
        f : float
            Forward rate/price.
        size_x : int, optional
            Grid size in spatial dimension (default 100).
        size_t : int, optional
            Grid size in time dimension (default 100).
        nd : float, optional
            Number of standard deviations for grid boundaries (default 3.0 or similar).
        shift : float, optional
            Shift parameter for Shifted SABR (default 0.0).
        """
        self.params = locals()
        self.params.pop("self")
        if variant.lower() == "af":
            # Sabr(maturity, forward, beta, alpha, nu, rho, shift)
            self._model_instance = native.Sabr(maturity, f, beta, alpha, nu, rho, shift)
            self._cpp_model = native.AfSabr(self._model_instance, size_x, size_t, nd)
        else:
            self._model_instance = native.Sabr(maturity, f, beta, alpha, nu, rho, shift)
            self._cpp_model = native.AntonovSabr(
                self._model_instance, size_x, size_t, nd
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
