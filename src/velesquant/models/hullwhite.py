from typing import Any

import numpy as np

from velesquant import (
    HullWhiteAnalyticEngine,
    OptionType,
    implied_vol,
)
from velesquant import (
    HullWhiteModel as NativeHullWhiteModel,
)

from ..instruments.bonds import Bond, CouponBond, ZeroCouponBond
from ..instruments.portfolio import Portfolio
from ..instruments.rates import Swaption
from ..market.container import Market
from ..market.curves import DiscountCurve
from .base import Model


class HullWhiteModel(Model):
    """
    Hull-White 1-Factor Model Wrapper.
    Uses native C++ HullWhiteModel and HullWhiteAnalyticEngine.
    """

    def __init__(self, kappa: float, sigma: float):
        self.kappa = kappa
        self.sigma = sigma
        self._cpp_model = None
        self._cached_native_objects = None

    def calibrate(self, instruments: list[Any], market_data: Any) -> "HullWhiteModel":
        # TODO: Implement calibration logic
        return self

    def to_dict(self) -> dict:
        return {"type": "HullWhiteModel", "kappa": self.kappa, "sigma": self.sigma}

    @classmethod
    def from_dict(cls, data: dict) -> "HullWhiteModel":
        return cls(kappa=data["kappa"], sigma=data["sigma"])

    def _create_native_objects(self, curve: DiscountCurve):
        """Helper to create native Model"""
        # Check cache validity
        if self._cached_native_objects:
            model, engine, cached_curve_id, cached_kappa, cached_sigma = (
                self._cached_native_objects
            )
            # We use id(curve) as a simple check, assuming curve doesn't mutate in place without id change
            # or simply check logic.
            # Ideally we check parameters.
            if (
                id(curve) == cached_curve_id
                and self.kappa == cached_kappa
                and self.sigma == cached_sigma
            ):
                return model, engine

        max_time = curve.times[-1] + 10.0  # Ensure coverage
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]

        # Core Model
        model = NativeHullWhiteModel(
            self.kappa, time_sigmas, sigmas, curve.times, curve.dfs
        )
        # Analytic Engine
        # pylint: disable=no-member
        engine = HullWhiteAnalyticEngine(model)  # type: ignore

        # Update cache
        self._cached_native_objects = (
            model,
            engine,
            id(curve),
            self.kappa,
            self.sigma,
        )
        return model, engine

    def simulate(self, times: list[float], curve: DiscountCurve) -> list[float]:
        """
        Simulate Hull-White paths.
        """
        model, _ = self._create_native_objects(curve)
        return model.simulate(times)

    def price_bond_option(
        self,
        expiry: float,
        maturity: float,
        strike: float,
        curve: DiscountCurve,
        option_type: str = "Call",
    ) -> float:
        """
        Price a European option on a Zero Coupon Bond.
        """
        _, engine = self._create_native_objects(curve)

        # Map string to enum
        otype = OptionType.Call
        if option_type.lower() == "put":
            otype = OptionType.Put

        return engine.option_bond(expiry, maturity, strike, otype)

    def get_swap_rate(
        self,
        expiry: float,
        tenor: float,
        curve: DiscountCurve,
        pay_frequency: float = 0.5,
    ) -> float:
        """
        Calculate Par Swap Rate.
        S = (P(Start) - P(End)) / Annuity
        """
        model, _ = self._create_native_objects(curve)

        start_time = expiry
        end_time = expiry + tenor

        # pylint: disable=no-member
        p_start = model.get_discount_factor(start_time)  # type: ignore
        p_end = model.get_discount_factor(end_time)

        # Annuity
        annuity = 0.0
        n_periods = int(round(tenor / pay_frequency))

        for i in range(1, n_periods + 1):
            t = start_time + i * pay_frequency
            annuity += pay_frequency * model.get_discount_factor(t)

        if annuity < 1e-12:
            return 0.0

        return (p_start - p_end) / annuity

    def swaption_implied_vol(
        self,
        expiry: float,
        tenor: float,
        price: float,
        curve: DiscountCurve,
        pay_frequency: float = 0.5,
    ) -> float:
        """
        Calculate implied volatility for a swaption.
        """
        # Get Swap Rate (Forward)
        swap_rate = self.get_swap_rate(expiry, tenor, curve, pay_frequency)

        # Calculate Annuity again (could optimize via caching or helper return)
        model, _ = self._create_native_objects(curve)
        start_time = expiry
        n_periods = int(round(tenor / pay_frequency))
        annuity = 0.0
        for i in range(1, n_periods + 1):
            t = start_time + i * pay_frequency
            # pylint: disable=no-member
            annuity += pay_frequency * model.get_discount_factor(t)  # type: ignore

        if annuity < 1e-12:
            return 0.0

        normalized_price = price / annuity

        # Assuming ATM Strike = Forward Swap Rate
        try:
            return implied_vol(expiry, swap_rate, swap_rate, normalized_price)
        except Exception:
            return 0.0

    def price(
        self,
        instrument: Swaption | Bond | Portfolio,
        market_data: DiscountCurve | Market,
        curve_name: str = "USD",
    ) -> float | np.ndarray | dict:
        """
        Price an instrument using the Hull-White model.
        """
        # Resolve Market Data
        if isinstance(market_data, Market):
            curve = market_data.get(curve_name, DiscountCurve)
        elif isinstance(market_data, DiscountCurve):
            curve = market_data
        else:
            raise TypeError(f"Unsupported market data type: {type(market_data)}")

        if isinstance(instrument, Portfolio):
            return self._price_portfolio(instrument, curve)

        # Vectorization check
        is_vectorized = False
        if hasattr(instrument, "expiry") and isinstance(
            instrument.expiry, (list, np.ndarray)
        ) or hasattr(instrument, "maturity") and isinstance(
            instrument.maturity, (list, np.ndarray)
        ):
            is_vectorized = True

        if is_vectorized:
            if isinstance(instrument, Swaption):
                vfunc = np.vectorize(
                    lambda e, t, s, p: self._price_swaption_scalar(e, t, s, p, curve)
                )
                return vfunc(
                    instrument.expiry,
                    instrument.tenor,
                    instrument.strike,
                    instrument.pay_frequency,
                )
            elif isinstance(instrument, Bond):
                vfunc = np.vectorize(
                    lambda m: self._price_bond_scalar(m, instrument, curve)
                )
                # For ZCB, we iterate over maturity. CouponBond is complex to vectorize simply by attribute
                # assuming ZCB vectorization by maturity for now as per test
                if hasattr(instrument, "maturity") and isinstance(
                    instrument.maturity, (list, np.ndarray)
                ):
                    return vfunc(instrument.maturity)
                else:
                    # Fallback if vectorization trigger unclear or other attrs
                    pass

        # Scalar Fallback
        if isinstance(instrument, Swaption):
            return self._price_swaption(instrument, curve)
        elif isinstance(instrument, Bond):
            return self._price_bond(instrument, curve)
        else:
            raise TypeError(f"Unsupported instrument type: {type(instrument)}")

    def _price_portfolio(self, portfolio: Portfolio, curve: DiscountCurve) -> float:
        total_npv = 0.0
        for inst in portfolio.instruments:
            # Simplest iteration for now
            total_npv += self.price(inst, curve)
        return total_npv

    def _price_swaption_scalar(self, expiry, tenor, strike, pay_freq, curve):
        _, engine = self._create_native_objects(curve)
        return engine.swaption(expiry, tenor, strike, pay_freq or 0.5)

    def _price_swaption(self, instrument: Swaption, curve: DiscountCurve) -> float:
        return self._price_swaption_scalar(
            instrument.expiry,
            instrument.tenor,
            instrument.strike,
            instrument.pay_frequency,
            curve,
        )

    def _price_bond_scalar(self, maturity, instrument, curve):
        # Helper that treats 'instrument' as a template but overrides maturity
        model, _ = self._create_native_objects(curve)

        if isinstance(instrument, ZeroCouponBond):
            # pylint: disable=no-member
            return model.get_discount_factor(maturity)  # type: ignore
        # CouponBond vectorization not fully supported in this simple scalar helper yet
        return 0.0

    def _price_bond(self, instrument: Bond, curve: DiscountCurve) -> float:
        model, _ = self._create_native_objects(curve)

        if isinstance(instrument, ZeroCouponBond):
            return model.get_discount_factor(instrument.maturity)

        elif isinstance(instrument, CouponBond):
            value = 0.0
            t = instrument.pay_frequency
            while t <= instrument.maturity + 1e-9:
                df = model.get_discount_factor(t)
                coupon = (
                    instrument.face_value
                    * instrument.coupon_rate
                    * instrument.pay_frequency
                )
                value += coupon * df
                t += instrument.pay_frequency
            value += instrument.face_value * model.get_discount_factor(
                instrument.maturity
            )
            return value
        else:
            raise TypeError(f"Unsupported bond type: {type(instrument)}")
