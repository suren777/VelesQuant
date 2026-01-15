from typing import Any, List, Union

import numpy as np

# We import the C++ native module.
# Depending on installation, it might be velesquant.native
from velesquant import native

from ..instruments.bonds import Bond, CouponBond, ZeroCouponBond
from ..instruments.portfolio import Portfolio
from ..instruments.rates import Swaption
from ..market.container import Market
from ..market.curves import DiscountCurve
from .base import Model


class HullWhiteModel(Model):
    """
    Hull-White 1-Factor Model Wrapper.
    """

    def __init__(self, kappa: float, sigma: float):
        self.kappa = kappa
        self.sigma = sigma
        self._cpp_model = None

    def calibrate(self, instruments: list[Any], market_data: Any) -> "HullWhiteModel":
        # TODO: Implement calibration logic calling native.HullWhite.calibrate
        return self

    def to_dict(self) -> dict:
        return {"type": "HullWhiteModel", "kappa": self.kappa, "sigma": self.sigma}

    @classmethod
    def from_dict(cls, data: dict) -> "HullWhiteModel":
        return cls(kappa=data["kappa"], sigma=data["sigma"])

    def simulate(self, times: List[float], curve: DiscountCurve) -> List[float]:
        """
        Simulate Hull-White paths.
        Requires a curve to construct the time-dependent drift.
        """
        # We need to construct the C++ object temporarily
        # Use a dummy max_time covering all simulation times?
        max_time = max(times) + 1.0
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]
        hw = native.HullWhite(self.kappa, time_sigmas, sigmas, curve.times, curve.dfs)
        return hw.simulation(times)

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
        max_time = max(expiry, maturity) + 1.0
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]
        hw = native.HullWhite(self.kappa, time_sigmas, sigmas, curve.times, curve.dfs)

        # Map string to enum
        otype = native.OptionType.Call
        if option_type.lower() == "put":
            otype = native.OptionType.Put

        return hw.optionBond(expiry, maturity, strike, otype)

    def price(
        self,
        instrument: Union[Swaption, Bond, Portfolio],
        market_data: Union[DiscountCurve, Market],
        curve_name: str = "USD",
    ) -> Union[float, np.ndarray, dict]:
        """
        Price an instrument using the Hull-White model.
        Supports vectorized pricing if instrument attributes are arrays.
        Supports Portfolio pricing by batching instruments.
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

        is_vectorized = False
        if hasattr(instrument, "expiry"):
            if isinstance(instrument.expiry, (list, np.ndarray)):
                is_vectorized = True
        elif hasattr(instrument, "maturity"):
            if isinstance(instrument.maturity, (list, np.ndarray)):
                is_vectorized = True

        if is_vectorized:
            # Use numpy vectorize to map the scalar pricing functions
            if isinstance(instrument, Swaption):
                # We need to vectorize over the fields of the instrument
                # Create a vectorized function that takes scalar components
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
                # For bonds, depends on type
                if isinstance(instrument, ZeroCouponBond):
                    vfunc = np.vectorize(lambda m: self._price_zc_scalar(m, curve))
                    return vfunc(instrument.maturity)
                elif isinstance(instrument, CouponBond):
                    vfunc = np.vectorize(
                        lambda m, pf, cr, fv: self._price_cb_scalar(
                            m, pf, cr, fv, curve
                        )
                    )
                    return vfunc(
                        instrument.maturity,
                        instrument.pay_frequency,
                        instrument.coupon_rate,
                        instrument.face_value,
                    )

        # Scalar Fallback
        if isinstance(instrument, Swaption):
            return self._price_swaption(instrument, curve)
        elif isinstance(instrument, Bond):
            return self._price_bond(instrument, curve)
        else:
            raise TypeError(f"Unsupported instrument type: {type(instrument)}")

    def _price_portfolio(self, portfolio: Portfolio, curve: DiscountCurve) -> dict:
        """
        Smart batching: Group instruments, vectorize them, and price.
        Returns a dict mapping {InstrumentType: np.ndarray of prices}
        or a flat list matching input order?
        For now, let's return a list of prices matching the input order.
        To do that, we need to be careful with re-ordering.

        Alternative: Return total NPV? Usually Portfolios want total risk/NPV.
        Let's return Total NPV for now as it's the simplest "pricing" of a portfolio.
        Future: Return detailed breakdown.
        """
        total_npv = 0.0
        groups = portfolio.group_by_type()

        for inst_type, instruments in groups.items():
            if not instruments:
                continue

            # Attempt to stack
            if inst_type == Swaption:
                # Stack
                expiries = np.array([i.expiry for i in instruments])
                tenors = np.array([i.tenor for i in instruments])
                strikes = np.array([i.strike for i in instruments])
                # Using default pay_freq if None
                freqs = np.array([i.pay_frequency or 0.5 for i in instruments])

                # Create vectorized instrument
                vec_inst = Swaption(
                    expiry=expiries,
                    tenor=tenors,
                    strike=strikes,
                    pay_frequency=freqs,
                )

                # Price
                prices = self.price(vec_inst, curve)
                total_npv += np.sum(prices)

            elif inst_type == ZeroCouponBond:
                maturities = np.array([i.maturity for i in instruments])
                vec_inst = ZeroCouponBond(maturity=maturities)
                prices = self.price(vec_inst, curve)
                total_npv += np.sum(prices)

            # TODO: Handle CouponBond stacking (tricky due to differing payment schedules?
            # Vectorization of CouponBond assumes SAME structure or broadcastable.
            # Different maturities/coupons might not broadcast well if they don't align in a grid.
            # Actually our scalar-vectorization simply iterates over the scalars.
            # So CouponBond vectorization works for arrays of scalars (all bonds differ).
            # So yes, we can stack CouponBonds too!)
            elif inst_type == CouponBond:
                maturities = np.array([i.maturity for i in instruments])
                freqs = np.array([i.pay_frequency for i in instruments])
                rates = np.array([i.coupon_rate for i in instruments])
                faces = np.array([i.face_value for i in instruments])

                vec_inst = CouponBond(
                    maturity=maturities,
                    pay_frequency=freqs,
                    coupon_rate=rates,
                    face_value=faces,
                )
                prices = self.price(vec_inst, curve)
                total_npv += np.sum(prices)

            else:
                # Fallback to loop
                for inst in instruments:
                    total_npv += self.price(inst, curve)

        return total_npv

    # Refactor existing _price methods to split setup vs calc if needed,
    # but for now we can rely on helper methods that take scalars.

    def _price_swaption_scalar(self, expiry, tenor, strike, pay_freq, curve):
        # Helper to call C++ with scalars
        max_time = expiry + tenor + 1.0
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]
        hw = native.HullWhite(self.kappa, time_sigmas, sigmas, curve.times, curve.dfs)
        return hw.swaption(expiry, tenor, strike, pay_freq)

    def _price_swaption(self, instrument: Swaption, curve: DiscountCurve) -> float:
        # Re-use scalar method
        return self._price_swaption_scalar(
            instrument.expiry,
            instrument.tenor,
            instrument.strike,
            instrument.pay_frequency or 0.5,
            curve,
        )  # Assuming default freq logic if missing? actually default is usually handled in dataclass or explicitly passed.

    def _price_bond(self, instrument: Bond, curve: DiscountCurve) -> float:
        # Construct C++ model for ZC pricing
        max_time = instrument.maturity + 1.0
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]

        hw = native.HullWhite(self.kappa, time_sigmas, sigmas, curve.times, curve.dfs)

        if isinstance(instrument, ZeroCouponBond):
            return hw.ZC(instrument.maturity)

        elif isinstance(instrument, CouponBond):
            value = 0.0
            t = instrument.pay_frequency
            while t <= instrument.maturity + 1e-9:
                df = hw.ZC(t)
                coupon = (
                    instrument.face_value
                    * instrument.coupon_rate
                    * instrument.pay_frequency
                )
                value += coupon * df
                t += instrument.pay_frequency

            value += instrument.face_value * hw.ZC(instrument.maturity)
            return value

        else:
            raise TypeError(f"Unsupported bond type: {type(instrument)}")

    def _price_zc_scalar(self, maturity, curve):
        max_time = maturity + 1.0
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]
        hw = native.HullWhite(self.kappa, time_sigmas, sigmas, curve.times, curve.dfs)
        return hw.ZC(maturity)

    def _price_cb_scalar(self, maturity, pay_freq, coupon_rate, face_value, curve):
        max_time = maturity + 1.0
        time_sigmas = [0.0, max_time]
        sigmas = [self.sigma, self.sigma]
        hw = native.HullWhite(self.kappa, time_sigmas, sigmas, curve.times, curve.dfs)

        value = 0.0
        t = pay_freq
        while t <= maturity + 1e-9:
            df = hw.ZC(t)
            coupon = face_value * coupon_rate * pay_freq
            value += coupon * df
            t += pay_freq

        value += face_value * hw.ZC(maturity)
        return value
