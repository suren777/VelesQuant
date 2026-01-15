import copy
from abc import ABC, abstractmethod
from typing import ClassVar, Optional, Type, Union

from ..instruments.base import Instrument
from ..market.base import MarketData, MarketDataContainer
from .enums import BumpType, ModelParam

# Type alias for market data: either a single data object or a container
MarketDataInput = Union[MarketData, MarketDataContainer]


class Model(ABC):
    """
    Abstract Base Class for pricing models.
    Models describe *how* market variables evolve.

    Subclasses should define a `Param` class attribute mapping to their
    specific ModelParam enum for type-safe sensitivity calculations.
    """

    # Subclasses should override with their specific Param enum
    # e.g., Param: ClassVar[Type[ModelParam]] = HullWhiteParam
    Param: ClassVar[Optional[Type[ModelParam]]] = None

    @abstractmethod
    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "Model":
        """
        Calibrate the model to a set of instruments.
        Returns a new calibrated model instance or self.
        """
        pass

    def with_params(self, **kwargs) -> "Model":
        """
        Return a copy of this model with updated parameters.
        Subclasses may override for more efficient copying.
        """
        cloned = copy.copy(self)
        for key, value in kwargs.items():
            if not hasattr(cloned, key):
                raise AttributeError(f"Model has no parameter '{key}'")
            setattr(cloned, key, value)
        return cloned

    # -------------------------------------------------------------------------
    # First-Order Greeks (with numerical defaults)
    # -------------------------------------------------------------------------

    def delta(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        First derivative with respect to underlying price/rate (dV/dS).

        Note: Requires bumping market data (spot/forward), not model params.
        Override in subclass or use market_sensitivity() when available.
        """
        raise NotImplementedError(
            "delta requires bumping market data (spot). Override in subclass."
        )

    def vega(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        First derivative with respect to volatility (dV/dσ).

        Default: Uses sensitivity() on SIGMA/ALPHA param if available.
        """
        vol_param = self._get_param("SIGMA") or self._get_param("ALPHA")
        if vol_param:
            return self.sensitivity(
                instrument, market_data, vol_param, bump_size, BumpType.ABSOLUTE
            )
        raise NotImplementedError(
            f"vega: No SIGMA or ALPHA parameter in {self.__class__.__name__}.Param"
        )

    def theta(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1 / 365,  # 1 day
    ) -> float:
        """
        First derivative with respect to time (dV/dt).

        Note: Requires bumping instrument expiry, not model params.
        Override in subclass.
        """
        raise NotImplementedError(
            "theta requires bumping instrument expiry. Override in subclass."
        )

    def rho(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        First derivative with respect to interest rate (dV/dr).

        Default: Uses sensitivity() on KAPPA param if available (for IR models).
        """
        kappa_param = self._get_param("KAPPA")
        if kappa_param:
            return self.sensitivity(
                instrument, market_data, kappa_param, bump_size, BumpType.ABSOLUTE
            )
        raise NotImplementedError(
            f"rho: No KAPPA parameter in {self.__class__.__name__}.Param"
        )

    # -------------------------------------------------------------------------
    # Second-Order Greeks (with numerical defaults)
    # -------------------------------------------------------------------------

    def gamma(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        Second derivative with respect to underlying price/rate (d²V/dS²).

        Note: Requires bumping market data. Override in subclass.
        """
        raise NotImplementedError(
            "gamma requires bumping market data. Override in subclass."
        )

    def vanna(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        Cross-derivative: sensitivity of delta to volatility (d²V/dSdσ).

        Default: Numerical approximation using delta at bumped vol.
        """
        vol_param = self._get_param("SIGMA") or self._get_param("ALPHA")
        if not vol_param:
            raise NotImplementedError(
                f"vanna: No SIGMA or ALPHA parameter in {self.__class__.__name__}.Param"
            )

        param_name = vol_param.value
        original_vol = getattr(self, param_name)

        model_up = self.with_params(**{param_name: original_vol + bump_size})
        model_down = self.with_params(**{param_name: original_vol - bump_size})

        delta_up = model_up.delta(instrument, market_data, bump_size)
        delta_down = model_down.delta(instrument, market_data, bump_size)

        return (delta_up - delta_down) / (2 * bump_size)

    def volga(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        Second derivative with respect to volatility (d²V/dσ²). Also called vomma.

        Default: Numerical second derivative via three-point stencil.
        """
        vol_param = self._get_param("SIGMA") or self._get_param("ALPHA")
        if not vol_param:
            raise NotImplementedError(
                f"volga: No SIGMA or ALPHA parameter in {self.__class__.__name__}.Param"
            )

        param_name = vol_param.value
        original_vol = getattr(self, param_name)

        model_up = self.with_params(**{param_name: original_vol + bump_size})
        model_mid = self
        model_down = self.with_params(**{param_name: original_vol - bump_size})

        price_up = model_up.price(instrument, market_data)
        price_mid = model_mid.price(instrument, market_data)
        price_down = model_down.price(instrument, market_data)

        return (price_up - 2 * price_mid + price_down) / (bump_size**2)

    def charm(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        bump_size: float = 1e-4,
    ) -> float:
        """
        Cross-derivative: sensitivity of delta to time (d²V/dSdt).

        Note: Requires bumping instrument expiry. Override in subclass.
        """
        raise NotImplementedError(
            "charm requires bumping instrument expiry. Override in subclass."
        )

    # -------------------------------------------------------------------------
    # Numerical Bump-and-Reprice
    # -------------------------------------------------------------------------

    def sensitivity(
        self,
        instrument: Instrument,
        market_data: MarketDataInput,
        param: ModelParam,
        bump_size: float = 1e-4,
        bump_type: BumpType = BumpType.ABSOLUTE,
    ) -> float:
        """
        Generic numerical sensitivity via bump-and-reprice (immutable).

        Creates temporary copies of the model with bumped parameters,
        prices with each, and computes the central difference.

        Args:
            instrument: The instrument to price.
            market_data: Market data container.
            param: The parameter to bump (from model's Param enum).
            bump_size: Size of the bump.
            bump_type: BumpType.ABSOLUTE or BumpType.RELATIVE.

        Returns:
            Numerical first derivative dV/d(param).
        """
        param_name = param.value  # Extract string value from enum

        if not hasattr(self, param_name):
            raise AttributeError(f"Model has no parameter '{param_name}'")

        original_value = getattr(self, param_name)

        if bump_type == BumpType.RELATIVE:
            bump = original_value * bump_size
        else:
            bump = bump_size

        # Create bumped copies (immutable - does not modify self)
        model_up = self.with_params(**{param_name: original_value + bump})
        model_down = self.with_params(**{param_name: original_value - bump})

        price_up = model_up.price(instrument, market_data)
        price_down = model_down.price(instrument, market_data)

        return (price_up - price_down) / (2 * bump)

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """
        Price an instrument. Subclasses should override this.
        Used by the generic sensitivity() method.
        """
        raise NotImplementedError

    # -------------------------------------------------------------------------
    # Helper Methods
    # -------------------------------------------------------------------------

    def _get_param(self, name: str) -> Optional[ModelParam]:
        """Get a parameter from this model's Param enum by name, or None if not found."""
        if self.Param is None:
            return None
        try:
            return self.Param[name]
        except KeyError:
            return None
