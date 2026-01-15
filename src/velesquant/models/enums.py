from enum import Enum, auto


class BumpType(Enum):
    """Type of bump for numerical sensitivity calculations."""

    ABSOLUTE = auto()  # bump_size is added/subtracted directly
    RELATIVE = auto()  # bump_size is multiplied by the current value


class ModelParam(str, Enum):
    """
    Base class for model parameter enums.
    Inherits from str to allow direct use as attribute names.
    """


class HullWhiteParam(ModelParam):
    """Parameters for Hull-White model."""

    KAPPA = "kappa"  # Mean reversion speed
    SIGMA = "sigma"  # Short rate volatility


class SabrParam(ModelParam):
    """Parameters for SABR model."""

    ALPHA = "alpha"  # Initial volatility
    BETA = "beta"  # CEV exponent
    NU = "nu"  # Vol of vol
    RHO = "rho"  # Correlation
    SHIFT = "shift"  # Shift for negative rates


class HestonParam(ModelParam):
    """Parameters for Heston model."""

    V0 = "v0"  # Initial variance
    KAPPA = "kappa"  # Mean reversion speed
    THETA = "theta"  # Long-run variance
    SIGMA = "sigma"  # Vol of vol
    RHO = "rho"  # Correlation
