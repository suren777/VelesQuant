from dataclasses import dataclass

from .base import Instrument


@dataclass
class Bond(Instrument):
    """
    Abstract Base Bond.
    """

    maturity: float


@dataclass
class ZeroCouponBond(Bond):
    """
    Zero Coupon Bond (pays 1.0 at maturity).
    """

    pass


@dataclass
class CouponBond(Bond):
    """
    Fixed Coupon Bond.
    """

    pay_frequency: float  # 0.5 = Semiannual
    coupon_rate: float  # 0.03 = 3%
    face_value: float = 1.0
