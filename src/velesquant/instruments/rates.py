from dataclasses import dataclass

from .base import Instrument


@dataclass
class Swaption(Instrument):
    """
    A Swaption (Option to enter a Swap).
    """

    expiry: float  # Time to expiry in years
    tenor: float  # Swap length in years
    strike: float  # Fixed rate
    pay_frequency: float = 0.5  # Default semi-annual
    receiver: bool = True  # True = Receiver Swaption (Call on Bond?), False = Payer
