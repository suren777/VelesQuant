from abc import ABC
from dataclasses import dataclass


@dataclass
class Instrument(ABC):
    """
    Abstract Base Class for all financial instruments.
    Instruments are pure data containers that describe *what* is traded.
    """

    pass
