from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field

from .base import Instrument


@dataclass
class Portfolio(Instrument):
    """
    A collection of instruments.
    """

    instruments: list[Instrument] = field(default_factory=list)

    def add(self, instrument: Instrument):
        self.instruments.append(instrument)

    def __iter__(self) -> Iterator[Instrument]:
        return iter(self.instruments)

    def __len__(self) -> int:
        return len(self.instruments)

    def group_by_type(self) -> dict[type[Instrument], list[Instrument]]:
        """
        Groups instruments by their concrete class type.
        Useful for batch vectorization.
        """
        groups = defaultdict(list)
        for inst in self.instruments:
            groups[type(inst)].append(inst)
        return groups
