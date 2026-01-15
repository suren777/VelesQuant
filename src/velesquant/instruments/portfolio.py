from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterator, List, Type


from .base import Instrument


@dataclass
class Portfolio(Instrument):
    """
    A collection of instruments.
    """

    instruments: List[Instrument] = field(default_factory=list)

    def add(self, instrument: Instrument):
        self.instruments.append(instrument)

    def __iter__(self) -> Iterator[Instrument]:
        return iter(self.instruments)

    def __len__(self) -> int:
        return len(self.instruments)

    def group_by_type(self) -> dict[Type[Instrument], List[Instrument]]:
        """
        Groups instruments by their concrete class type.
        Useful for batch vectorization.
        """
        groups = defaultdict(list)
        for inst in self.instruments:
            groups[type(inst)].append(inst)
        return groups
