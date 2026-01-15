from typing import Protocol, Type, TypeVar, runtime_checkable

T = TypeVar("T")


@runtime_checkable
class MarketData(Protocol):
    """
    Protocol for market data containers.
    Any object that can be used as market data in pricing should implement this.
    """

    def to_dict(self) -> dict:
        """Serialize to dictionary."""


@runtime_checkable
class MarketDataContainer(Protocol):
    """
    Protocol for containers that hold multiple market data objects (e.g., Market).
    """

    def get(self, name: str, kind: Type[T]) -> T:
        """Retrieve a specific piece of market data by name and type."""

    def add(self, name: str, data: MarketData) -> None:
        """Add market data to the container."""
