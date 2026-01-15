from typing import Any, Type, TypeVar

T = TypeVar("T")


class Market:
    """
    A container for heterogeneous market data.
    Stores data keyed by (Name, Type).
    """

    def __init__(self):
        # Key: (name, type), Value: data object
        self._data: dict[tuple[str, type], Any] = {}

    def add(self, name: str, data: Any):
        """
        Add a market data object.
        The type key is inferred from the data object's class.
        """
        key = (name, type(data))
        self._data[key] = data

    def get(self, name: str, kind: Type[T]) -> T:
        """
        Retrieve market data by name and expected type.
        Raises KeyError if not found.
        """
        key = (name, kind)
        if key not in self._data:
            raise KeyError(f"Market data not found: {name} of type {kind.__name__}")
        return self._data[key]

    def __contains__(self, item: tuple[str, type]) -> bool:
        return item in self._data

    def to_dict(self) -> dict:
        """
        Serialize market contents.
        Keys are serialized as strings "Name:Type".
        Values must support to_dict().
        """
        serialized = {}
        for (name, kind), data in self._data.items():
            key_str = f"{name}:{kind.__name__}"
            if hasattr(data, "to_dict"):
                serialized[key_str] = data.to_dict()
            else:
                # Fallback or error? For now fallback to simple data types if possible
                serialized[key_str] = data
        return serialized

    @classmethod
    def from_dict(cls, data: dict) -> "Market":
        market = cls()
        # Need a registry of supported types to deserialize correctly
        # For now, let's hardcode support for DiscountCurve as we discover it
        from .curves import DiscountCurve

        type_map = {"DiscountCurve": DiscountCurve}

        for key_str, value_data in data.items():
            try:
                name, type_name = key_str.split(":")
            except ValueError:
                continue  # Skip malformed keys

            if type_name in type_map:
                obj_class = type_map[type_name]
                if hasattr(obj_class, "from_dict"):
                    obj = obj_class.from_dict(value_data)
                    market.add(name, obj)
            else:
                # If we don't know the type, maybe we can't reconstruct perfectly
                # Or we just store the raw dict?
                pass

        return market
