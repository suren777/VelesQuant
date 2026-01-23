from velesquant import CTree, ExerciseStyle, OptionType, TreeType

from ..instruments.base import Instrument
from .base import MarketDataInput, Model


class TreeModel(Model):
    """
    Binomial and Trinomial Tree Pricing Model.

    Supports American, European, and Bermudan exercise styles.
    """

    def __init__(
        self,
        spot: float,
        times: list[float],
        forwards: list[float],
        ivs: list[float],
        r: list[float] | None = None,
        q: list[float] | None = None,
    ):
        """
        Initialize Tree model.
        """
        self._spot = spot
        self._times = times
        self._forwards = forwards
        self._ivs = ivs
        self._r = r if r is not None else []
        self._q = q if q is not None else []

        self._cpp_model = CTree(
            spot, self._times, self._forwards, self._ivs, self._r, self._q
        )

    def calibrate(
        self, instruments: list[Instrument], market_data: MarketDataInput
    ) -> "TreeModel":
        """Calibration not implemented for general instrument list."""
        return self

    def _get_enum(self, enum_cls, value: str):
        """Helper to get enum member with validation."""
        try:
            return getattr(enum_cls, value)
        except AttributeError:
            members = [m.name for m in enum_cls.__members__.values()]
            raise ValueError(
                f"Invalid {enum_cls.__name__}: '{value}'. Expected one of {members}"
            ) from None

    def price_binomial(
        self,
        strike: float,
        maturity: float,
        n_nodes: int = 100,
        exercise_style: str = "American",
        option_type: str = "Call",
        tree_type: str = "Recombining",
    ) -> float:
        """
        Price an option using a Binomial Tree.
        """
        # Get enums from native
        style = self._get_enum(ExerciseStyle, exercise_style)
        otype = self._get_enum(OptionType, option_type)
        ttype = self._get_enum(TreeType, tree_type)

        return self._cpp_model.calculate_binomial(
            strike, maturity, n_nodes, style, otype, ttype
        )

    def price_trinomial(
        self,
        strike: float,
        maturity: float,
        n_nodes: int = 100,
        exercise_style: str = "American",
        option_type: str = "Call",
        tree_type: str = "Recombining",
    ) -> float:
        """Price an option using a Trinomial Tree."""
        style = self._get_enum(ExerciseStyle, exercise_style)
        otype = self._get_enum(OptionType, option_type)
        ttype = self._get_enum(TreeType, tree_type)

        return self._cpp_model.calculate_trinomial(
            strike, maturity, n_nodes, style, otype, ttype
        )

    def price(self, instrument: Instrument, market_data: MarketDataInput) -> float:
        """General price method (to be implemented)."""
        raise NotImplementedError("General price method not implemented for TreeModel")

    def to_dict(self) -> dict:
        """Serialize model state."""
        return {
            "type": "TreeModel",
            "spot": self._spot,
            "times": self._times,
            "forwards": self._forwards,
            "ivs": self._ivs,
            "r": self._r,
            "q": self._q,
        }
