from dataclasses import dataclass

import numpy as np


@dataclass
class DiscountCurve:
    """
    A Yield Term Structure based on Discount Factors.
    """

    times: list[float]
    dfs: list[float]

    def __post_init__(self):
        if len(self.times) != len(self.dfs):
            raise ValueError("Times and DFs must have the same length")

        self.log_dfs = np.log(self.dfs)

    def discount(self, t: float) -> float:
        """
        Return the discount factor at time t.
        Uses log-linear interpolation on discount factors (linear on log-dfs).
        """
        # Linear interpolation on log discount factors
        log_df = np.interp(
            t,
            self.times,
            self.log_dfs,
            left=self.log_dfs[0],  # Extrapolate flat log-DF (dummy)
            right=self.log_dfs[-1],  # Extrapolate flat log-DF (dummy)
        )
        return float(np.exp(log_df))

    def to_dict(self) -> dict:
        return {"type": "DiscountCurve", "times": self.times, "dfs": self.dfs}

    @classmethod
    def from_dict(cls, data: dict) -> "DiscountCurve":
        return cls(times=data["times"], dfs=data["dfs"])
