from typing import List, Union

import numpy as np

# Import native module
from velesquant import native

from .base import Model


class SabrModel(Model):
    """
    SABR Stochastic Volatility Model Wrapper.
    """

    def __init__(
        self,
        maturity: float,
        forward: float,
        beta: float = 0.5,
        alpha: float = 0.2,
        nu: float = 0.5,
        rho: float = -0.5,
        shift: float = 0.0,
    ):
        self.maturity = maturity
        self.forward = forward
        self.alpha = alpha
        self.beta = beta
        self.nu = nu
        self.rho = rho
        self.shift = shift

        # Construct the underlying C++ model immediately
        if shift == 0.0:
            self._cpp_model = native.Sabr(maturity, forward, beta, alpha, nu, rho)
        else:
            self._cpp_model = native.Sabr(
                maturity, forward, beta, alpha, nu, rho, shift
            )

    def implied_vol(self, strike: float) -> float:
        """
        Calculate implied volatility for a given strike.
        """
        return self._cpp_model.impliedVol(strike)

    def calibrate(
        self,
        strikes: Union[List[float], np.ndarray],
        quotes: Union[List[float], np.ndarray],
        calibration_target: str = "Volatility",
    ) -> "SabrModel":
        """
        Calibrate the SABR model to a smile.

        Args:
            strikes: List or array of strike prices
            quotes: List or array of market quotes (Vols or Prices)
            calibration_target: "Volatility" or "Price"

        Returns:
            Computed (new) SabrModel with calibrated parameters.
        """
        # Convert inputs to Eigen-compatible numpy arrays (MatrixXd)
        # The binding expects slices or 1D arrays, ensuring they are 2D (Nx1) for MatrixXd conversion might be safest
        # logic in binding: loop rows/cols.

        s_mat = np.array(strikes).reshape(-1, 1)
        q_mat = np.array(quotes).reshape(-1, 1)

        target = native.CalibrationTarget.Volatility
        if calibration_target.lower() == "price":
            target = native.CalibrationTarget.Price

        # Call C++ calibration
        # The C++ method returns a dict of parameters and updates internal state
        params = self._cpp_model.calibrate(s_mat, q_mat, target)

        # Create a new instance with calibrated params
        # Note: The C++ object 'self._cpp_model' is already updated in place!
        # But to be 'immutable-ish' Pythonic, we might return a new wrapper?
        # For now, let's update self and return self.

        self.alpha = params["alpha"]
        self.nu = params["nu"]
        self.rho = params["rho"]

        return self
