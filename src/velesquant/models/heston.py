from typing import List, Union

import numpy as np

from velesquant import native

from .base import Model


class HestonModel(Model):
    """
    Heston Stochastic Volatility Model Wrapper.
    """

    def __init__(
        self,
        spot: float,
        var0: float,
        kappa: float,
        theta: float,
        xi: float,
        rho: float,
        seed: int = 42,
    ):
        self.spot = spot
        self.var0 = var0
        self.kappa = kappa
        self.theta = theta
        self.xi = xi
        self.rho = rho
        self.seed = seed

        self._cpp_model = native.Heston(spot, var0, kappa, theta, xi, rho, seed)

    def price_option(
        self, maturity: float, forward: float, strike: float, option_type: str = "call"
    ) -> float:
        """
        Price a European option using Heston semi-analytical formula.
        """
        return self._cpp_model.hestonPrice(
            maturity, forward, strike, option_type.lower()
        )

    def simulate(self, times: List[float], forwards: List[float]) -> List[float]:
        """
        Simulate Heston paths.
        """
        return self._cpp_model.simulationHeston(times, forwards)

    def calibrate(
        self,
        maturities: Union[List[float], np.ndarray],
        forwards: Union[List[float], np.ndarray],
        strikes: Union[List[float], np.ndarray],
        quotes: Union[List[float], np.ndarray],
        calibration_target: str = "Price",
    ) -> "HestonModel":
        """
        Calibrate Heston model to market quotes.

        Note: The C++ binding expects matrices. We will reshape inputs to vectors (Nx1 matrices).
        """
        # Convert to numpy arrays and reshape to (-1, 1) for MatrixXd compatibility
        m_mat = np.array(maturities).reshape(-1, 1)
        f_mat = np.array(forwards).reshape(-1, 1)
        s_mat = np.array(strikes).reshape(-1, 1)
        q_mat = np.array(quotes).reshape(-1, 1)

        # C++ expects std::string for quote type
        params = self._cpp_model.calibrate(
            m_mat, f_mat, s_mat, q_mat, calibration_target
        )

        # Update python attributes
        self.var0 = params["var0"]
        self.kappa = params["kappa"]
        self.theta = params["theta"]
        self.xi = params["xi"]
        self.rho = params["rho"]

        return self
