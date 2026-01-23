
import numpy as np

from velesquant import CalibrationTarget, Heston, OptionType

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

        self._cpp_model = Heston(spot, var0, kappa, theta, xi, rho, seed)

    def price_option(
        self, maturity: float, forward: float, strike: float, option_type: str = "call"
    ) -> float:
        """
        Price a European option using Heston semi-analytical formula.
        """
        native_type = (
            OptionType.Call if option_type.lower() == "call" else OptionType.Put
        )
        return self._cpp_model.price(maturity, forward, strike, native_type)

    def simulate(self, times: list[float], forwards: list[float]) -> list[float]:
        """
        Simulate Heston paths.
        """
        return self._cpp_model.simulate(times, forwards)

    def calibrate(
        self,
        maturities: list[float] | np.ndarray,
        forwards: list[float] | np.ndarray,
        strikes: list[float] | np.ndarray,
        quotes: list[float] | np.ndarray,
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

        # Map string to Enum
        native_target = (
            CalibrationTarget.Price
            if calibration_target.capitalize() == "Price"
            else CalibrationTarget.Volatility
        )
        params = self._cpp_model.calibrate(m_mat, f_mat, s_mat, q_mat, native_target)

        # Update python attributes
        self.var0 = params["var0"]
        self.kappa = params["kappa"]
        self.theta = params["theta"]
        self.xi = params["xi"]
        self.rho = params["rho"]

        return self
