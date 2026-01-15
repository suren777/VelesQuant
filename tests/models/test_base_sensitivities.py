from typing import Any

import pytest

from velesquant.models.base import Model


class ConcreteModel(Model):
    def calibrate(self, instruments: list[Any], market_data: Any) -> "Model":
        return self


def test_base_sensitivities_raise_not_implemented():
    model = ConcreteModel()

    # Test that base methods raise NotImplementedError
    with pytest.raises(NotImplementedError):
        model.delta(None, None)

    with pytest.raises(NotImplementedError):
        model.gamma(None, None)

    with pytest.raises(NotImplementedError):
        model.vega(None, None)
