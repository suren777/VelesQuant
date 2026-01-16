from velesquant.models import HybridHWModel


def test_hybrid_hw_wrapper_init():
    # s0, v0, r0, kappa, eta, rho, sigma1, sigma2, a
    model = HybridHWModel(100.0, 0.04, 0.02, 2.0, 0.3, -0.5, 0.2, 0.01, 0.1)
    assert model._cpp_model is not None
    assert model.fair_value(1.0, 100.0) > 0

    d = model.to_dict()
    assert d["type"] == "HybridHWModel"
    assert d["s0"] == 100.0
    assert d["rho"] == -0.5
