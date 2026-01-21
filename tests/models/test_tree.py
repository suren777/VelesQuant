"""Tests for CTree (Binomial/Trinomial Tree) model."""

from velesquant import native
from velesquant.models import TreeModel


def test_ctree_binding():
    """Test that CTree can be instantiated and called."""
    # Market data: spot, times, forwards, implied vols
    S = 100.0
    T = [0.5, 1.0]
    F = [100.5, 101.0]
    IV = [0.2, 0.22]
    r = [0.01, 0.01]
    q = [0.002, 0.002]

    tree = native.CTree(S, T, F, IV, r, q)
    # European call using binomial tree
    price = tree.calculate_binomial(
        100.0,
        1.0,
        100,
        native.ExerciseStyle.European,
        native.OptionType.Call,
        native.TreeType.Recombining,
    )
    assert isinstance(price, float)
    assert price >= 0.0


def test_ctree_american_put():
    """Test American put pricing - should be >= European put."""
    S = 100.0
    T = [0.5, 1.0]
    F = [100.5, 101.0]
    IV = [0.2, 0.22]

    tree = native.CTree(S, T, F, IV)

    euro_put = tree.calculate_binomial(
        100.0,
        1.0,
        100,
        native.ExerciseStyle.European,
        native.OptionType.Put,
        native.TreeType.Recombining,
    )
    amer_put = tree.calculate_binomial(
        100.0,
        1.0,
        100,
        native.ExerciseStyle.American,
        native.OptionType.Put,
        native.TreeType.Recombining,
    )

    # American put should be >= European put
    assert amer_put >= euro_put


def test_ctree_trinomial():
    """Test trinomial tree pricing."""
    S = 100.0
    T = [0.5, 1.0]
    F = [100.5, 101.0]
    IV = [0.2, 0.22]

    tree = native.CTree(S, T, F, IV)
    price = tree.calculate_trinomial(
        100.0,
        1.0,
        50,
        native.ExerciseStyle.European,
        native.OptionType.Call,
        native.TreeType.Recombining,
    )
    assert isinstance(price, float)
    assert price >= 0.0


def test_trees_wrapper_init():
    model = TreeModel(100.0, [0.5, 1.0], [102.0, 105.0], [0.2, 0.2])
    price = model.price_binomial(105.0, 1.0, n_nodes=50)
    assert price > 0

    price_tri = model.price_trinomial(105.0, 1.0, n_nodes=50)
    assert price_tri > 0

    d = model.to_dict()
    assert d["type"] == "TreeModel"
    assert d["spot"] == 100.0
