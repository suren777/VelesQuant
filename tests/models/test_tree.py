"""Tests for CTree (Binomial/Trinomial Tree) model."""

from velesquant import CTree, ExerciseStyle, OptionType, TreeType


def test_ctree_binding():
    """Test that CTree can be instantiated and called."""
    # Market data: spot, times, forwards, implied vols
    S = 100.0
    T = [0.5, 1.0]
    F = [100.5, 101.0]
    IV = [0.2, 0.22]
    r = [0.01, 0.01]
    q = [0.002, 0.002]

    tree = CTree(S, T, F, IV, r, q)
    # European call using binomial tree
    price = tree.calculateBinomial(
        100.0, 1.0, 100, ExerciseStyle.European, OptionType.Call, TreeType.Recombining
    )
    assert isinstance(price, float)
    assert price >= 0.0


def test_ctree_american_put():
    """Test American put pricing - should be >= European put."""
    S = 100.0
    T = [0.5, 1.0]
    F = [100.5, 101.0]
    IV = [0.2, 0.22]

    tree = CTree(S, T, F, IV)

    euro_put = tree.calculateBinomial(
        100.0, 1.0, 100, ExerciseStyle.European, OptionType.Put, TreeType.Recombining
    )
    amer_put = tree.calculateBinomial(
        100.0, 1.0, 100, ExerciseStyle.American, OptionType.Put, TreeType.Recombining
    )

    # American put should be >= European put
    assert amer_put >= euro_put


def test_ctree_trinomial():
    """Test trinomial tree pricing."""
    S = 100.0
    T = [0.5, 1.0]
    F = [100.5, 101.0]
    IV = [0.2, 0.22]

    tree = CTree(S, T, F, IV)
    price = tree.calculateTrinomial(
        100.0, 1.0, 50, ExerciseStyle.European, OptionType.Call, TreeType.Recombining
    )
    assert isinstance(price, float)
    assert price >= 0.0
