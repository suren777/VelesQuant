import velesquant.native as m


def test_log_basket_init_and_sim():
    # 2 Assets
    spots = [100.0, 100.0]
    strikes = [100.0, 100.0]
    maturities = [1.0, 1.0]

    # 2x1 matrix provided as list of lists? No, it's vector<vector<double>>
    # If M assets, looks like Forwards might be matrix of [asset][time]?
    # Or just inputs... Header says Vdoub Spot, Vdoub Strike, Vdoub Maturities, Mdoub Forwards, Mdoub IV, Mdoub correlation
    # Let's assume simplest case: 1 time step?

    # Looking at usage, Forwards/IV might be term structures.
    # Let's try simple inputs.

    forwards = [[100.0], [100.0]]
    iv = [[0.2], [0.2]]
    correlation = [[1.0, 0.5], [0.5, 1.0]]

    lb = m.LogNormalBasket(spots, strikes, maturities, forwards, iv, correlation)

    assert lb.get_nassets() == 2

    schedule = [1.0]
    paths = lb.simulate_basket(schedule)
    # Returns Mdoub: [sim][asset]? Or [asset][sim]?
    # Assuming list of lists
    assert isinstance(paths, list)
    assert len(paths) > 0
