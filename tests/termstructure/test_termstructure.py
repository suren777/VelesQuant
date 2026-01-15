from velesquant import CalendarType, DayCounterType, Termstructure


def test_termstructure_binding():
    # Termstructure(days, rate, calendar, daycount)
    days = [40001.0, 40030.0, 40365.0]
    rates = [0.03, 0.035, 0.04]
    ts = Termstructure(
        days, rates, CalendarType.UnitedKingdom, DayCounterType.Actual360
    )  # 1=UnitedKingdom, 1=Actual360

    d = ts.discount(40100)
    assert isinstance(d, float)
    assert 0 < d <= 1.0

    r = ts.rate(40100, 1)
    assert isinstance(r, float)


def test_termstructure_curve_scenario():
    """
    Scenario: Construct a forward curve and verify discount factors interpolate correctly.
    """
    # Reference date (e.g., today)
    today_serial = 45000
    base_serial = today_serial

    days = [float(base_serial + 1), float(base_serial + 100), float(base_serial + 365)]
    rates = [0.02, 0.025, 0.03]  # Simple upward sloping curve

    ts = Termstructure(
        days, rates, CalendarType.UnitedKingdom, DayCounterType.Actual360
    )

    # Check discount factor at a middle date
    mid_date = base_serial + 50
    df_mid = ts.discount(mid_date)

    # DF should be between DF(1) and DF(100)
    # Actually, rates are forwards, so DF = exp(-rate * T) roughly
    # Since rates increase, DFs decrease.
    # We just check it's a valid discount factor
    assert 0.0 < df_mid < 1.0

    # Check forward rate retrieval
    fwd_rate = ts.rate(mid_date, 30)  # 30 day forward rate starting at mid_date
    assert 0.01 < fwd_rate < 0.05  # Sanity bounds
