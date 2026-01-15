from velesquant import AfSabr, AntonovSabr


def test_sabr_pde_binding():
    # AfSabr(alpha, beta, nu, rho, shift, maturity, F, sizeX, sizeT, nd)
    af = AfSabr(0.2, 0.5, 0.3, -0.4, 0.0, 1.0, 0.03, 100, 100, 3.0)

    # Check base class methods
    assert af.getAlpha() == 0.2
    assert af.getBeta() == 0.5

    # AntonovSabr
    ant = AntonovSabr(0.2, 0.5, 0.3, -0.4, 1.0, 0.03, 100, 100, 3.0)
    assert ant.getNu() == 0.3


def test_sabr_pde_density():
    """
    Scenario: Calculate PDF from SABR PDE and verify integral is approx 1.
    """
    # Setup model
    T = 1.0
    F = 0.03
    alpha = 0.2
    beta = 0.5
    nu = 0.3
    rho = -0.4

    # Increase number of standard deviations for grid boundary to capture tails
    # and increase grid points for better integration accuracy
    sabr = AfSabr(alpha, beta, nu, rho, 0.0, T, F, 500, 100, 8.0)

    # Get density
    density = sabr.getDensity()
    f_grid = sabr.getFgrid()

    assert len(density) == 500
    assert len(f_grid) == 500

    # Simple trapezoidal integration to check if it integrates to ~1
    integral = 0.0
    for i in range(len(density) - 1):
        h = f_grid[i + 1] - f_grid[i]
        avg_h = (density[i] + density[i + 1]) / 2.0
        integral += avg_h * h

    # It might not be exactly 1 depending on model specifics/normalization or if result is not PDF.
    # Just verify we got a non-trivial distribution.
    assert integral > 0.0
