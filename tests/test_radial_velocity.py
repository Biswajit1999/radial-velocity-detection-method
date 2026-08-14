"""Executable checks on the Kepler solver and RV physics, including a
regression guard against the earlier 1/sqrt(Mstar) scaling-law bug
(the correct limiting scaling is Mstar^-2/3, not Mstar^-1/2).
"""

import numpy as np
import radial_velocity_demo as rv


def test_kepler_equation_residual_is_tiny():
    mean_anomaly = np.linspace(0, 2 * np.pi, 50)
    for ecc in (0.0, 0.3, 0.6, 0.9):
        ecc_anomaly = rv.solve_kepler(mean_anomaly.copy(), ecc)
        residual = ecc_anomaly - ecc * np.sin(ecc_anomaly) - mean_anomaly
        assert np.max(np.abs(residual)) < 1e-10, f"ecc={ecc}: max residual {np.max(np.abs(residual))}"


def test_circular_orbit_reduces_to_sinusoid():
    # For e=0 and omega=0, true anomaly = mean anomaly, so the Keplerian
    # RV curve must reduce exactly to a plain cosine.
    period, t0, k, gamma, omega = 10.0, 0.0, 50.0, 0.0, 0.0
    time = np.linspace(0, 20, 200)
    rv_vals = rv.keplerian_rv(time, period, t0, k, 0.0, omega, gamma)
    expected = k * np.cos(2 * np.pi * (time - t0) / period)
    assert np.allclose(rv_vals, expected, atol=1e-8)


def test_minimum_mass_round_trip():
    # Independently re-derive K from the recovered minimum mass (using a
    # forward calculation written directly in this test, not the module's
    # internal residual function) and confirm it reproduces the injected K.
    cases = [(4.23, 92.0, 0.04, 1.02), (10.0, 5.0, 0.0, 0.5), (365.25, 0.09, 0.0167, 1.0)]
    for period, k, ecc, star_mass in cases:
        mp_earth = rv.minimum_mass_mearth(period, k, ecc, star_mass)
        period_s = period * rv.DAY_S
        mp_kg = mp_earth * 5.972e24
        star_kg = star_mass * rv.M_SUN
        predicted_k = (2 * np.pi * rv.G / period_s) ** (1 / 3) * mp_kg / (star_kg + mp_kg) ** (2 / 3) / np.sqrt(1 - ecc**2)
        assert abs(predicted_k - k) / k < 1e-6, f"{(period, k, ecc, star_mass)}: predicted K {predicted_k}"


def test_mass_scales_as_mstar_two_thirds_not_one_half():
    # Regression guard: an earlier version of this repo's README claimed
    # K scales as 1/sqrt(Mstar); the correct limiting relation (Mp << Mstar)
    # is Mstar^-2/3. At fixed K, minimum mass should scale as Mstar^(2/3),
    # clearly distinguishable from the incorrect Mstar^(1/2) scaling.
    k, period, ecc = 5.0, 10.0, 0.0
    mp_low = rv.minimum_mass_mearth(period, k, ecc, 0.5)
    mp_high = rv.minimum_mass_mearth(period, k, ecc, 2.0)
    ratio = mp_high / mp_low

    expected_two_thirds = (2.0 / 0.5) ** (2 / 3)
    expected_one_half_bug = (2.0 / 0.5) ** 0.5

    assert abs(ratio - expected_two_thirds) / expected_two_thirds < 0.01
    assert abs(ratio - expected_one_half_bug) / expected_one_half_bug > 0.1


def test_51_pegasi_b_worked_example_is_self_consistent():
    # The real archive solution for 51 Peg b: P=4.2307969 d, Mp sini=0.464 MJ,
    # e=0.0042, Mstar=1.069 Msun, K=55.73 m/s (Rosenthal et al., NASA
    # Exoplanet Archive). Since the archive mass was itself derived from
    # this same K, this is an internal-consistency check, not independent
    # validation -- exactly as the README states.
    period, ecc, star_mass, published_k = 4.2307969, 0.0042, 1.069, 55.73
    mp_earth_published = 0.464 * 317.8  # Jupiter masses -> Earth masses
    period_s = period * rv.DAY_S
    mp_kg = mp_earth_published * 5.972e24
    star_kg = star_mass * rv.M_SUN
    predicted_k = (2 * np.pi * rv.G / period_s) ** (1 / 3) * mp_kg / (star_kg + mp_kg) ** (2 / 3) / np.sqrt(1 - ecc**2)
    assert abs(predicted_k - published_k) / published_k < 0.01
