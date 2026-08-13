"""Radial velocity method demonstration: inject a realistic Keplerian
stellar wobble signal, sampled the way a real ground-based spectrograph
campaign actually observes (irregular cadence, seasonal gaps), with
real-world-level Doppler noise, then recover the orbital period via a
Lomb-Scargle periodogram and fit the full Keplerian orbit to recover
the velocity semi-amplitude K and the planet's minimum mass.

This is a PEDAGOGICAL DEMONSTRATION with simulated data, not a specific
real target's raw archival spectra (see README.md for why, and see this
portfolio's *-exoplanet-report repos for 11 planets analyzed directly
from real archival JWST/HST/Spitzer/ground-based data). The injected
semi-amplitude, period, and noise level are drawn from real, published
regimes (a hot-Jupiter-class signal and HARPS's own published ~1 m/s
single-measurement precision for a bright, quiet star), so the recovery
statistics below are a genuine, physically grounded test of the
method's real-world sensitivity.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.signal import lombscargle

plt.style.use(["science", "no-latex"])

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

rng = np.random.default_rng(seed=7)

G = 6.674e-11
M_SUN = 1.989e30
DAY_S = 86400.0

# Injected "ground truth" orbit (realistic hot-Jupiter regime).
TRUE_PERIOD_DAYS = 4.23
TRUE_K_MS = 92.0
TRUE_ECC = 0.04
TRUE_OMEGA_RAD = 1.1
TRUE_T0_DAYS = 3.5
GAMMA_MS = 0.0
STAR_MASS_MSUN = 1.02  # real-like Sun-like host

N_OBS = 60
BASELINE_DAYS = 200.0
RV_NOISE_MS = 1.0  # real HARPS-class single-measurement precision for a bright quiet star
JITTER_MS = 1.5  # real-like stellar jitter added in quadrature


def solve_kepler(mean_anomaly: np.ndarray, ecc: float) -> np.ndarray:
    ecc_anomaly = mean_anomaly.copy()
    for _ in range(50):
        ecc_anomaly -= (ecc_anomaly - ecc * np.sin(ecc_anomaly) - mean_anomaly) / (1 - ecc * np.cos(ecc_anomaly))
    return ecc_anomaly


def keplerian_rv(time: np.ndarray, period: float, t0: float, k: float, ecc: float, omega: float, gamma: float) -> np.ndarray:
    mean_anomaly = 2 * np.pi * ((time - t0) / period % 1.0)
    ecc_anomaly = solve_kepler(mean_anomaly, ecc)
    true_anomaly = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(ecc_anomaly / 2), np.sqrt(1 - ecc) * np.cos(ecc_anomaly / 2))
    return k * (np.cos(true_anomaly + omega) + ecc * np.cos(omega)) + gamma


def minimum_mass_mearth(period_days: float, k_ms: float, ecc: float, star_mass_msun: float) -> float:
    period_s = period_days * DAY_S
    star_mass_kg = star_mass_msun * M_SUN

    def residual(mp_kg: float) -> float:
        predicted_k = (2 * np.pi * G / period_s) ** (1 / 3) * mp_kg / (star_mass_kg + mp_kg) ** (2 / 3) / np.sqrt(1 - ecc**2)
        return predicted_k - k_ms

    mp_kg = brentq(residual, 1e22, 1e29)
    return mp_kg / 5.972e24


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    # Irregular real-campaign-like sampling: random nights within observing
    # windows, with seasonal gaps (no data for ~40% of the baseline).
    candidate_days = np.sort(rng.uniform(0, BASELINE_DAYS, N_OBS * 3))
    in_season = (candidate_days % 100) < 60
    time = np.sort(rng.choice(candidate_days[in_season], size=min(N_OBS, in_season.sum()), replace=False))

    true_rv = keplerian_rv(time, TRUE_PERIOD_DAYS, TRUE_T0_DAYS, TRUE_K_MS, TRUE_ECC, TRUE_OMEGA_RAD, GAMMA_MS)
    total_noise = np.sqrt(RV_NOISE_MS**2 + JITTER_MS**2)
    rv = true_rv + rng.normal(0, total_noise, size=time.size)
    rv_err = np.full_like(rv, total_noise)

    # Lomb-Scargle periodogram to find the orbital period.
    freq_grid = np.linspace(2 * np.pi / 20.0, 2 * np.pi / 1.0, 20000)
    power = lombscargle(time, rv - rv.mean(), freq_grid, normalize=True)
    best_period = 2 * np.pi / freq_grid[np.argmax(power)]

    # Refine with a full Keplerian least-squares fit seeded at the LS period.
    # Bounds are enforced directly on the optimizer's parameters (not just
    # inside the model function), so the fitted eccentricity used later for
    # the mass calculation can't land outside the physical range even if
    # the optimizer's search path briefly considers points outside it.
    def model(t, period, t0, k, ecc, omega, gamma):
        return keplerian_rv(t, period, t0, k, ecc, omega, gamma)

    p0 = [best_period, time[np.argmax(rv)], (rv.max() - rv.min()) / 2, 0.05, 0.0, 0.0]
    bounds_lo = [0.5, 0.0, 0.0, 0.0, -2 * np.pi, -50.0]
    bounds_hi = [50.0, BASELINE_DAYS, 500.0, 0.9, 2 * np.pi, 50.0]
    popt, pcov = curve_fit(model, time, rv, p0=p0, sigma=rv_err, bounds=(bounds_lo, bounds_hi), maxfev=20000)
    fit_period, fit_t0, fit_k, fit_ecc, fit_omega, fit_gamma = popt
    perr = np.sqrt(np.diag(pcov))

    mp_recovered = minimum_mass_mearth(fit_period, fit_k, fit_ecc, STAR_MASS_MSUN)
    mp_true = minimum_mass_mearth(TRUE_PERIOD_DAYS, TRUE_K_MS, TRUE_ECC, STAR_MASS_MSUN)

    period_error_pct = abs(fit_period - TRUE_PERIOD_DAYS) / TRUE_PERIOD_DAYS * 100
    k_error_pct = abs(fit_k - TRUE_K_MS) / TRUE_K_MS * 100
    mass_error_pct = abs(mp_recovered - mp_true) / mp_true * 100

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "injected", "recovered", "error_pct"])
        writer.writerow(["period_days", TRUE_PERIOD_DAYS, f"{fit_period:.4f} +/- {perr[0]:.4f}", f"{period_error_pct:.2f}"])
        writer.writerow(["K_ms", TRUE_K_MS, f"{fit_k:.2f} +/- {perr[2]:.2f}", f"{k_error_pct:.2f}"])
        writer.writerow(["eccentricity", TRUE_ECC, f"{fit_ecc:.3f}", "-"])
        writer.writerow(["Mp_sini_mearth", f"{mp_true:.1f}", f"{mp_recovered:.1f}", f"{mass_error_pct:.2f}"])
        writer.writerow(["n_observations", N_OBS, len(time), "target vs. actual fitted sample size"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    trial_periods = 2 * np.pi / freq_grid
    axes[0].plot(trial_periods, power, color="#2f6f4f", lw=0.7)
    axes[0].axvline(TRUE_PERIOD_DAYS, color="#a8431f", ls="--", lw=1.2, label=f"Injected period = {TRUE_PERIOD_DAYS} d")
    axes[0].set_xlim(1, 20)
    axes[0].set_xlabel("Trial period [days]")
    axes[0].set_ylabel("Lomb-Scargle power")
    axes[0].set_title("Periodogram")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.25)

    phase = ((time - fit_t0) / fit_period) % 1.0
    t_model = np.linspace(0, fit_period, 500)
    rv_model = keplerian_rv(t_model, fit_period, fit_t0, fit_k, fit_ecc, fit_omega, fit_gamma)
    order = np.argsort(t_model % fit_period / fit_period)
    axes[1].errorbar(phase, rv, yerr=rv_err, fmt="o", ms=4, color="#1f4e79", capsize=2, label="Simulated RV data")
    axes[1].plot((t_model % fit_period) / fit_period, rv_model, ".", ms=1.5, color="#a8431f", label="Fitted Keplerian orbit")
    axes[1].set_xlabel("Orbital phase")
    axes[1].set_ylabel("Radial velocity [m/s]")
    axes[1].set_title(f"Phase-folded RV curve (K = {fit_k:.1f} m/s)")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.25)

    fig.suptitle("Radial velocity: recovering an injected Keplerian orbit")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "radial_velocity_recovery.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'radial_velocity_recovery.png'}")
    print(f"Injected period {TRUE_PERIOD_DAYS} d -> recovered {fit_period:.4f} d ({period_error_pct:.2f}% error)")
    print(f"Injected K {TRUE_K_MS} m/s -> recovered {fit_k:.2f} m/s ({k_error_pct:.2f}% error)")
    print(f"Injected Mp sin i {mp_true:.1f} Mearth -> recovered {mp_recovered:.1f} Mearth ({mass_error_pct:.2f}% error)")


if __name__ == "__main__":
    main()
