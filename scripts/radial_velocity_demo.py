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
import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.signal import lombscargle

plt.rcParams.update(
    {
        "axes.spines.right": False,
        "axes.spines.top": False,
        "figure.facecolor": "white",
        "font.size": 9,
        "savefig.bbox": "tight",
    }
)

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

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
FREQUENCY_GRID = np.linspace(2 * np.pi / 20.0, 2 * np.pi / 1.0, 6000)
NULL_TRIALS = 1000
RECOVERY_TRIALS = 100
RECOVERY_AMPLITUDES_MS = (1.5, 3.0, 5.0, 10.0, 20.0, 92.0)


def solve_kepler(mean_anomaly: np.ndarray, ecc: float) -> np.ndarray:
    ecc_anomaly = mean_anomaly.copy()
    for _ in range(50):
        ecc_anomaly -= (ecc_anomaly - ecc * np.sin(ecc_anomaly) - mean_anomaly) / (1 - ecc * np.cos(ecc_anomaly))
    return ecc_anomaly


def keplerian_rv(
    time: np.ndarray,
    period: float,
    t0: float,
    k: float,
    ecc: float,
    omega: float,
    gamma: float,
) -> np.ndarray:
    mean_anomaly = 2 * np.pi * ((time - t0) / period % 1.0)
    ecc_anomaly = solve_kepler(mean_anomaly, ecc)
    true_anomaly = 2 * np.arctan2(
        np.sqrt(1 + ecc) * np.sin(ecc_anomaly / 2),
        np.sqrt(1 - ecc) * np.cos(ecc_anomaly / 2),
    )
    return k * (np.cos(true_anomaly + omega) + ecc * np.cos(omega)) + gamma


def minimum_mass_mearth(
    period_days: float, k_ms: float, ecc: float, star_mass_msun: float
) -> float:
    period_s = period_days * DAY_S
    star_mass_kg = star_mass_msun * M_SUN

    def residual(mp_kg: float) -> float:
        predicted_k = (
            (2 * np.pi * G / period_s) ** (1 / 3)
            * mp_kg
            / (star_mass_kg + mp_kg) ** (2 / 3)
            / np.sqrt(1 - ecc**2)
        )
        return predicted_k - k_ms

    mp_kg = brentq(residual, 1e22, 1e29)
    return mp_kg / 5.972e24


def make_cadence(rng: np.random.Generator) -> np.ndarray:
    """Draw the frozen 60-visit two-season teaching cadence."""
    candidate_days = np.sort(rng.uniform(0, BASELINE_DAYS, N_OBS * 3))
    in_season = (candidate_days % 100) < 60
    return np.sort(rng.choice(candidate_days[in_season], size=N_OBS, replace=False))


def periodogram(time: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Return normalized power, maximum power, and its trial period."""
    power = lombscargle(time, values - values.mean(), FREQUENCY_GRID, normalize=True)
    best = int(np.argmax(power))
    return power, float(power[best]), float(2 * np.pi / FREQUENCY_GRID[best])


def wilson_interval(
    successes: int, trials: int, z: float = 1.959963984540054
) -> tuple[float, float]:
    fraction = successes / trials
    denominator = 1 + z**2 / trials
    centre = (fraction + z**2 / (2 * trials)) / denominator
    half_width = z * np.sqrt(fraction * (1 - fraction) / trials + z**2 / (4 * trials**2)) / denominator
    return max(0.0, centre - half_width), min(1.0, centre + half_width)


def calibrate_null(
    time: np.ndarray,
    noise_ms: float,
    rng: np.random.Generator,
    trials: int = NULL_TRIALS,
) -> np.ndarray:
    """Maximum LS power under white Gaussian noise on the exact cadence."""
    maxima = np.empty(trials)
    for index in range(trials):
        _, maxima[index], _ = periodogram(time, rng.normal(0, noise_ms, len(time)))
    return maxima


def recovery_experiment(
    time: np.ndarray,
    noise_ms: float,
    threshold: float,
    rng: np.random.Generator,
    trials: int = RECOVERY_TRIALS,
) -> list[dict[str, float | int | str]]:
    """Fixed-cadence period recovery conditional on this injection family."""
    rows = []
    for amplitude in RECOVERY_AMPLITUDES_MS:
        injected = keplerian_rv(
            time, TRUE_PERIOD_DAYS, TRUE_T0_DAYS, amplitude, TRUE_ECC, TRUE_OMEGA_RAD, GAMMA_MS
        )
        successes = 0
        recovered_periods = []
        for _ in range(trials):
            _, peak_power, recovered_period = periodogram(
                time, injected + rng.normal(0, noise_ms, len(time))
            )
            recovered_periods.append(recovered_period)
            if peak_power >= threshold and abs(recovered_period / TRUE_PERIOD_DAYS - 1) <= 0.01:
                successes += 1
        low, high = wilson_interval(successes, trials)
        rows.append(
            {
                "K_ms": amplitude,
                "trials": trials,
                "recoveries": successes,
                "recovery_fraction": successes / trials,
                "wilson95_low": low,
                "wilson95_high": high,
                "median_best_period_days": float(np.median(recovered_periods)),
                "criterion": (
                    "max power >= cadence null 99th percentile and "
                    "|Pbest/Pinj-1| <= 1%"
                ),
            }
        )
    return rows


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(seed=7)

    # Irregular real-campaign-like sampling: random nights within observing
    # windows, with seasonal gaps (no data for ~40% of the baseline).
    time = make_cadence(rng)

    true_rv = keplerian_rv(
        time,
        TRUE_PERIOD_DAYS,
        TRUE_T0_DAYS,
        TRUE_K_MS,
        TRUE_ECC,
        TRUE_OMEGA_RAD,
        GAMMA_MS,
    )
    total_noise = np.sqrt(RV_NOISE_MS**2 + JITTER_MS**2)
    rv = true_rv + rng.normal(0, total_noise, size=time.size)
    rv_err = np.full_like(rv, total_noise)

    # Lomb-Scargle periodogram to find the orbital period.
    power, observed_peak_power, best_period = periodogram(time, rv)

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
    popt, pcov = curve_fit(
        model,
        time,
        rv,
        p0=p0,
        sigma=rv_err,
        absolute_sigma=True,
        bounds=(bounds_lo, bounds_hi),
        maxfev=20000,
    )
    fit_period, fit_t0, fit_k, fit_ecc, fit_omega, fit_gamma = popt
    perr = np.sqrt(np.diag(pcov))

    mp_recovered = minimum_mass_mearth(fit_period, fit_k, fit_ecc, STAR_MASS_MSUN)
    mp_true = minimum_mass_mearth(TRUE_PERIOD_DAYS, TRUE_K_MS, TRUE_ECC, STAR_MASS_MSUN)

    period_error_pct = abs(fit_period - TRUE_PERIOD_DAYS) / TRUE_PERIOD_DAYS * 100
    k_error_pct = abs(fit_k - TRUE_K_MS) / TRUE_K_MS * 100
    mass_error_pct = abs(mp_recovered - mp_true) / mp_true * 100

    null_maxima = calibrate_null(time, total_noise, rng)
    null_threshold_99 = float(np.quantile(null_maxima, 0.99, method="higher"))
    null_exceedances = int(np.sum(null_maxima >= observed_peak_power))
    fap_plus_one = (null_exceedances + 1) / (NULL_TRIALS + 1)
    recovery_rows = recovery_experiment(time, total_noise, null_threshold_99, rng)

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "injected", "recovered", "error_pct"])
        writer.writerow(["period_days", TRUE_PERIOD_DAYS, f"{fit_period:.4f} +/- {perr[0]:.4f}", f"{period_error_pct:.2f}"])
        writer.writerow(["K_ms", TRUE_K_MS, f"{fit_k:.2f} +/- {perr[2]:.2f}", f"{k_error_pct:.2f}"])
        writer.writerow(["eccentricity", TRUE_ECC, f"{fit_ecc:.3f}", "-"])
        writer.writerow(["Mp_sini_mearth", f"{mp_true:.1f}", f"{mp_recovered:.1f}", f"{mass_error_pct:.2f}"])
        writer.writerow(["n_observations", N_OBS, len(time), "target vs. actual fitted sample size"])
        writer.writerow(["periodogram_peak_power", "noise-only null", f"{observed_peak_power:.8f}", "-"])
        writer.writerow(["null_99pct_max_power", f"{NULL_TRIALS} white-noise trials", f"{null_threshold_99:.8f}", "-"])
        writer.writerow(["mc_false_alarm_exceedances", NULL_TRIALS, null_exceedances, "-"])
        writer.writerow(["mc_fap_plus_one", "cadence + known white noise", f"{fap_plus_one:.8f}", "-"])

    sensitivity_path = FIG_DIR / "amplitude_recovery.csv"
    with sensitivity_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=recovery_rows[0])
        writer.writeheader()
        writer.writerows(recovery_rows)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))

    trial_periods = 2 * np.pi / FREQUENCY_GRID
    axes[0, 0].plot(trial_periods, power, color="#2f6f4f", lw=0.7)
    axes[0, 0].axhline(
        null_threshold_99,
        color="#596675",
        ls=":",
        lw=1.1,
        label="99% noise-only max-power threshold",
    )
    axes[0, 0].axvline(
        TRUE_PERIOD_DAYS,
        color="#a8431f",
        ls="--",
        lw=1.2,
        label=f"Injected period = {TRUE_PERIOD_DAYS} d",
    )
    axes[0, 0].set(
        xlim=(1, 20),
        xlabel="Trial period [days]",
        ylabel="Lomb–Scargle power",
        title="Periodogram and cadence-calibrated threshold",
    )
    axes[0, 0].legend(fontsize=7)
    axes[0, 0].grid(alpha=0.25)

    phase = ((time - fit_t0) / fit_period) % 1.0
    phase_model = np.linspace(0, 1, 500)
    t_model = fit_t0 + phase_model * fit_period
    rv_model = keplerian_rv(t_model, fit_period, fit_t0, fit_k, fit_ecc, fit_omega, fit_gamma)
    axes[0, 1].errorbar(
        phase,
        rv,
        yerr=rv_err,
        fmt="o",
        ms=4,
        color="#1f4e79",
        capsize=2,
        label="Simulated RV data",
    )
    axes[0, 1].plot(
        phase_model,
        rv_model,
        ".",
        ms=1.5,
        color="#a8431f",
        label="Fitted Keplerian orbit",
    )
    axes[0, 1].set(
        xlabel="Orbital phase",
        ylabel="Radial velocity [m/s]",
        title=f"Phase-folded fit (K = {fit_k:.1f} m/s)",
    )
    axes[0, 1].legend(fontsize=7)
    axes[0, 1].grid(alpha=0.25)

    window_power = np.abs(np.exp(1j * np.outer(FREQUENCY_GRID, time)).mean(axis=1)) ** 2
    axes[1, 0].plot(trial_periods, window_power, color="#735c8f", lw=0.7)
    axes[1, 0].set(
        xlim=(1, 20),
        xlabel="Trial period [days]",
        ylabel="Spectral-window power",
        title="Sampling window (alias structure)",
    )
    axes[1, 0].grid(alpha=0.25)

    amplitudes = np.array([row["K_ms"] for row in recovery_rows])
    fractions = np.array([row["recovery_fraction"] for row in recovery_rows])
    lower = np.array([row["wilson95_low"] for row in recovery_rows])
    upper = np.array([row["wilson95_high"] for row in recovery_rows])
    axes[1, 1].errorbar(
        amplitudes,
        fractions,
        yerr=[fractions - lower, upper - fractions],
        fmt="o-",
        color="#2f6f4f",
        capsize=3,
    )
    axes[1, 1].axvline(total_noise, color="#596675", ls=":", label="Per-visit σ")
    axes[1, 1].set(
        xscale="log",
        ylim=(-0.04, 1.04),
        xlabel="Injected K [m/s]",
        ylabel="Period recovery fraction",
        title="Fixed-cadence amplitude sensitivity",
    )
    axes[1, 1].legend(fontsize=7)
    axes[1, 1].grid(alpha=0.25)

    fig.suptitle("Radial velocity: recovering an injected Keplerian orbit")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "radial_velocity_recovery.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {sensitivity_path}")
    print(f"Wrote {FIG_DIR / 'radial_velocity_recovery.png'}")
    print(f"Injected period {TRUE_PERIOD_DAYS} d -> recovered {fit_period:.4f} d ({period_error_pct:.2f}% error)")
    print(f"Injected K {TRUE_K_MS} m/s -> recovered {fit_k:.2f} m/s ({k_error_pct:.2f}% error)")
    print(f"Injected Mp sin i {mp_true:.1f} Mearth -> recovered {mp_recovered:.1f} Mearth ({mass_error_pct:.2f}% error)")
    print(f"Cadence-calibrated null: {null_exceedances}/{NULL_TRIALS} maxima reached the observed "
          f"power; plus-one FAP={fap_plus_one:.4g}")
    print("Recovery fractions: " + ", ".join(
        f"K={row['K_ms']:g}: {row['recoveries']}/{row['trials']}" for row in recovery_rows
    ))


if __name__ == "__main__":
    main()
