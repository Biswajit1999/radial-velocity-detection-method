# Radial Velocity — Exoplanet Detection Method

The method that found the first exoplanet around a Sun-like star: watch
a star's spectral lines shift back and forth as an orbiting planet's
gravity tugs it in a small, periodic wobble. This repo explains the
physics and implements a real Lomb-Scargle periodogram plus a full
Keplerian orbit fit in Python from scratch, validated by injecting a
known signal and recovering it.

## The physics

A planet doesn't orbit a fixed star — star and planet both orbit their
common center of mass. The star's resulting reflex velocity is
periodically Doppler-shifted toward and away from Earth, imprinted as a
tiny shift in the star's absorption-line wavelengths. The velocity
semi-amplitude of that wobble is:

$$K = \left(\frac{2\pi G}{P}\right)^{1/3} \frac{M_p \sin i}{(M_\star + M_p)^{2/3}} \frac{1}{\sqrt{1-e^2}}$$

where $P$ is the orbital period, $M_\star$ and $M_p$ the star and
planet masses, $e$ the orbital eccentricity, and $i$ the (usually
unknown) orbital inclination — which is why radial velocity alone gives
only a **minimum mass**, $M_p \sin i$, not the true mass. Jupiter
induces about 12.5 m/s on the Sun; Earth induces about 9 cm/s — far
below what any instrument can currently measure for an Earth twin
around a Sun-like star.

## Why this method matters

Radial velocity was how 51 Pegasi b, the first exoplanet found around a
Sun-like star, was discovered in 1995 — and unlike transit photometry,
it works regardless of orbital inclination (as long as it's not
perfectly face-on) and directly constrains mass rather than radius,
making it the essential complement to transit surveys for measuring
real planet densities and bulk compositions.

**Real limitation:** the measured wobble scales with $1/\sqrt{M_\star}$
and favors massive, close-in planets around bright, quiet, slowly
rotating stars — real stellar "jitter" from spots, plages, and granulation
(often 1-5 m/s or more) is a real, fundamental noise floor that no
amount of instrumental precision alone can remove.

## What this repo's code does

`scripts/radial_velocity_demo.py`:

1. Injects a known Keplerian orbit (period, semi-amplitude K,
   eccentricity) sampled the way a **real ground-based spectrograph
   campaign actually observes** — irregular nightly cadence with real
   seasonal observing-window gaps, not a smooth, evenly sampled curve.
2. Adds noise combining **HARPS's own published ~1 m/s single-
   measurement photon-noise precision** with a realistic ~1.5 m/s
   stellar jitter term added in quadrature — both real, published
   regimes.
3. Recovers the orbital period from a Lomb-Scargle periodogram, then
   refines period, K, eccentricity, and argument of periastron with a
   full Keplerian least-squares fit (solving Kepler's equation via
   Newton-Raphson).
4. Converts the recovered K into a minimum mass $M_p \sin i$ using the
   real formula above and a real-like host star mass, and reports the
   error against the known injected "ground truth."

Run it yourself:

```bash
pip install -r requirements.txt
python scripts/radial_velocity_demo.py
```

## Result

| Quantity | Injected | Recovered | Error |
|---|---|---|---|
| Period | 4.23 days | 4.2303 days | 0.01% |
| K (semi-amplitude) | 92.0 m/s | 92.00 m/s | 0.01% |
| Mp sin i | 235.7 Earth masses | 235.8 Earth masses | 0.00% |

With realistic HARPS-class noise and only 60 irregularly sampled
nights, the periodogram cleanly and unambiguously identifies the true
period, and the phase-folded Keplerian fit recovers the orbit to
well under 1% error on every parameter.

## Why this repo uses simulated (not raw archival) data

This repo demonstrates the *method itself* — its sensitivity to
sampling, noise, and orbital geometry — which is best shown with a
known "ground truth" to validate recovery against. This portfolio's
companion `*-exoplanet-report` repositories instead each analyze one
real target's actual archival JWST/HST/Spitzer/ground-based spectra
directly, with zero simulated data. Both approaches are stated plainly
here rather than blurring the two.

## Repository structure

```text
scripts/radial_velocity_demo.py   Keplerian RV model + periodogram + fit + injection-recovery test
figures/                          generated plot + summary_statistics.csv
```

## References

1. Mayor, M. and Queloz, D., 1995. A Jupiter-mass companion to a
   solar-type star. *Nature*, 378, pp.355-359 — the first exoplanet
   discovered around a Sun-like star (51 Pegasi b).
2. Marcy, G.W. and Butler, R.P., 1998. Detection of Extrasolar Giant
   Planets. *Annual Review of Astronomy and Astrophysics*, 36,
   pp.57-97.
3. Lovis, C. and Fischer, D., 2010. Radial Velocity Techniques for
   Exoplanets, in *Exoplanets*, ed. S. Seager, University of Arizona
   Press, pp.27-53.
4. Lomb, N.R., 1976. Least-squares frequency analysis of unequally
   spaced data. *Astrophysics and Space Science*, 39, pp.447-462.
5. Mayor, M. et al., 2003. Setting New Standards with HARPS. *The
   Messenger*, 114, pp.20-24 — real ~1 m/s instrumental precision.
6. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
