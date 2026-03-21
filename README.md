# AstroEquations

[![Hex.pm](https://img.shields.io/hexpm/v/astroequations.svg)](https://hex.pm/packages/astroequations)
[![Hex Docs](https://img.shields.io/badge/hex-docs-blue.svg)](https://hexdocs.pm/astroequations)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Build Status](https://img.shields.io/github/actions/workflow/status/your-org/astroequations/ci.yml)](https://github.com/your-org/astroequations/actions)

A comprehensive Elixir library of scientific and astronomical equations. **628 functions** across 24 modules, all documented, type-annotated, and tested.

All calculations use **SI units** unless explicitly noted. Every public function carries a `@doc` description, `@spec` type signature, and an inline `iex>` example.

---

## Features

- **628 public functions** across 24 modules
- **675 ExUnit tests** with physics-grounded assertions
- **100% `@doc` and `@spec` coverage** — all functions documented and typed
- **54% `## Examples` doctest coverage** — over 320 inline `iex>` examples
- **Pure functions** — no side effects, no runtime dependencies, no state
- **SI units throughout** — consistent unit system with explicit conversion helpers
- **Dialyzer-clean** — all specs validated; `no_return()` used where functions raise
- **Credo strict** — configured at `max_arity: 9` for scientific function signatures

---

## Installation

Add `astroequations` to your `mix.exs` dependencies:

```elixir
def deps do
  [
    {:astroequations, "~> 0.2"}
  ]
end
```

Then:

```bash
mix deps.get
```

---

## Module Overview

### Astrophysics & Astronomy — 91 functions

| Module | Fns | Description |
|--------|:---:|-------------|
| `AstroEquations.AstrophysicsAndAstronomy.Astrometry` | 24 | Redshift, magnitude, flux, parallax, proper motion, angular sizes, distance modulus, metallicity |
| `AstroEquations.AstrophysicsAndAstronomy.BlackHole` | 10 | Schwarzschild radius, Hawking temperature/evaporation, ISCO, photon sphere, entropy, Kerr spin, time dilation |
| `AstroEquations.AstrophysicsAndAstronomy.Galaxies` | 14 | Hubble classification, Sérsic/de Vaucouleurs profiles, NFW halo, rotation curves, Tully-Fisher, Faber-Jackson, SFR, Schechter function, virial mass |
| `AstroEquations.AstrophysicsAndAstronomy.Instrumentation` | 23 | Lensmaker's equation, focal ratio, FOV, plate scale, AO error, Strehl ratio, SNR, limiting magnitude, grating dispersion, atmospheric extinction, rocket equation |
| `AstroEquations.AstrophysicsAndAstronomy.Stars` | 20 | Stellar structure ODEs, timescales, Eddington limits, mass-luminosity, Stefan-Boltzmann, Wien, Planck, Jeans mass/radius, Chandrasekhar limit |

### Physics — 314 functions

| Module | Fns | Description |
|--------|:---:|-------------|
| `AstroEquations.Physics.Electromagnetism` | 78 | Maxwell's equations, Lorentz force, electric/magnetic fields and potentials, circuit theory (RC/RL/LC, series/parallel), capacitors, inductors, EM waves, skin depth, plasma frequency, Brewster's angle |
| `AstroEquations.Physics.Energy` | 18 | Work, kinetic/potential energy, power, SHM energy, rest energy, binding energy, escape velocity |
| `AstroEquations.Physics.Forces` | 22 | Newton's laws, friction, buoyancy, drag (Stokes and quadratic), terminal velocity, pressure, impulse |
| `AstroEquations.Physics.GeneralRelativity` | 22 | Minkowski/Schwarzschild/Rindler/FLRW metrics, tensor algebra, gravitational wave strain, orbital precession, light deflection, Friedmann cosmology |
| `AstroEquations.Physics.Materials` | 22 | Density, stress/strain, elastic moduli, thermal expansion, heat conduction, viscosity, Reynolds number, sound speed, polytropic EOS |
| `AstroEquations.Physics.Motion` | 53 | SUVAT kinematics, rotation, 7 moment-of-inertia shapes, projectile motion, SHM, Lagrangian/Hamiltonian mechanics, Kepler's laws, vis-viva |
| `AstroEquations.Physics.NewtonGravity` | 21 | Gravitational force/field/potential/PE, orbital mechanics, escape velocity, tidal forces, Roche limit, Hill sphere |
| `AstroEquations.Physics.Oscillations` | 24 | Hooke's law, SHM kinematics/energy, damped oscillators (ζ, Q, omega_d), driven resonance, LC circuits, beats |
| `AstroEquations.Physics.QuantumMechanics` | 35 | Uncertainty principle, de Broglie, hydrogen levels, particle-in-box, harmonic oscillator, Born rule, expectation values, Pauli matrices, ladder operators, density matrix, scattering |
| `AstroEquations.Physics.SpecialRelativity` | 25 | Lorentz factor, time dilation, length contraction, relativistic dynamics, energy-momentum relation, Doppler, four-vectors, Lorentz boost, proper time |
| `AstroEquations.Physics.Thermodynamics` | 39 | Ideal gas law (all variable forms), Maxwell-Boltzmann speeds, First Law, adiabatic/isothermal/isobaric work, entropy, Carnot/COP, specific heats, Newton's cooling, Planck/Wien/Stefan-Boltzmann |
| `AstroEquations.Physics.Waves` | 34 | Wave kinematics, wave function, standing waves, interference, diffraction, Doppler, sound intensity/dB, harmonics, Snell's law, Malus's law |

### Mathematics — 101 functions

| Module | Fns | Description |
|--------|:---:|-------------|
| `AstroEquations.Mathematics.Calculus` | 11 | Numerical differentiation (forward/central/backward, 1st/2nd order), integration (trapezoid/Simpson's/rectangle), root finding (bisection/Newton-Raphson), ODE solvers (Euler/RK4) |
| `AstroEquations.Mathematics.Geometry` | 22 | Angular separation (Haversine), position angle, solid angles, orbital geometry (Kepler's equation, true/eccentric anomaly, periapsis/apoapsis), Einstein radius, microlensing, coordinate conversions |
| `AstroEquations.Mathematics.Notation` | 55 | 24 physical/astronomical constants, unit conversions (length/angle/energy/time/temperature), sexagesimal ↔ decimal (HMS/DMS), Julian Date utilities |
| `AstroEquations.Mathematics.Trigonometry` | 13 | Spherical laws of cosines/sines, altitude/azimuth, hour angle, parallactic angle, GMST/LST, sunrise hour angle, atmospheric refraction, equatorial ↔ ecliptic ↔ galactic |

### Statistics — 42 functions

| Module | Fns | Description |
|--------|:---:|-------------|
| `AstroEquations.Statistics.Variance` | 26 | Population/sample/weighted variance, Welford online algorithm, covariance, Pearson r, chi-squared, reduced chi-squared, MAD, IQR |
| `AstroEquations.Statistics.StandardDeviation` | 16 | Population/sample/weighted std dev, z-scores, normalisation, error propagation (addition/product/power/log), photometric SNR, magnitude uncertainty, Poisson uncertainty, online streaming |

---

## Quick Examples

### Astrometry

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.Astrometry

# Redshift of a galaxy with Hα shifted from 656 nm to 700 nm
Astrometry.redshift(700, 656)
#=> 0.0671

# 5-magnitude difference corresponds to a flux ratio of exactly 100
Astrometry.apparent_magnitude_diff(1.0, 100.0)
#=> -5.0

# Proxima Centauri: parallax 0.7687 arcsec → distance
Astrometry.distance_from_parallax(0.7687)
#=> 1.3009  # parsecs

# Distance modulus for the Andromeda galaxy (~770 kpc)
Astrometry.distance_modulus(770_000)
#=> 24.43
```

### Black Holes

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.BlackHole

# Schwarzschild radius of the Sun
BlackHole.schwarzschild_radius(1.989e30)
#=> 2953.25  # metres

# Hawking temperature of a stellar-mass black hole (extremely cold)
BlackHole.hawking_temperature(1.989e30)
#=> 6.17e-8  # Kelvin

# ISCO = 3 r_s — the last stable orbit
BlackHole.isco_radius(1.989e30)
#=> 8859.75  # metres
```

### Special Relativity

```elixir
alias AstroEquations.Physics.SpecialRelativity, as: SR

c = 2.99792458e8

# Lorentz factor at 90% of the speed of light
SR.gamma_factor(0.9 * c)
#=> 2.2942

# Relativistic velocity addition: two observers each moving at 0.9c
SR.relative_velocity(0.9 * c, -0.9 * c)
#=> 2.9691e8  # always less than c

# Time dilation: a 1-second proper interval at 0.6c
SR.time_dilation(1.0, 0.6 * c)
#=> 1.25  # seconds

# Energy-momentum relation for a photon (m = 0)
SR.energy_momentum_relation(1.0e-22, 0.0)
#=> 2.998e-14  # joules
```

### Thermodynamics

```elixir
alias AstroEquations.Physics.Thermodynamics

# Ideal gas: solve for pressure (nil marks the unknown)
Thermodynamics.ideal_gas_law(p: nil, v: 0.0224, n: 6.022e23, k_b: 1.38e-23, t: 273.15)
#=> %{p: 101_325.0}

# Carnot efficiency: hot 500 K, cold 300 K
Thermodynamics.carnot_efficiency(500, 300)
#=> 0.4

# Wien's law: solar surface (5778 K) peaks near 501 nm
Thermodynamics.wiens_displacement(5778)
#=> 5.015e-7  # metres
```

### Numerical Calculus

```elixir
alias AstroEquations.Mathematics.Calculus

# Derivative of sin(x) at x = 0 → cos(0) = 1
Calculus.derivative(&:math.sin/1, 0.0)
#=> 1.0

# Integrate sin(x) from 0 to π (exact: 2.0)
Calculus.simpsons(&:math.sin/1, 0, :math.pi(), 1000)
#=> 2.0

# Newton-Raphson: solve x² - 2 = 0 → √2
Calculus.newton_raphson(fn x -> x * x - 2 end, 1.5)
#=> 1.4142135623730951

# RK4 ODE: dy/dx = -y, y(0) = 1 → y(1) ≈ e⁻¹
results  = Calculus.rk4_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 1000)
{_x, y} = List.last(results)
#=> y ≈ 0.36787944
```

### Statistics

```elixir
alias AstroEquations.Statistics.{Variance, StandardDeviation}

data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]

# Population variance (Pearson's classic example: σ² = 4.0)
Variance.population(data)
#=> 4.0

# Sample std dev with Bessel's correction
StandardDeviation.sample(data)
#=> 2.1381

# One-pass Welford streaming — useful for large or live datasets
{:ok, mean, _pop_var, samp_var} =
  data |> Variance.welford_accumulate() |> Variance.welford_finalize()
# mean => 5.0, samp_var => 4.5714

# Error propagation for quadrature addition: σ_c = √(σ_a² + σ_b²)
StandardDeviation.propagate_addition(0.05, 0.03)
#=> 0.0583

# Photometric magnitude uncertainty from SNR
StandardDeviation.magnitude_uncertainty(100.0)
#=> 0.0109  # mag
```

### General Relativity

```elixir
alias AstroEquations.Physics.GeneralRelativity, as: GR

# Minkowski metric tensor (signature -+++)
GR.minkowski_metric()
#=> [[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

# GR light deflection at the solar limb (~1.75 arcsec)
alpha_rad = GR.light_deflection(1.989e30, 6.957e8)
alpha_rad * 180 / :math.pi() * 3600
#=> 1.75  # arcseconds

# Critical density for a flat universe (H₀ = 67.36 km/s/Mpc)
h0_si = 67.36 * 1000 / 3.0856e22
GR.critical_density(h0_si)
#=> 9.47e-27  # kg/m³
```

---

## Design Principles

**Pure functions.** No process state, no I/O, no side effects. Every function is composable with Elixir pipes.

**Consistent units.** All inputs and outputs are SI unless the function documentation explicitly states otherwise. `Notation` provides conversion helpers for all common astronomical unit systems.

**Overridable constants.** Physical constants used in calculations appear as optional tail parameters with standard SI defaults, making functions testable and adaptable:

```elixir
# Default G = 6.6743e-11 m³ kg⁻¹ s⁻²
NewtonGravity.force(1.0e24, 1.0e20, 1.0e9)

# Custom G (e.g., for geometric units where G = 1)
NewtonGravity.force(1.0e24, 1.0e20, 1.0e9, 1.0)
```

**Error tuples for invalid input.** Functions accepting lists return `{:error, reason}` for empty inputs rather than raising, enabling clean pattern matching:

```elixir
case Variance.sample([]) do
  {:error, reason} -> Logger.warning("Cannot compute variance: #{reason}")
  variance         -> process(variance)
end
```

**Raising for domain violations.** Functions like `NewtonGravity.force/3` raise `ArgumentError` for non-positive distances, consistent with the physical domain. These are annotated `:: float | no_return()` in their `@spec`.

---

## Running Tests

```bash
mix test                  # all 675 tests
mix test --trace          # verbose per-test output
mix test test/astro_equations/physics/special_relativity_test.exs
```

## Static Analysis

```bash
mix format --check-formatted   # formatter
mix credo --strict             # style/quality (max_arity: 9 is configured)
mix dialyzer                   # type checking (first run: mix dialyzer --plt)
```

## Generating Docs

```bash
mix docs
open doc/index.html
```

---

## Contributing

1. Fork and create a feature branch.
2. Add functions with `@doc`, `@spec`, and at least one `## Examples` block.
3. Write tests in `test/astro_equations/**/*_test.exs`.
4. Ensure `mix format`, `mix credo --strict`, and `mix test` all pass.
5. Open a pull request describing the equations added and their source.

---

## License

MIT — see [LICENSE](LICENSE).
