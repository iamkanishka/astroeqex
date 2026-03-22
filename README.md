# AstroEquations

[![Hex.pm](https://img.shields.io/hexpm/v/astroequations.svg)](https://hex.pm/packages/astroequations)
[![Hex Docs](https://img.shields.io/badge/hex-docs-lightgreen.svg)](https://hexdocs.pm/astroequations)
[![CI](https://github.com/your-org/astroequations/actions/workflows/ci.yml/badge.svg)](https://github.com/your-org/astroequations/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Elixir](https://img.shields.io/badge/elixir-~%3E%201.18-purple.svg)](https://elixir-lang.org)

A comprehensive Elixir library of scientific and astronomical equations.

**696+ pure functions** across 23 modules — all in SI units, fully typed,
fully documented, zero runtime dependencies.

---

## Modules

| Namespace | Module | Fns | Description |
|---|---|---:|---|
| Astrophysics | `Astrometry` | 32 | Redshift, magnitudes, flux, parallax, distances |
| | `BlackHole` | 15 | Schwarzschild radius, Hawking, ISCO, Kerr, ergosphere |
| | `Galaxies` | 22 | Profiles, rotation curves, NFW, SFR, Schechter |
| | `Instrumentation` | 30 | Optics, AO, SNR, spectroscopy, plate scale |
| | `Stars` | 27 | Structure, timescales, Eddington, Jeans, pulsars |
| Physics | `Electromagnetism` | 85 | Maxwell, Lorentz, circuits, EM waves, Hall effect |
| | `Energy` | 18 | Work, KE, PE, power, SHM energy, E = mc² |
| | `Forces` | 22 | Newton's laws, drag, buoyancy, friction, pressure |
| | `GeneralRelativity` | 25 | Metrics, geodesics, gravitational waves, cosmology |
| | `Materials` | 28 | Stress, strain, viscosity, sound speed, EOS |
| | `Motion` | 54 | Kinematics, rotation, moments of inertia, Lagrangian |
| | `NewtonGravity` | 18 | Gravity, orbits, Roche limit, Hill sphere |
| | `Oscillations` | 31 | SHM, damping, resonance, LC circuit, beats |
| | `QuantumMechanics` | 37 | Uncertainty, de Broglie, hydrogen atom, operators |
| | `SpecialRelativity` | 26 | Lorentz factor, four-vectors, proper time |
| | `Thermodynamics` | 42 | Ideal gas, entropy, Carnot, Planck, Maxwell-Boltzmann |
| | `Waves` | 39 | Diffraction, Doppler, standing waves, Snell's Law |
| Mathematics | `Calculus` | 16 | Differentiation, integration, root finding, ODEs |
| | `Geometry` | 25 | Angular separation, orbital geometry, lensing |
| | `Notation` | 58 | Physical constants, unit conversions, Julian dates |
| | `Trigonometry` | 14 | Spherical trig, altitude/azimuth, sidereal time |
| Statistics | `StandardDeviation` | 18 | Std dev, z-score, error propagation, photometric SNR |
| | `Variance` | 14 | Variance, Welford, Pearson/Spearman r, χ², IQR |

---

## Installation

```elixir
# mix.exs
def deps do
  [
    {:astroequations, "~> 0.3"}
  ]
end
```

---

## Quick Start

```elixir
alias AstroEquations.Physics.{Energy, SpecialRelativity, Thermodynamics}
alias AstroEquations.AstrophysicsAndAstronomy.{Astrometry, BlackHole, Stars}
alias AstroEquations.Mathematics.Notation

# Classical mechanics
Energy.kinetic_energy(1.0, 10.0)              #=> 50.0  (J)
Energy.escape_velocity(5.972e24, 6.371e6)     #=> 11_185.7  (m/s)

# Special Relativity
c = Notation.speed_of_light()
SpecialRelativity.gamma_factor(0.9 * c)       #=> 2.294
SpecialRelativity.time_dilation(1.0, 0.9 * c) #=> 2.294  (seconds)

# Thermodynamics
Thermodynamics.carnot_efficiency(800, 300)    #=> 0.625
Thermodynamics.wiens_displacement(5778)       #=> 5.015e-7  (m — Sun's peak)

# Astrophysics
BlackHole.schwarzschild_radius(1.989e30)      #=> 2953.25  (m)
BlackHole.hawking_temperature(1.989e30)       #=> ~6.17e-8  (K)
Stars.nuclear_timescale(1.0)                  #=> 1.0e10  (yr — Sun's lifetime)
Astrometry.luminosity_distance(1.0)           #=> 857.6  (Mpc at z = 1)

# Unit conversions
Notation.pc_to_m(1)       #=> 3.085_677_581e16  (m)
Notation.ev_to_j(13.6)    #=> 2.179e-18          (J)
Notation.solar_mass()     #=> 1.988_92e30         (kg)
```

---

## Design

### Pure functions
Every function is a pure mathematical transformation — no state, no side
effects. Safe in any concurrent or distributed context.

### SI units throughout
All inputs and outputs use SI base units unless explicitly stated otherwise
in the function documentation (e.g. `nuclear_timescale` returns years,
`recession_velocity` uses km/s for conventional astronomy usage).

### Input guards
Physics-critical functions carry `when is_positive(mass)` style guards that
reject nonsensical inputs at the boundary with a `FunctionClauseError` rather
than silently producing NaN or incorrect results downstream.

```elixir
# These raise FunctionClauseError immediately:
Energy.kinetic_energy(-1.0, 10.0)      # negative mass
BlackHole.schwarzschild_radius(0)      # zero mass
Oscillations.spring_period(-5.0, 1.0) # negative spring constant
```

### Composable guards
Every module exports three reusable guards you can import into your own code:

```elixir
import AstroEquations.Physics.Energy,
  only: [is_positive: 1, is_non_negative: 1, is_real: 1]

def specific_intensity(flux, area) when is_positive(flux) and is_positive(area) do
  flux / area
end
```

### Fully typed
Every public function has `@spec`, `@doc`, `@doc since:`, and domain-typed
parameters. Every module exports `@type` definitions with `@typedoc`. Dialyzer
and ExDoc both have full coverage.

---

## Quality and CI

```bash
# Full quality check (format + credo + dialyzer + tests)
mix check

# Auto-fix formatting then strict lint
mix lint

# Run tests with coverage report
mix coveralls.html

# Tests for a single namespace
mix test test/physics/

# Only fast unit tests (skip @tag :property)
mix test --exclude property

# CI pipeline (same as check but also verifies no unused deps)
mix ci
```

---

## Contributing

1. Fork and create a feature branch.
2. Follow the conventions in existing modules:
   - Every new public function needs `@spec`, `@doc`, `@doc since:`, at least
     one `## Examples` block, and appropriate input guards.
   - Large number literals must use underscore separators (`2.997_924_58e8`).
   - Zero-arity functions must not use parentheses in their `def` line.
   - Predicate functions must not start with `is_` (use `foo?` not `is_foo?`).
3. Run `mix check` — all checks must pass.
4. Update `CHANGELOG.md` under `[Unreleased]`.

---

## License

Apache-2.0