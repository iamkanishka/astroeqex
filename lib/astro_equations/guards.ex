defmodule AstroEquations.Guards do
  @moduledoc """
  Canonical guard predicates for the entire AstroEquations library.

  This module is the **single source of truth** for all guard macros.
  Every AstroEquations module imports from here — there are no inline
  `defguard` blocks elsewhere in the codebase.

  ## Usage

  ### Library modules (recommended)

      defmodule AstroEquations.Physics.MyModule do
        use AstroEquations.Guards

        def speed_of_light(c) when is_positive(c), do: c
      end

  ### Test files and scripts

      import AstroEquations.Guards

      assert is_positive(9.81)
      assert is_subluminal(1.0e8, 2.997_924_58e8)

  ## Guard groups

  | Group | Guards |
  |---|---|
  | Numeric primitives | `is_real/1`, `is_positive/1`, `is_non_negative/1`, `is_negative/1`, `is_valid_fraction/1` |
  | Integer | `is_positive_integer/1`, `is_non_neg_integer/1` |
  | Collection | `is_non_empty_list/1` |
  | Classical mechanics | `is_mass/1`, `is_length/1`, `is_distance/1`, `is_speed/1`, `is_pressure/1` |
  | Thermodynamics | `is_temperature/1`, `is_entropy/1` |
  | Waves & optics | `is_wavelength/1`, `is_frequency/1`, `is_angle/1` |
  | Electromagnetism | `is_permittivity/1`, `is_charge/1` |
  | Special relativity | `is_subluminal/2` |
  | Quantum mechanics | `is_quantum_number/1`, `is_probability/1` |
  | Astrophysics | `is_redshift/1`, `is_eccentricity/1`, `is_bound_orbit/1`, `is_spin_param/1`, `is_snr/1`, `is_airmass/1`, `is_magnitude/1` |
  | Multi-arg | `are_positive/2`, `are_positive/3`, `are_non_negative/2` |

  ## Design notes

  All guards are composed exclusively from Erlang's primitive guard BIFs
  (`is_number/1`, `is_integer/1`, `is_float/1`, `is_list/1`, and comparison
  operators). This ensures they are legal in `when` clauses everywhere:
  function heads, `case`, `cond`, `receive`, and `with`.

  Guards do **not** call any user-defined functions — doing so would make
  them illegal in `when` clauses and cause a `CompileError`.

  Domain-specific guards (`is_mass/1`, `is_temperature/1`, etc.) are semantic
  aliases for the underlying numeric constraint, providing self-documenting
  code without additional runtime cost.
  """

  # ---------------------------------------------------------------------------
  # __using__ macro — `use AstroEquations.Guards` in any module
  # ---------------------------------------------------------------------------

  @doc false
  defmacro __using__(_opts) do
    quote do
      import AstroEquations.Guards
    end
  end

  # ===========================================================================
  # 1. Numeric primitives
  # ===========================================================================

  @doc """
  Guards that `x` is a real number — either an integer or a float.

  This is the most permissive numeric guard; prefer `is_positive/1` or
  `is_non_negative/1` when the domain requires it.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_real(3)
      true
      iex> is_real(3.14)
      true
      iex> is_real("3")
      false
  """

  defguard is_real(x) when is_integer(x) or is_float(x)

  @doc """
  Guards that `x` is a positive real number (`x > 0`).

  The most commonly used guard in the library. Rejects zero as well as
  negative values, making it appropriate for all physical quantities that
  cannot be zero: mass, radius, temperature (Kelvin), wavelength, etc.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_positive(9.81)
      true
      iex> is_positive(0.0)
      false
      iex> is_positive(-1)
      false
  """

  defguard is_positive(x) when is_number(x) and x > 0

  @doc """
  Guards that `x` is a non-negative real number (`x >= 0`).

  Use for quantities that may legitimately be zero: displacement, height,
  speed (at rest), elapsed time, amplitude, redshift.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_non_negative(0.0)
      true
      iex> is_non_negative(5.0)
      true
      iex> is_non_negative(-0.1)
      false
  """

  defguard is_non_negative(x) when is_number(x) and x >= 0

  @doc """
  Guards that `x` is a strictly negative real number (`x < 0`).

  Use for quantities that are always negative by definition: gravitational
  potential energy, binding energy, work done against a force.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_negative(-1.0)
      true
      iex> is_negative(0.0)
      false
  """

  defguard is_negative(x) when is_number(x) and x < 0

  @doc """
  Guards that `x` is a dimensionless fraction in the closed interval `[0, 1]`.

  Use for efficiencies, quantum efficiencies, probabilities, and any ratio
  that is bounded between 0 and 1 inclusive.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_valid_fraction(0.5)
      true
      iex> is_valid_fraction(0.0)
      true
      iex> is_valid_fraction(1.0)
      true
      iex> is_valid_fraction(1.001)
      false
  """

  defguard is_valid_fraction(x) when is_number(x) and x >= 0 and x <= 1

  # ===========================================================================
  # 2. Integer guards
  # ===========================================================================

  @doc """
  Guards that `n` is a strictly positive integer (`n ≥ 1`).

  Use for principal quantum numbers, harmonic orders, diffraction orders,
  and any discrete physical quantity starting at 1.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_positive_integer(1)
      true
      iex> is_positive_integer(0)
      false
      iex> is_positive_integer(1.0)
      false
  """

  defguard is_positive_integer(n) when is_integer(n) and n > 0

  @doc """
  Guards that `n` is a non-negative integer (`n ≥ 0`).

  Use for quantum harmonic oscillator levels (`n = 0, 1, 2, …`), ladder
  operator states, and counting quantities that may be zero.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_non_neg_integer(0)
      true
      iex> is_non_neg_integer(5)
      true
      iex> is_non_neg_integer(-1)
      false
  """

  defguard is_non_neg_integer(n) when is_integer(n) and n >= 0

  # ===========================================================================
  # 3. Collection guards
  # ===========================================================================

  @doc """
  Guards that `xs` is a non-empty list.

  Use as the first precondition on all statistical functions that cannot
  operate on empty datasets (variance, standard deviation, median, etc.).

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_non_empty_list([1, 2, 3])
      true
      iex> is_non_empty_list([])
      false
      iex> is_non_empty_list("not a list")
      false
  """

  defguard is_non_empty_list(xs) when is_list(xs) and xs != []

  # ===========================================================================
  # 4. Classical mechanics
  # ===========================================================================

  @doc """
  Guards that `m` is a valid mass (`m > 0`, in kilograms).

  Mass is always strictly positive; zero or negative mass has no physical
  meaning in classical or relativistic mechanics.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_mass(9.109_383_70e-31)   # electron mass
      true
      iex> is_mass(0)
      false
  """

  defguard is_mass(m) when is_number(m) and m > 0

  @doc """
  Guards that `l` is a valid length or distance (`l > 0`, in metres).

  Lengths and distances are strictly positive: zero separation implies
  the same point, and negative lengths are not physical.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_length(6.371e6)   # Earth radius
      true
      iex> is_length(0)
      false
  """

  defguard is_length(l) when is_number(l) and l > 0

  @doc """
  Guards that `r` is a valid radial distance (`r > 0`, in metres).

  Alias for `is_length/1` with semantics specific to radial coordinates
  and separation distances.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_distance(3.844e8)   # Earth–Moon distance
      true
  """

  defguard is_distance(r) when is_number(r) and r > 0

  @doc """
  Guards that `v` is a valid speed (`v >= 0`, in m/s).

  Speed is non-negative; use `is_subluminal/2` to additionally enforce
  the relativistic speed limit.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_speed(343.0)   # speed of sound
      true
      iex> is_speed(0.0)     # at rest — valid
      true
      iex> is_speed(-1.0)
      false
  """

  defguard is_speed(v) when is_number(v) and v >= 0

  @doc """
  Guards that `p` is a valid pressure (`p > 0`, in pascals).

  Absolute pressure is always positive; gauge pressures are handled
  separately as signed quantities.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_pressure(101_325.0)   # standard atmospheric pressure
      true
      iex> is_pressure(0)
      false
  """

  defguard is_pressure(p) when is_number(p) and p > 0

  # ===========================================================================
  # 5. Thermodynamics
  # ===========================================================================

  @doc """
  Guards that `t` is a valid thermodynamic temperature (`t > 0`, in kelvin).

  Absolute temperature on the Kelvin scale is always strictly positive.
  Zero kelvin (absolute zero) is a theoretical limit, never achieved.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_temperature(5778.0)   # solar effective temperature
      true
      iex> is_temperature(0)
      false
      iex> is_temperature(-10)
      false
  """

  defguard is_temperature(t) when is_number(t) and t > 0

  @doc """
  Guards that `s` is a valid entropy value (`s >= 0`, in J/K).

  By the Third Law of Thermodynamics, entropy is non-negative.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_entropy(0.0)    # ground state
      true
      iex> is_entropy(42.5)
      true
      iex> is_entropy(-1.0)
      false
  """

  defguard is_entropy(s) when is_number(s) and s >= 0

  # ===========================================================================
  # 6. Waves & optics
  # ===========================================================================

  @doc """
  Guards that `lambda` is a valid wavelength (`lambda > 0`, in metres).

  Wavelength is strictly positive. Valid across the full electromagnetic
  spectrum: radio (m), IR (μm), visible (nm), X-ray (pm).

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_wavelength(500.0e-9)   # 500 nm — green light
      true
      iex> is_wavelength(0)
      false
  """

  defguard is_wavelength(lambda) when is_number(lambda) and lambda > 0

  @doc """
  Guards that `f` is a valid frequency (`f > 0`, in hertz).

  Frequency is strictly positive. A frequency of zero corresponds to a
  static (DC) field, not an oscillation.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_frequency(440.0)       # concert A
      true
      iex> is_frequency(5.997e14)    # visible light
      true
      iex> is_frequency(0)
      false
  """

  defguard is_frequency(f) when is_number(f) and f > 0

  @doc """
  Guards that `theta` is a valid angle (any real number, in radians).

  Angles in radians are unbounded — they can be negative, zero, or larger
  than 2π for accumulated rotations and phase differences.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_angle(0.0)
      true
      iex> is_angle(-:math.pi())
      true
      iex> is_angle(1000.0)
      true
  """

  defguard is_angle(theta) when is_number(theta)

  # ===========================================================================
  # 7. Electromagnetism
  # ===========================================================================

  @doc """
  Guards that `epsilon` is a valid permittivity (`epsilon > 0`, in F/m).

  Electric permittivity is always positive. The vacuum value is
  ε₀ ≈ 8.854_187_812_8e-12 F/m.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_permittivity(8.854_187_812_8e-12)
      true
      iex> is_permittivity(0)
      false
  """

  defguard is_permittivity(epsilon) when is_number(epsilon) and epsilon > 0

  @doc """
  Guards that `q` is a valid electric charge (any real number, in coulombs).

  Charge is signed: positive for protons, negative for electrons,
  and zero for neutral particles.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_charge(1.602_176_634e-19)    # proton charge
      true
      iex> is_charge(-1.602_176_634e-19)   # electron charge
      true
      iex> is_charge(0.0)                  # neutral
      true
  """

  defguard is_charge(q) when is_number(q)

  # ===========================================================================
  # 8. Special relativity
  # ===========================================================================

  @doc """
  Guards that speed `v` is strictly subluminal (`|v| < |c|`).

  Prevents singularities in the Lorentz factor γ = 1/√(1 − v²/c²).
  At v = c, γ diverges; at v > c, γ becomes imaginary.

  Both `v` and `c` should be in consistent units (both m/s, both km/s, etc.).
  The condition checks v² < c² to avoid a square root, supporting both
  positive and negative velocities (e.g. approaching vs receding sources).

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_subluminal(1.0e8, 2.997_924_58e8)
      true
      iex> is_subluminal(3.0e8, 2.997_924_58e8)
      false
      iex> is_subluminal(-1.0e8, 2.997_924_58e8)   # receding: OK
      true
  """

  defguard is_subluminal(v, c)
           when is_number(v) and is_number(c) and c > 0 and v * v < c * c

  # ===========================================================================
  # 9. Quantum mechanics
  # ===========================================================================

  @doc """
  Guards that `n` is a valid principal quantum number (positive integer).

  Principal quantum numbers start at 1: the hydrogen ground state is n = 1.
  For harmonic oscillator levels, which start at 0, use `is_non_neg_integer/1`.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_quantum_number(1)    # hydrogen ground state
      true
      iex> is_quantum_number(0)    # not a valid principal number
      false
      iex> is_quantum_number(1.0)  # must be an integer
      false
  """

  defguard is_quantum_number(n) when is_integer(n) and n > 0

  @doc """
  Guards that `p` is a valid probability (fraction in `[0, 1]`).

  A probability of 0 means impossible; 1 means certain.
  Used for Born rule probabilities, quantum efficiencies, and
  measurement outcomes.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_probability(0.5)
      true
      iex> is_probability(0.0)    # impossible event
      true
      iex> is_probability(1.0)    # certain event
      true
      iex> is_probability(1.001)
      false
  """

  defguard is_probability(p) when is_number(p) and p >= 0 and p <= 1

  # ===========================================================================
  # 10. Astrophysics
  # ===========================================================================

  @doc """
  Guards that `z` is a valid cosmological redshift (`z >= 0`).

  Redshift is non-negative: z = 0 means the object is at rest relative to
  the observer; z > 0 means it is receding; z < 0 (blueshift) is handled
  by the relativistic Doppler function separately.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_redshift(0.0)    # at rest
      true
      iex> is_redshift(1.0)    # distant galaxy
      true
      iex> is_redshift(-0.1)   # blueshifted — not a redshift
      false
  """

  defguard is_redshift(z) when is_number(z) and z >= 0

  @doc """
  Guards that `e` is a valid orbital eccentricity (`e >= 0`).

  Valid for all conic sections:
  - `e = 0.0` — circular orbit
  - `0 < e < 1` — elliptical orbit
  - `e = 1.0` — parabolic trajectory
  - `e > 1.0` — hyperbolic trajectory

  For the elliptical-only constraint, use `is_bound_orbit/1`.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_eccentricity(0.206)   # Mercury's eccentricity
      true
      iex> is_eccentricity(1.5)     # hyperbolic flyby
      true
      iex> is_eccentricity(-0.1)
      false
  """

  defguard is_eccentricity(e) when is_number(e) and e >= 0

  @doc """
  Guards that `e` is a valid eccentricity for a gravitationally **bound**
  elliptical orbit (`0 <= e < 1`).

  Parabolic (e = 1) and hyperbolic (e > 1) trajectories are not bound.
  Use `is_eccentricity/1` for the general case including unbound trajectories.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_bound_orbit(0.0)     # circular
      true
      iex> is_bound_orbit(0.999)   # highly elliptical
      true
      iex> is_bound_orbit(1.0)     # parabolic — not bound
      false
  """

  defguard is_bound_orbit(e) when is_number(e) and e >= 0 and e < 1

  @doc """
  Guards that `a` is a valid Kerr dimensionless spin parameter (`0 <= a* <= 1`).

  The dimensionless spin parameter a* = Jc/(GM²) is bounded:
  - `a* = 0` — Schwarzschild (non-rotating) black hole
  - `0 < a* < 1` — sub-extremal Kerr black hole
  - `a* = 1` — extremal Kerr black hole (theoretical maximum)

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_spin_param(0.0)    # Schwarzschild
      true
      iex> is_spin_param(0.998)  # near-extremal (like M87*)
      true
      iex> is_spin_param(1.0)    # extremal
      true
      iex> is_spin_param(1.001)
      false
  """

  defguard is_spin_param(a) when is_number(a) and a >= 0 and a <= 1

  @doc """
  Guards that `snr` is a valid signal-to-noise ratio (`snr >= 0`).

  SNR is non-negative: SNR = 0 means the signal is completely buried in
  noise. In practice, detection thresholds are typically SNR ≥ 3–5.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_snr(27.4)
      true
      iex> is_snr(0.0)    # pure noise
      true
      iex> is_snr(-1.0)
      false
  """

  defguard is_snr(snr) when is_number(snr) and snr >= 0

  @doc """
  Guards that `x` is a valid airmass value (`x >= 1`).

  Airmass is the optical path length through Earth's atmosphere relative
  to zenith. At zenith (z = 0°), X = 1 (minimum). It increases with
  zenith angle and diverges at the horizon.

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_airmass(1.0)    # zenith
      true
      iex> is_airmass(2.0)    # z = 60°
      true
      iex> is_airmass(0.9)    # physically impossible
      false
  """

  defguard is_airmass(x) when is_number(x) and x >= 1

  @doc """
  Guards that `m` is a valid astronomical magnitude (any real number).

  The Pogson magnitude scale is unbounded:
  - Large positive values (e.g. +30) — extremely faint objects
  - Small or negative values (e.g. −26.7) — very bright objects (Sun)
  - There is no physical upper or lower limit

  ## Examples

      iex> import AstroEquations.Guards
      iex> is_magnitude(-1.46)    # Sirius apparent magnitude
      true
      iex> is_magnitude(4.83)     # Sun's absolute magnitude
      true
      iex> is_magnitude(30.0)     # faint galaxy
      true
  """

  defguard is_magnitude(m) when is_number(m)

  # ===========================================================================
  # 11. Multi-argument conveniences
  # ===========================================================================

  @doc """
  Guards that both `a` and `b` are positive real numbers (`a > 0` and `b > 0`).

  A convenience guard for functions that take two independent positive
  physical quantities, avoiding verbose multi-clause `when` expressions.

  ## Examples

      iex> import AstroEquations.Guards
      iex> are_positive(3.0, 4.0)
      true
      iex> are_positive(0.0, 4.0)
      false
      iex> are_positive(3.0, -1.0)
      false
  """

  defguard are_positive(a, b)
           when is_number(a) and a > 0 and is_number(b) and b > 0

  @doc """
  Guards that all three of `a`, `b`, and `c` are positive real numbers.

  A convenience guard for functions with three independent positive
  physical quantities (e.g. mass, radius, temperature — or m1, m2, r).

  ## Examples

      iex> import AstroEquations.Guards
      iex> are_positive(1.0, 2.0, 3.0)
      true
      iex> are_positive(1.0, 0.0, 3.0)
      false
  """

  defguard are_positive(a, b, c)
           when is_number(a) and a > 0 and is_number(b) and b > 0 and is_number(c) and c > 0

  @doc """
  Guards that both `a` and `b` are non-negative real numbers (`a >= 0` and `b >= 0`).

  Use for pairs of quantities where zero is valid: velocities, heights,
  amplitudes, angles.

  ## Examples

      iex> import AstroEquations.Guards
      iex> are_non_negative(0.0, 5.0)
      true
      iex> are_non_negative(-1.0, 5.0)
      false
  """

  defguard are_non_negative(a, b)
           when is_number(a) and a >= 0 and is_number(b) and b >= 0
end
