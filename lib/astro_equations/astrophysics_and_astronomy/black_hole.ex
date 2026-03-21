defmodule AstroEquations.AstrophysicsAndAstronomy.BlackHole do
  @moduledoc """
  A collection of functions for calculating black hole properties.

  Covers:
  - Schwarzschild radius (event horizon) — standard and solar-mass forms
  - Hawking temperature and radiation
  - Evaporation timescale
  - Photon sphere radius
  - ISCO (Innermost Stable Circular Orbit) — Schwarzschild and Kerr
  - Bekenstein-Hawking entropy
  - Kerr spin parameter and ergosphere
  - Gravitational time dilation at the event horizon
  - Black hole shadow radius
  - Accretion luminosity and Bondi radius
  """

  # m³ kg⁻¹ s⁻²
  @gravitational_constant 6.67430e-11
  # m/s
  @speed_of_light 2.99792458e8
  # J·s (reduced Planck)
  @hbar 1.054571817e-34
  # J/K
  @boltzmann 1.380649e-23
  # kg
  @solar_mass 1.989e30

  # ---------------------------------------------------------------------------
  # Schwarzschild Radius
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the Schwarzschild radius (event horizon radius) for a given mass.

  r_s = 2 G M / c²

  ## Parameters
    - mass: Mass in kilograms

  ## Returns
    Schwarzschild radius in meters

  ## Examples
      iex> BlackHole.schwarzschild_radius(1.989e30) |> Float.round(2)
      2953.25
  """
  @spec schwarzschild_radius(number) :: float
  def schwarzschild_radius(mass) do
    2 * @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
  end

  @doc """
  Calculates the Schwarzschild radius for a mass given in solar masses.

  ## Parameters
    - solar_masses: Mass in units of solar masses (M☉)

  ## Returns
    Schwarzschild radius in meters

  ## Examples
      iex> BlackHole.schwarzschild_radius_solar(1) |> Float.round(2)
      2953.25
  """
  @spec schwarzschild_radius_solar(number) :: float
  def schwarzschild_radius_solar(solar_masses) do
    schwarzschild_radius(solar_masses * @solar_mass)
  end

  # ---------------------------------------------------------------------------
  # Hawking Radiation
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the Hawking temperature of a Schwarzschild black hole.

  T_H = ħ c³ / (8π G M k_B)

  Smaller (less massive) black holes are hotter.

  ## Parameters
    - mass: Mass in kilograms

  ## Returns
    Hawking temperature in Kelvin

  ## Examples
      iex> BlackHole.hawking_temperature(1.989e30) > 0
      true
  """
  @spec hawking_temperature(number) :: float
  def hawking_temperature(mass) do
    @hbar * :math.pow(@speed_of_light, 3) /
      (8 * :math.pi() * @gravitational_constant * mass * @boltzmann)
  end

  @doc """
  Calculates the Hawking evaporation timescale of a black hole.

  t_evap = 5120 π G² M³ / (ħ c⁴)

  ## Parameters
    - mass: Initial mass in kilograms

  ## Returns
    Evaporation timescale in seconds

  ## Examples
      iex> BlackHole.evaporation_time(1.0e10) > 0
      true
  """
  @spec evaporation_time(number) :: float
  def evaporation_time(mass) do
    5120 * :math.pi() * :math.pow(@gravitational_constant, 2) * :math.pow(mass, 3) /
      (@hbar * :math.pow(@speed_of_light, 4))
  end

  @doc """
  Estimates the minimum initial mass of a primordial black hole that has not yet
  fully evaporated by time t.

  M_min = (ħ c⁴ t / (5120 π G²))^(1/3)

  ## Parameters
    - time: Age of the universe or elapsed time (s)

  ## Returns
    Minimum surviving mass in kg

  ## Examples
      iex> BlackHole.min_surviving_pbh_mass(4.35e17) > 0
      true
  """
  @spec min_surviving_pbh_mass(number) :: float
  def min_surviving_pbh_mass(time) do
    :math.pow(
      @hbar * :math.pow(@speed_of_light, 4) * time /
        (5120 * :math.pi() * :math.pow(@gravitational_constant, 2)),
      1 / 3
    )
  end

  # ---------------------------------------------------------------------------
  # Photon Sphere & ISCO
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the photon sphere radius for a Schwarzschild black hole.

  r_ph = 3 G M / c² = 1.5 × r_s

  Light in circular orbits at this radius is unstable.

  ## Examples
      iex> BlackHole.photon_sphere_radius(1.989e30) |> Float.round(2)
      4429.88
  """
  @spec photon_sphere_radius(number) :: float
  def photon_sphere_radius(mass) do
    3 * @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
  end

  @doc """
  Calculates the ISCO radius for a Schwarzschild black hole.

  r_ISCO = 6 G M / c² = 3 × r_s

  The innermost stable circular orbit for a non-rotating black hole.

  ## Examples
      iex> BlackHole.isco_radius(1.989e30) |> Float.round(2)
      8859.75
  """
  @spec isco_radius(number) :: float
  def isco_radius(mass) do
    6 * @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
  end

  @doc """
  Calculates the ISCO radius for a Kerr (spinning) black hole.

  The prograde ISCO (co-rotating orbit) shrinks toward r_s as a* → 1,
  while the retrograde orbit increases to 9 r_g.

  r_ISCO = r_g × Z₂ ∓ √((3Z₁ − Z₂)(3Z₁ + Z₂))

  where r_g = GM/c², Z₁ and Z₂ are functions of the spin parameter.

  Uses the Bardeen et al. (1972) analytic formula.

  ## Parameters
    - mass:         Mass in kilograms
    - spin_param:   Dimensionless spin a* in [0, 1]
    - prograde:     true for co-rotating orbit (default), false for retrograde

  ## Returns
    ISCO radius in meters

  ## Examples
      iex> BlackHole.kerr_isco_radius(1.989e30, 0.0) |> Float.round(2)
      8859.75
  """
  @spec kerr_isco_radius(number, number, boolean) :: float
  def kerr_isco_radius(mass, spin_param, prograde \\ true) do
    r_g = @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
    a = spin_param

    z1 =
      1 +
        :math.pow(1 - a * a, 1 / 3) *
          (:math.pow(1 + a, 1 / 3) + :math.pow(1 - a, 1 / 3))

    z2 = :math.sqrt(3 * a * a + z1 * z1)
    sign = if prograde, do: -1, else: 1
    r_g * (3 + z2 + sign * :math.sqrt((3 - z1) * (3 + z1 + 2 * z2)))
  end

  # ---------------------------------------------------------------------------
  # Entropy
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the Bekenstein-Hawking entropy of a black hole.

  S = k_B c³ A / (4 G ħ),  where A = 4π r_s²

  The entropy is proportional to the horizon area, not the volume.

  ## Returns
    Entropy in J/K

  ## Examples
      iex> BlackHole.bekenstein_hawking_entropy(1.989e30) > 0
      true
  """
  @spec bekenstein_hawking_entropy(number) :: float
  def bekenstein_hawking_entropy(mass) do
    r_s = schwarzschild_radius(mass)
    area = 4 * :math.pi() * :math.pow(r_s, 2)
    @boltzmann * :math.pow(@speed_of_light, 3) * area / (4 * @gravitational_constant * @hbar)
  end

  # ---------------------------------------------------------------------------
  # Kerr Black Hole
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the dimensionless Kerr spin parameter a* = J c / (G M²), clamped to [0, 1].

  ## Parameters
    - mass:             Mass in kilograms
    - angular_momentum: Angular momentum J in kg·m²/s

  ## Returns
    Dimensionless spin parameter in [0.0, 1.0]

  ## Examples
      iex> BlackHole.kerr_spin_parameter(1.989e30, 1.0e47) |> Float.round(4)
      0.1195
  """
  @spec kerr_spin_parameter(number, number) :: float
  def kerr_spin_parameter(mass, angular_momentum) do
    a_star =
      angular_momentum * @speed_of_light /
        (@gravitational_constant * :math.pow(mass, 2))

    min(max(a_star, 0.0), 1.0)
  end

  @doc """
  Outer ergosphere radius of a Kerr black hole at polar angle θ.

  r_ergo = G M/c² + √((G M/c²)² − a² cos²θ)

  where a = J/(Mc) is the specific angular momentum (in metres, not dimensionless).

  At the equator (θ = π/2), cos θ = 0 and r_ergo = r_s regardless of spin.
  At the poles (θ = 0 or π), r_ergo = r_s when a = 0, and they coincide with
  the event horizon for a* = 1.

  ## Parameters
    - mass:        Mass in kilograms
    - spin_param:  Dimensionless spin a* in [0, 1]
    - theta:       Polar angle from spin axis (radians, default: π/2 = equator)

  ## Returns
    Ergosphere radius in meters

  ## Examples
      iex> BlackHole.ergosphere_radius(1.989e30, 0.5) |> Float.round(2)
      2953.25
  """
  @spec ergosphere_radius(number, number, number) :: float
  def ergosphere_radius(mass, spin_param, theta \\ :math.pi() / 2) do
    r_g = @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
    # specific angular momentum in metres
    a_m = spin_param * r_g
    cos_theta = :math.cos(theta)
    r_g + :math.sqrt(max(r_g * r_g - a_m * a_m * cos_theta * cos_theta, 0.0))
  end

  # ---------------------------------------------------------------------------
  # Black Hole Shadow
  # ---------------------------------------------------------------------------

  @doc """
  Shadow radius (apparent radius of the black hole silhouette) for a Schwarzschild BH.

  r_shadow = 3√3 G M / c² = (3√3 / 2) × r_s ≈ 2.598 × r_s

  This is the photon capture cross-section projected at infinity.

  ## Returns
    Shadow radius in meters

  ## Examples
      iex> BlackHole.bh_shadow_radius(1.989e30) > BlackHole.photon_sphere_radius(1.989e30)
      true
  """
  @spec bh_shadow_radius(number) :: float
  def bh_shadow_radius(mass) do
    3 * :math.sqrt(3) * @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
  end

  # ---------------------------------------------------------------------------
  # Accretion
  # ---------------------------------------------------------------------------

  @doc """
  Accretion luminosity from mass accretion rate: L = η ṁ c²

  η is the radiative efficiency (≈ 0.06 for Schwarzschild, ≈ 0.42 for maximal Kerr).

  ## Parameters
    - accretion_rate: Mass accretion rate ṁ in kg/s
    - efficiency:     Radiative efficiency η (default: 0.1)

  ## Returns
    Luminosity in Watts

  ## Examples
      iex> BlackHole.accretion_luminosity(1.0e20, 0.1) > 0
      true
  """
  @spec accretion_luminosity(number, number) :: float
  def accretion_luminosity(accretion_rate, efficiency \\ 0.1) do
    efficiency * accretion_rate * :math.pow(@speed_of_light, 2)
  end

  @doc """
  Bondi accretion radius: r_B = 2 G M / v_∞²

  The radius within which the black hole's gravity dominates the surrounding gas flow.

  ## Parameters
    - mass:      Black hole mass (kg)
    - v_inf:     Velocity at infinity (sound speed or bulk velocity in m/s)
    - sound_speed: Additional quadrature term (m/s, default: same as v_inf for simplicity)

  ## Returns
    Bondi radius in meters

  ## Examples
      iex> BlackHole.bondi_radius(1.989e30, 1.0e5) > 0
      true
  """
  @spec bondi_radius(number, number) :: float
  def bondi_radius(mass, v_inf) do
    2 * @gravitational_constant * mass / (v_inf * v_inf)
  end

  # ---------------------------------------------------------------------------
  # Gravitational Time Dilation
  # ---------------------------------------------------------------------------

  @doc """
  Gravitational time dilation relative to an observer at infinity.

  dt_local / dt_inf = √(1 − r_s / r)

  Approaches 0 at the event horizon and 1 far from the black hole.

  ## Parameters
    - mass:   Mass in kilograms
    - radius: Radial coordinate r in meters (must be > r_s)

  ## Returns
    Time dilation factor in [0, 1]

  ## Examples
      iex> BlackHole.gravitational_time_dilation(1.989e30, 6.957e8) |> Float.round(6)
      0.999999
  """
  @spec gravitational_time_dilation(number, number) :: float
  def gravitational_time_dilation(mass, radius) do
    r_s = schwarzschild_radius(mass)
    :math.sqrt(max(1 - r_s / radius, 0.0))
  end
end
