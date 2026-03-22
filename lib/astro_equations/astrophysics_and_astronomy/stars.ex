defmodule AstroEquations.AstrophysicsAndAstronomy.Stars do
  @moduledoc """
  A collection of functions for calculating stellar structure and evolution properties.

  Covers:
  - Fundamental equations of stellar structure (hydrostatic equilibrium, mass conservation,
    energy generation, radiative and convective transport)
  - Timescales (Kelvin-Helmholtz, nuclear, dynamical)
  - Gravitational potential energy
  - Eddington luminosity, mass, and mass-loss rate
  - Reimers mass-loss formula (stellar winds)
  - Mass-luminosity relationship
  - Stefan-Boltzmann luminosity and effective temperature
  - Stellar radius from luminosity and temperature
  - Main-sequence lifetime
  - Wien's displacement law and Planck function
  - Jeans mass and radius
  - Chandrasekhar mass limit
  - Convection criterion (Schwarzschild)
  - Neutron star radius estimate
  - Pulsar spin-down luminosity
  - Electron-scattering opacity

  Solar units (M☉, L☉, R☉) are used where stated.
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Stellar mass in solar masses (M☉). Must be positive."
  @type solar_mass_unit :: float()

  @typedoc "Stellar luminosity in solar luminosities (L☉). Must be positive."
  @type solar_lum_unit :: float()

  @typedoc "Temperature in kelvin (K). Must be positive."
  @type temperature :: float()

  @typedoc "Stellar radius in solar radii (R☉). Must be positive."
  @type solar_radius_unit :: float()

  @typedoc "Timescale in years (yr) or seconds (s) — see function doc."
  @type timescale :: float()

  # m³ kg⁻¹ s⁻²
  @gravitational_constant 6.674_30e-11
  # J m⁻³ K⁻⁴
  @radiation_constant 7.565_723e-16
  # m/s
  @speed_of_light 2.997_924_58e8
  # kg
  @solar_mass 1.989e30
  # m
  @solar_radius 6.957e8
  # W
  @solar_luminosity 3.828e26
  # kg
  @proton_mass 1.672_621_9e-27
  # m²
  @thomson_cross_section 6.652_458_7e-29
  # J/K
  @boltzmann 1.380_649e-23
  # W m⁻² K⁻⁴
  @stefan_boltzmann 5.670_374_419e-8
  # m·K
  @wien_b 2.897_771_955e-3
  # typical accretion efficiency
  @efficiency 0.1
  # s/yr
  @seconds_per_year 3.155_76e7

  # ---------------------------------------------------------------------------
  # Equations of Stellar Structure
  # ---------------------------------------------------------------------------

  @doc """
  Pressure gradient from hydrostatic equilibrium: dP/dr = −G M_r ρ / r²

  ## Parameters
    - mass_interior: M_r — mass interior to radius r (kg)
    - density:       ρ(r) — local density (kg/m³)
    - radius:        r — radial coordinate (m)

  ## Returns
    dP/dr in Pa/m

  ## Examples
      iex> Stars.hydrostatic_equilibrium(1.989e30, 1.408e3, 6.957e8) |> Float.round(0)
      -34446.0
  """

  @spec hydrostatic_equilibrium(number, number, number) :: float
  def hydrostatic_equilibrium(mass_interior, density, radius)
      when is_positive(mass_interior) and is_positive(density) and is_positive(radius) do
    -@gravitational_constant * mass_interior * density / :math.pow(radius, 2)
  end

  @doc """
  Mass gradient from continuity (mass conservation): dM/dr = 4π r² ρ

  ## Returns
    dM/dr in kg/m

  ## Examples
      iex> Stars.mass_conservation(1408.0, 6.957e8) > 0
      true
  """

  @spec mass_conservation(number, number) :: float
  def mass_conservation(density, radius) when is_positive(density) and is_positive(radius) do
    4 * :math.pi() * :math.pow(radius, 2) * density
  end

  @doc """
  Luminosity gradient from energy generation: dL/dr = 4π r² ρ ε

  ## Parameters
    - density:    ρ(r) (kg/m³)
    - energy_gen: ε — specific energy generation rate (W/kg)
    - radius:     r (m)

  ## Returns
    dL/dr in W/m

  ## Examples
      iex> Stars.energy_equation(1408.0, 1.934e-7, 6.957e8) > 0
      true
  """

  @spec energy_equation(number, number, number) :: float
  def energy_equation(density, energy_gen, radius) when is_positive(density) do
    4 * :math.pi() * :math.pow(radius, 2) * density * energy_gen
  end

  @doc """
  Temperature gradient for radiative energy transport.

  dT/dr = −(3 κ ρ L_r) / (16 a c π T³ r²)

  ## Parameters
    - opacity:     κ (m²/kg)
    - density:     ρ (kg/m³)
    - temperature: T (K)
    - luminosity:  L_r (W)
    - radius:      r (m)

  ## Returns
    dT/dr in K/m

  ## Examples
      iex> Stars.radiative_transport(0.2, 1408, 5778, 3.828e26, 6.957e8) < 0
      true
  """

  @spec radiative_transport(number, number, number, number, number) :: float
  def radiative_transport(opacity, density, temperature, luminosity, radius)
      when is_positive(density) do
    numerator = 3 * opacity * density * luminosity

    denominator =
      16 * @radiation_constant * @speed_of_light *
        :math.pi() * :math.pow(temperature, 3) * :math.pow(radius, 2)

    -numerator / denominator
  end

  @doc """
  Adiabatic temperature gradient for convective transport (ideal ionised plasma).

  dT/dr|_ad = −(γ−1)/γ × (μ m_H / k_B) × (G M_r / r²)

  Approximated for a fully ionised hydrogen-helium plasma: γ ≈ 5/3, μ ≈ 0.62.

  ## Parameters
    - mass_interior: M_r (kg)
    - radius:        r (m)
    - mu:            Mean molecular weight (default: 0.62)
    - gamma:         Adiabatic index (default: 5/3)

  ## Returns
    Adiabatic temperature gradient in K/m

  ## Examples
      iex> Stars.adiabatic_gradient(1.989e30, 6.957e8) < 0
      true
  """

  @spec adiabatic_gradient(number, number, number, number) :: float
  def adiabatic_gradient(mass_interior, radius, mu \\ 0.62, gamma \\ 5 / 3)
      when is_positive(radius) do
    m_h = 1.6735575e-27

    -(gamma - 1) / gamma * (mu * m_h / @boltzmann) *
      (@gravitational_constant * mass_interior / :math.pow(radius, 2))
  end

  @doc """
  Schwarzschild convection criterion: checks if a layer is convectively unstable.

  A layer is unstable to convection when:
  |dT/dr|_actual > |dT/dr|_adiabatic

  ## Parameters
    - grad_actual:    Actual local temperature gradient |dT/dr| (K/m, positive magnitude)
    - grad_adiabatic: Adiabatic gradient magnitude (K/m, positive magnitude)

  ## Returns
    true if convectively unstable (Schwarzschild criterion satisfied)

  ## Examples
      iex> Stars.convectively_unstable?(2.0e-4, 1.0e-4)
      true
  """

  @spec convectively_unstable?(number, number) :: boolean
  def convectively_unstable?(grad_actual, grad_adiabatic) do
    grad_actual > grad_adiabatic
  end

  @doc """
  Rosseland mean opacity for electron scattering (fully ionised plasma).

  κ_es = 0.2 × (1 + X) m²/kg

  where X is the hydrogen mass fraction (X ≈ 0.70 for solar composition).

  ## Parameters
    - hydrogen_fraction: Mass fraction of hydrogen X (default: 0.70)

  ## Returns
    Electron-scattering opacity in m²/kg

  ## Examples
      iex> Stars.opacity_electron_scattering() |> Float.round(4)
      0.034
  """

  @spec opacity_electron_scattering(number) :: float
  def opacity_electron_scattering(hydrogen_fraction \\ 0.70) do
    0.2 * (1 + hydrogen_fraction)
  end

  # ---------------------------------------------------------------------------
  # Stellar Timescales
  # ---------------------------------------------------------------------------

  @doc """
  Kelvin-Helmholtz (thermal) timescale: τ_KH = G M² / (R L)

  The time for a star to radiate away its gravitational potential energy.

  ## Parameters
    - mass:       In solar masses
    - radius:     In solar radii
    - luminosity: In solar luminosities

  ## Returns
    Timescale in years

  ## Examples
      iex> Stars.kelvin_helmholtz_timescale(1, 1, 1) |> round()
      31_484_441
  """

  @spec kelvin_helmholtz_timescale(number, number, number) :: float
  def kelvin_helmholtz_timescale(mass, radius, luminosity) when is_positive(mass) do
    m = mass * @solar_mass
    r = radius * @solar_radius
    l = luminosity * @solar_luminosity

    @gravitational_constant * :math.pow(m, 2) / (r * l) / @seconds_per_year
  end

  @doc """
  Nuclear (main-sequence) timescale: τ_nuc ≈ 10¹⁰ × M⁻³ yr

  Uses the mass-luminosity relation L ∝ M⁴ for solar-type stars.
  For M = 1 M☉ this gives ~10¹⁰ yr, consistent with the Sun's expected
  main-sequence lifetime.

  ## Parameters
    - mass: Stellar mass in solar masses

  ## Returns
    Timescale in years

  ## Examples
      iex> Stars.nuclear_timescale(1.0) |> Float.round(-9)
      1.0e10
  """

  @spec nuclear_timescale(number) :: float
  def nuclear_timescale(mass) when is_positive(mass), do: 1.0e10 * :math.pow(mass, -3)

  @doc """
  Dynamical (free-fall) timescale: τ_dyn = √(R³ / (G M))

  ## Parameters
    - mass:   In solar masses
    - radius: In solar radii

  ## Returns
    Timescale in seconds

  ## Examples
      iex> Stars.dynamical_timescale(1.0, 1.0) > 0
      true
  """

  @spec dynamical_timescale(number, number) :: float
  def dynamical_timescale(mass, radius) when is_positive(mass) do
    m = mass * @solar_mass
    r = radius * @solar_radius
    :math.sqrt(:math.pow(r, 3) / (@gravitational_constant * m))
  end

  @doc """
  Main-sequence lifetime from mass-luminosity scaling:
  t_MS ≈ t_☉ × (M/M☉) / (L/L☉)

  ## Parameters
    - mass:       Stellar mass in solar masses
    - t_sun_gyr:  Solar main-sequence lifetime in Gyr (default: 10.0)

  ## Returns
    Main-sequence lifetime in Gyr

  ## Examples
      iex> Stars.main_sequence_lifetime(1.0) |> Float.round(4)
      10.0
  """

  @spec main_sequence_lifetime(number, number) :: float
  def main_sequence_lifetime(mass, t_sun_gyr \\ 10.0) when is_positive(mass) do
    luminosity = mass_luminosity(mass)
    t_sun_gyr * mass / luminosity
  end

  # ---------------------------------------------------------------------------
  # Energetics
  # ---------------------------------------------------------------------------

  @doc """
  Gravitational potential energy: U = −G M² / R

  ## Returns
    Potential energy in joules

  ## Examples
      iex> Stars.gravitational_potential_energy(1.0, 1.0) < 0
      true
  """

  @spec gravitational_potential_energy(number, number) :: float
  def gravitational_potential_energy(mass, radius) when is_positive(mass) do
    m = mass * @solar_mass
    r = radius * @solar_radius
    -@gravitational_constant * :math.pow(m, 2) / r
  end

  # ---------------------------------------------------------------------------
  # Eddington Limits
  # ---------------------------------------------------------------------------

  @doc """
  Eddington luminosity limit: L_Edd = 4π G M m_p c / σ_T

  ## Parameters
    - mass: Stellar mass in solar masses

  ## Returns
    Luminosity limit in solar luminosities

  ## Examples
      iex> Stars.eddington_luminosity(1) |> round()
      32_000
  """

  @spec eddington_luminosity(number) :: float
  def eddington_luminosity(mass) when is_positive(mass) do
    4 * :math.pi() * @gravitational_constant * mass * @solar_mass *
      @proton_mass * @speed_of_light / @thomson_cross_section / @solar_luminosity
  end

  @doc """
  Eddington mass from luminosity (approximate inverse): M ≈ 3.1×10⁻⁵ × L/L☉

  ## Examples
      iex> Stars.eddington_mass(32_000.0) |> Float.round(2)
      1.0
  """

  @spec eddington_mass(number) :: float
  def eddington_mass(luminosity), do: 3.1e-5 * luminosity

  @doc """
  Eddington mass-loss rate: ṁ = L_Edd / (η c²)

  ## Returns
    Mass-loss rate in M☉/yr

  ## Examples
      iex> Stars.eddington_mass_loss_rate(1.0) > 0
      true
  """

  @spec eddington_mass_loss_rate(number) :: float
  def eddington_mass_loss_rate(mass) when is_positive(mass) do
    l_edd = eddington_luminosity(mass) * @solar_luminosity

    l_edd / (@efficiency * :math.pow(@speed_of_light, 2)) *
      @seconds_per_year / @solar_mass
  end

  @doc """
  Reimers mass-loss formula for cool giant stars.

  ṁ = η_R × L R / (M)  (Reimers 1975 empirical form)
  ṁ [M☉/yr] = 4×10⁻¹³ × η_R × (L/L☉)(R/R☉)/(M/M☉)

  Typical η_R ≈ 0.4–0.5 for RGB/AGB stars.

  ## Parameters
    - mass:      Stellar mass in solar masses
    - radius:    Stellar radius in solar radii
    - luminosity: Stellar luminosity in solar luminosities
    - eta_r:     Reimers efficiency parameter (default: 0.5)

  ## Returns
    Mass-loss rate in M☉/yr

  ## Examples
      iex> Stars.stellar_wind_mass_loss(1.0, 100.0, 1000.0) > 0
      true
  """

  @spec stellar_wind_mass_loss(number, number, number, number) :: float
  def stellar_wind_mass_loss(mass, radius, luminosity, eta_r \\ 0.5) when is_positive(mass) do
    4.0e-13 * eta_r * luminosity * radius / mass
  end

  # ---------------------------------------------------------------------------
  # Stellar Luminosity & Radius
  # ---------------------------------------------------------------------------

  @doc """
  Luminosity from Stefan-Boltzmann law: L = 4π R² σ T_eff⁴

  ## Returns
    Luminosity in Watts

  ## Examples
      iex> Stars.stefan_boltzmann_luminosity(6.957e8, 5778) > 0
      true
  """

  @spec stefan_boltzmann_luminosity(number, number) :: float
  def stefan_boltzmann_luminosity(radius, temperature)
      when is_positive(radius) and is_positive(temperature) do
    4 * :math.pi() * :math.pow(radius, 2) *
      @stefan_boltzmann * :math.pow(temperature, 4)
  end

  @doc """
  Stellar radius from luminosity and effective temperature: R = √(L / (4π σ T⁴))

  ## Examples
      iex> Stars.radius_from_luminosity_temperature(3.828e26, 5778) > 0
      true
  """

  @spec radius_from_luminosity_temperature(number, number) :: float
  def radius_from_luminosity_temperature(luminosity, temperature) when is_positive(temperature) do
    :math.sqrt(
      luminosity /
        (4 * :math.pi() * @stefan_boltzmann * :math.pow(temperature, 4))
    )
  end

  @doc """
  Mass-luminosity relation (piecewise power-law fit).

  - M < 0.43 M☉  : L = 0.23 M^2.3
  - 0.43–2 M☉   : L = M^4
  - 2–20 M☉     : L = 1.5 M^3.5
  - 20–55 M☉    : smooth interpolation from 20 M☉ reference
  - M > 55 M☉   : L ≈ 32_000 M  (Eddington-limited)

  ## Parameters
    - mass: Stellar mass in solar masses

  ## Returns
    Luminosity in solar luminosities

  ## Examples
      iex> Stars.mass_luminosity(1) |> Float.round(4)
      1.0
      iex> Stars.mass_luminosity(10) |> Float.round(4)
      3162.2777
  """

  @spec mass_luminosity(number) :: float
  def mass_luminosity(mass) when is_positive(mass) do
    cond do
      mass < 0.43 -> 0.23 * :math.pow(mass, 2.3)
      mass < 2 -> 1.0 * :math.pow(mass, 4.0)
      mass < 20 -> 1.5 * :math.pow(mass, 3.5)
      mass < 55 -> 1.5 * :math.pow(20, 3.5) * :math.pow(mass / 20, 2.3)
      true -> 32_000.0 * mass
    end
  end

  # ---------------------------------------------------------------------------
  # Radiation
  # ---------------------------------------------------------------------------

  @doc """
  Wien's displacement law: λ_max = b / T

  ## Returns
    Peak wavelength in meters

  ## Examples
      iex> Stars.wien_peak_wavelength(5778) |> Float.round(11)
      5.015e-7
  """

  @spec wien_peak_wavelength(number) :: float
  def wien_peak_wavelength(temperature) when is_positive(temperature), do: @wien_b / temperature

  @doc """
  Planck spectral radiance B_λ(T) in W m⁻² sr⁻¹ m⁻¹.

  ## Examples
      iex> Stars.planck_function(500.0e-9, 5778) > 0
      true
  """

  @spec planck_function(number, number) :: float
  def planck_function(wavelength, temperature) when is_positive(wavelength) do
    h = 6.626_070_15e-34
    c = @speed_of_light
    kb = @boltzmann

    2 * h * c * c /
      (:math.pow(wavelength, 5) * (:math.exp(h * c / (wavelength * kb * temperature)) - 1))
  end

  # ---------------------------------------------------------------------------
  # Star Formation (Jeans Criterion)
  # ---------------------------------------------------------------------------

  @doc """
  Jeans mass: M_J = (5 k_B T / (G m_H μ))^(3/2) × (3 / (4π ρ))^(1/2)

  ## Parameters
    - temperature: Cloud temperature in Kelvin
    - density:     Mass density in kg/m³
    - mu:          Mean molecular weight (default: 2.0 for H₂)

  ## Returns
    Jeans mass in kg

  ## Examples
      iex> Stars.jeans_mass(10, 1.0e-17) > 0
      true
  """

  @spec jeans_mass(number, number, number) :: float
  def jeans_mass(temperature, density, mu \\ 2.0)
      when is_positive(temperature) and is_positive(density) do
    m_h = 1.6735575e-27

    term1 =
      :math.pow(
        5 * @boltzmann * temperature /
          (@gravitational_constant * m_h * mu),
        3 / 2
      )

    term2 = :math.sqrt(3 / (4 * :math.pi() * density))
    term1 * term2
  end

  @doc """
  Jeans radius: R_J = √(15 k_B T / (4π G m_H μ ρ))

  ## Examples
      iex> Stars.jeans_radius(10, 1.0e-17) > 0
      true
  """

  @spec jeans_radius(number, number, number) :: float
  def jeans_radius(temperature, density, mu \\ 2.0) when is_positive(temperature) do
    m_h = 1.6735575e-27

    :math.sqrt(
      15 * @boltzmann * temperature /
        (4 * :math.pi() * @gravitational_constant * m_h * mu * density)
    )
  end

  # ---------------------------------------------------------------------------
  # Compact Remnants
  # ---------------------------------------------------------------------------

  @doc """
  Chandrasekhar mass limit for a white dwarf: M_Ch ≈ 5.83 / μ_e² M☉

  ## Parameters
    - mu_e: Mean molecular weight per electron (default: 2.0 for C/O WD)

  ## Returns
    Chandrasekhar mass in solar masses

  ## Examples
      iex> Stars.chandrasekhar_mass() |> Float.round(4)
      1.4575
  """

  @spec chandrasekhar_mass(number) :: float
  def chandrasekhar_mass(mu_e \\ 2.0), do: 5.83 / :math.pow(mu_e, 2)

  @doc """
  Approximate neutron star (NS) radius from the nuclear saturation density.

  R_NS ≈ (3 M / (4π ρ_nuc))^(1/3)

  Uses ρ_nuc ≈ 2.3×10¹⁷ kg/m³ (≈ 4× nuclear saturation density).

  ## Parameters
    - mass: Neutron star mass in solar masses (default: 1.4)

  ## Returns
    Approximate NS radius in meters

  ## Examples
      iex> Stars.neutron_star_radius() > 0
      true
  """

  @spec neutron_star_radius(number) :: float
  def neutron_star_radius(mass \\ 1.4) when is_positive(mass) do
    rho_nuc = 2.3e17
    m_kg = mass * @solar_mass
    :math.pow(3 * m_kg / (4 * :math.pi() * rho_nuc), 1 / 3)
  end

  @doc """
  Pulsar spin-down luminosity from period and period derivative.

  L_sd = −I Ω Ω̇ = 4π² I Ṗ / P³

  Assuming a canonical neutron star moment of inertia I = 10⁴⁵ g·cm² = 10³⁸ kg·m².

  ## Parameters
    - period:          Spin period P in seconds
    - period_deriv:    Period derivative Ṗ (dimensionless, s/s)
    - moment_inertia:  I in kg·m² (default: 1.0×10³⁸)

  ## Returns
    Spin-down luminosity in Watts

  ## Examples
      iex> Stars.pulsar_spindown_luminosity(0.033, 4.2e-13) > 0
      true
  """

  @spec pulsar_spindown_luminosity(number, number, number) :: float
  def pulsar_spindown_luminosity(period, period_deriv, moment_inertia \\ 1.0e38) do
    4 * :math.pow(:math.pi(), 2) * moment_inertia * period_deriv / :math.pow(period, 3)
  end

  @doc """
  Characteristic pulsar age (spin-down age): τ_c = P / (2 Ṗ)

  Assumes a braking index of n = 3 (magnetic dipole radiation) and
  initial spin period P₀ ≪ P.

  ## Parameters
    - period:       Spin period (s)
    - period_deriv: Period derivative (s/s)

  ## Returns
    Characteristic age in seconds

  ## Examples
      iex> Stars.pulsar_characteristic_age(0.033, 4.2e-13) > 0
      true
  """

  @spec pulsar_characteristic_age(number, number) :: float

  def pulsar_characteristic_age(period, period_deriv) do
    period / (2 * period_deriv)
  end
end
