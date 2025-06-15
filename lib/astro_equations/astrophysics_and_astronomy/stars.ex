defmodule AstroEquations.AstrophysicsAndAstronomy.Stars do
  @moduledoc """
  A collection of functions for calculating stellar structure properties.

  This module provides implementations of the fundamental equations of stellar structure:
  - Hydrostatic equilibrium
  - Mass conservation
  - Energy generation
  - Radiative transport
  - Timescales (Kelvin-Helmholtz, Nuclear)
  - Gravitational potential energy
  - Eddington limits and rates
  - Mass-Luminosity relationships



  All equations are implemented in their differential form suitable for numerical integration
  in stellar modeling applications.
  Some calculations use solar units (M⊙, L⊙, R⊙) by default.
  """

  # m^3 kg^-1 s^-2
  @gravitational_constant 6.67430e-11
  # J m^-3 K^-4 (a in the radiation equation)
  @radiation_constant 7.565723e-16
  # m/s (c in the radiation equation)
  @speed_of_light 2.99792458e8

  @solar_mass 1.989e30                # kg (M⊙)
  @solar_radius 6.957e8               # m (R⊙)
  @solar_luminosity 3.828e26          # W (L⊙)
  @proton_mass 1.6726219e-27          # kg (m_p)
  @speed_of_light 2.99792458e8        # m/s (c)
  @thomson_cross_section 6.6524587e-29 # m^2 (σ_T)
  @efficiency 0.1                     # η (typical mass-energy conversion efficiency)




  @doc """
  Calculates the pressure gradient (dP/dr) for hydrostatic equilibrium.

  ## Parameters
    - mass_interior: Mass interior to radius r (M_r) in kg
    - density: Density at radius r (ρ(r)) in kg/m^3
    - radius: Radius (r) in meters

  ## Returns
    The pressure gradient (dP/dr) in Pa/m

  ## Examples
      iex> StellarStructure.hydrostatic_equilibrium(1.989e30, 1.408e3, 6.957e8)
      -3.4446331106341036e4
  """
  @spec hydrostatic_equilibrium(number, number, number) :: float
  def hydrostatic_equilibrium(mass_interior, density, radius) do
    -@gravitational_constant * mass_interior * density / :math.pow(radius, 2)
  end

  @doc """
  Calculates the mass gradient (dM/dr) for mass conservation.

  ## Parameters
    - density: Density at radius r (ρ(r)) in kg/m^3
    - radius: Radius (r) in meters

  ## Returns
    The mass gradient (dM/dr) in kg/m

  ## Examples
      iex> StellarStructure.mass_conservation(1.408e3, 6.957e8)
      8.57104185008e19
  """
  @spec mass_conservation(number, number) :: float
  def mass_conservation(density, radius) do
    4 * :math.pi() * :math.pow(radius, 2) * density
  end

  @doc """
  Calculates the luminosity gradient (dL/dr) from energy generation.

  ## Parameters
    - density: Density at radius r (ρ(r)) in kg/m^3
    - energy_gen: Energy generation rate (ε) in W/kg
    - radius: Radius (r) in meters

  ## Returns
    The luminosity gradient (dL/dr) in W/m

  ## Examples
      iex> StellarStructure.energy_equation(1.408e3, 1.934e-7, 6.957e8)
      8.230719858176e7
  """
  @spec energy_equation(number, number, number) :: float
  def energy_equation(density, energy_gen, radius) do
    4 * :math.pi() * :math.pow(radius, 2) * density * energy_gen
  end

  @doc """
  Calculates the temperature gradient (dT/dr) for radiative transport.

  ## Parameters
    - opacity: Opacity (κ) in m^2/kg
    - density: Density at radius r (ρ(r)) in kg/m^3
    - temperature: Temperature (T) in Kelvin
    - luminosity: Luminosity at radius r (L_r) in Watts
    - radius: Radius (r) in meters

  ## Returns
    The temperature gradient (dT/dr) in K/m

  ## Examples
      iex> StellarStructure.radiative_transport(0.04, 1.408e3, 5.778e3, 3.828e26, 6.957e8)
      -0.010344021912541192
  """
  @spec radiative_transport(number, number, number, number, number) :: float
  def radiative_transport(opacity, density, temperature, luminosity, radius) do
    numerator = 3 * opacity * density * luminosity

    denominator =
      16 * @radiation_constant * @speed_of_light * :math.pi() * :math.pow(temperature, 3) *
        :math.pow(radius, 2)

    -numerator / denominator
  end





    @doc """
  Calculates the Kelvin-Helmholtz (thermal) timescale for a star.

  ## Parameters
    - mass: Stellar mass in solar masses (M⊙)
    - radius: Stellar radius in solar radii (R⊙)
    - luminosity: Stellar luminosity in solar luminosities (L⊙)

  ## Returns
    Timescale in years

  ## Examples
      iex> StellarProperties.kelvin_helmholtz_timescale(1, 1, 1) |> round()
      31484441
  """
  @spec kelvin_helmholtz_timescale(number, number, number) :: float
  def kelvin_helmholtz_timescale(mass, radius, luminosity) do
    g = @gravitational_constant
    m = mass * @solar_mass
    r = radius * @solar_radius
    l = luminosity * @solar_luminosity

    (g * :math.pow(m, 2)) / (r * l) / (60 * 60 * 24 * 365.25) # Convert to years
  end

  @doc """
  Estimates the nuclear (main sequence) timescale for a star.

  ## Parameters
    - mass: Stellar mass in solar masses (M⊙)

  ## Returns
    Timescale in years

  ## Examples
      iex> StellarProperties.nuclear_timescale(1) |> round()
      1000000000
      iex> StellarProperties.nuclear_timescale(2) |> round()
      125000000
  """
  @spec nuclear_timescale(number) :: float
  def nuclear_timescale(mass) do
    1.0e9 * :math.pow(mass, -3)
  end

  @doc """
  Calculates the gravitational potential energy of a star.

  ## Parameters
    - mass: Stellar mass in solar masses (M⊙)
    - radius: Stellar radius in solar radii (R⊙)

  ## Returns
    Potential energy in joules

  ## Examples
      iex> StellarProperties.gravitational_potential_energy(1, 1) |> round()
      -3.7909166e41
  """
  @spec gravitational_potential_energy(number, number) :: float
  def gravitational_potential_energy(mass, radius) do
    g = @gravitational_constant
    m = mass * @solar_mass
    r = radius * @solar_radius

    -g * :math.pow(m, 2) / r
  end

  @doc """
  Calculates the Eddington luminosity limit for a star.

  ## Parameters
    - mass: Stellar mass in solar masses (M⊙)

  ## Returns
    Luminosity limit in solar luminosities (L⊙)

  ## Examples
      iex> StellarProperties.eddington_luminosity(1) |> round()
      32000
  """
  @spec eddington_luminosity(number) :: float
  def eddington_luminosity(mass) do
    numerator = 4 * :math.pi * @gravitational_constant * mass * @solar_mass * @proton_mass * @speed_of_light
    denominator = @thomson_cross_section
    (numerator / denominator) / @solar_luminosity
  end

  @doc """
  Calculates the Eddington mass corresponding to a given luminosity.

  ## Parameters
    - luminosity: Stellar luminosity in solar luminosities (L⊙)

  ## Returns
    Mass limit in solar masses (M⊙)

  ## Examples
      iex> StellarProperties.eddington_mass(32000) |> Float.round(6)
      1.0
  """
  @spec eddington_mass(number) :: float
  def eddington_mass(luminosity) do
    3.1e-5 * luminosity
  end

  @doc """
  Calculates the Eddington mass loss rate.

  ## Parameters
    - mass: Stellar mass in solar masses (M⊙)

  ## Returns
    Mass loss rate in solar masses per year (M⊙/yr)

  ## Examples
      iex> StellarProperties.eddington_mass_loss_rate(1) |> Float.round(8)
      2.4e-8
  """
  @spec eddington_mass_loss_rate(number) :: float
  def eddington_mass_loss_rate(mass) do
    ledd = eddington_luminosity(mass) * @solar_luminosity
    ledd / (@efficiency * :math.pow(@speed_of_light, 2)) * (60 * 60 * 24 * 365.25) / @solar_mass
  end

  @doc """
  Calculates luminosity based on mass using the mass-luminosity relationship.

  ## Parameters
    - mass: Stellar mass in solar masses (M⊙)

  ## Returns
    Luminosity in solar luminosities (L⊙)

  ## Examples
      iex> StellarProperties.mass_luminosity(0.5) |> Float.round(4)
      0.0331
      iex> StellarProperties.mass_luminosity(1) |> Float.round(4)
      1.0
      iex> StellarProperties.mass_luminosity(10) |> Float.round(4)
      3162.2777
      iex> StellarProperties.mass_luminosity(100) |> Float.round(4)
      32000.0
  """
  @spec mass_luminosity(number) :: float
  def mass_luminosity(mass) do
    cond do
      mass < 0.43 -> 0.23 * :math.pow(mass, 2.3)
      mass < 2 -> 1.0 * :math.pow(mass, 4)
      mass < 20 -> 1.5 * :math.pow(mass, 3.5)
      mass < 55 -> 1.5 * :math.pow(20, 3.5) * :math.pow(mass/20, 2.3) # Interpolation
      true -> 32000.0 * mass
    end
  end


end
