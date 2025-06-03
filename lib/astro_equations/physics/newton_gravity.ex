defmodule AstroEquations.Physics.NewtonGravity do
  @moduledoc """
  Provides functions for calculating Newtonian gravitational interactions.

  This module implements gravitational force, potential, field, potential energy,
  and Kepler's Third Law calculations.
  """
  # m³ kg⁻¹ s⁻²
  @gravitational_constant 6.67430e-11

  @doc """
  Calculates the magnitude of gravitational force between two masses.

  ## Parameters
    - `m`: First mass (kg)
    - `M`: Second mass (kg)
    - `r`: Distance between masses (m)

  ## Examples
      iex> AstroEquations.Physics.NewtonGravity.force(5.972e24, 7.348e22, 3.844e8) # Earth-Moon system
      1.9817963747937598e20
  """
  @spec force(number, number, number) :: float
  def force(first_mass, second_mass, r)
      when is_number(first_mass) and is_number(second_mass) and is_number(r) and r > 0 do
    @gravitational_constant * first_mass * second_mass / :math.pow(r, 2)
  end

  def force(_m, _M, _r), do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Calculates gravitational potential at a point due to a mass.

  ## Parameters
    - `M`: Source mass (kg)
    - `r`: Distance from mass (m)

  ## Examples
      iex> Physics.Gravity.potential(5.972e24, 6.371e6) # Earth's surface
      -62565145.91113584
  """
  @spec potential(number, number) :: float
  def potential(source_mass, r) when r > 0, do: -@gravitational_constant * source_mass / r
  def potential(_source_mass, _r), do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Calculates gravitational field strength at a distance from a mass.

  ## Parameters
    - `M`: Source mass (kg)
    - `r`: Distance from mass (m)

  ## Examples
      iex> Physics.Gravity.field(5.972e24, 6.371e6) # Earth's surface
      9.819649653686391
  """
  @spec field(number, number) :: float
  def field(source_mass, r) when r > 0,
    do: @gravitational_constant * source_mass / :math.pow(r, 2)

  def field(_source_mass, _r), do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Calculates gravitational potential energy between two masses.

  ## Parameters
    - `m`: First mass (kg)
    - `M`: Second mass (kg)
    - `r`: Distance between masses (m)

  ## Examples
      iex> Physics.Gravity.potential_energy(1000, 5.972e24, 6.371e6) # 1 ton on Earth's surface
      -6.256514591113584e10
  """
  @spec potential_energy(number, number, number) :: float
  def potential_energy(first_mass, second_mass, r) when r > 0,
    do: -@gravitational_constant * first_mass * second_mass / r

  def potential_energy(_first_mass, _second_mass, _r),
    do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Approximates gravitational potential energy near a planet's surface (mgh).

  ## Parameters
    - `m`: Object mass (kg)
    - `g`: Gravitational acceleration (m/s²)
    - `h`: Height above surface (m)

  ## Examples
      iex> Physics.Gravity.approximate_potential_energy(1000, 9.81, 100)
      981000.0
  """
  @spec approximate_potential_energy(number, number, number) :: float
  def approximate_potential_energy(m, g, h), do: m * g * h

  @doc """
  Calculates the orbital period according to Kepler's Third Law.

  ## Parameters
    - `r`: Semi-major axis (m)
    - `m`: Primary mass (kg)
    - `M`: Secondary mass (kg)

  ## Examples
      iex> Physics.Gravity.orbital_period(3.844e8, 5.972e24, 7.348e22) # Earth-Moon system
      2360449.1175632763
  """
  @spec orbital_period(number, number, number) :: float
  def orbital_period(r, primary_mass, secondary_mass) when r > 0 do
    numerator = 4 * :math.pow(:math.pi(), 2) * :math.pow(r, 3)
    denominator = @gravitational_constant * (primary_mass + secondary_mass)
    :math.sqrt(numerator / denominator)
  end

  def orbital_period(_r, _m, _M), do: raise(ArgumentError, "Distance must be positive")
end
