defmodule AstroEquations.AstrophysicsAndAstronomy.BlackHole do
  @moduledoc """
  A collection of functions for calculating black hole properties.

  This module provides implementations of fundamental black hole equations,
  starting with the Schwarzschild radius calculation.
  """

  @gravitational_constant 6.67430e-11  # m^3 kg^-1 s^-2
  @speed_of_light 2.99792458e8        # m/s

  @doc """
  Calculates the Schwarzschild radius (event horizon) for a given mass.

  ## Parameters
    - mass: The mass of the object in kilograms

  ## Returns
    The Schwarzschild radius in meters

  ## Examples
      iex> BlackHolePhysics.schwarzschild_radius(1.989e30) # For the Sun
      2953.2500761002446
  """
  @spec schwarzschild_radius(number) :: float
  def schwarzschild_radius(mass) do
    2 * @gravitational_constant * mass / :math.pow(@speed_of_light, 2)
  end

  @doc """
  Calculates the Schwarzschild radius in solar radii units.

  ## Parameters
    - solar_masses: Mass in units of solar masses (M☉)

  ## Returns
    The Schwarzschild radius in meters

  ## Examples
      iex> BlackHolePhysics.schwarzschild_radius_solar(1) # For 1 solar mass
      2953.2500761002446
  """
  @spec schwarzschild_radius_solar(number) :: float
  def schwarzschild_radius_solar(solar_masses) do
    schwarzschild_radius(solar_masses * 1.989e30)
  end
end
