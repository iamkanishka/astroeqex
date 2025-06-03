defmodule AstroEquations.Physics.Materials do
  @moduledoc """
  Provides functions for calculating material properties.

  This module implements density calculations for both discrete and continuous mass distributions.
  """

  @doc """
  Calculates density as mass per unit volume (discrete version).

  ## Parameters
    - `mass`: Total mass of the object (kg)
    - `volume`: Total volume of the object (m³)

  ## Examples
      iex> Physics.Materials.density(10, 2)
      5.0

      iex> Physics.Materials.density(5.5, 2.2)
      2.5
  """
  @spec density(number, number) :: float
  def density(mass, volume) when volume != 0, do: mass / volume
  def density(_mass, 0), do: raise(ArgumentError, "Volume cannot be zero")

  @doc """
  Calculates density as the derivative of mass with respect to volume (continuous version).

  ## Parameters
    - `mass_function`: Function that describes how mass changes with volume
    - `volume`: Volume at which to calculate the density (m³)
    - `delta`: Small change in volume for numerical differentiation (default: 1.0e-6)

  ## Examples
      iex> linear_mass = fn v -> 2.5 * v end
      iex> Physics.Materials.continuous_density(linear_mass, 2.0)
      2.5
  """
  @spec continuous_density((number -> number), number, number) :: float
  def continuous_density(mass_function, volume, delta \\ 1.0e-6) do
    dm = mass_function.(volume + delta) - mass_function.(volume - delta)
    dV = 2 * delta
    dm / dV
  end
end
