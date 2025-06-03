defmodule AstroEquations.Physics.Oscillations do
  @moduledoc """
  Provides functions for calculating spring-related physics properties.

  This module implements Hooke's Law for springs, potential energy in springs,
  and angular frequency of a spring-mass system.
  """

  @doc """
  Calculates the force exerted by a spring using Hooke's Law.

  ## Parameters
    - `k_s`: Spring constant (N/m)
    - `x`: Displacement from equilibrium position (m)

  ## Examples
      iex> Physics.Spring.force(10, 0.5)
      -5.0
  """
  @spec force(number, number) :: float
  def force(k_s, x), do: -k_s * x

  @doc """
  Calculates the potential energy stored in a spring.

  ## Parameters
    - `k_s`: Spring constant (N/m)
    - `x`: Displacement from equilibrium position (m)

  ## Examples
      iex> Physics.Spring.potential_energy(10, 0.5)
      1.25
  """
  @spec potential_energy(number, number) :: float
  def potential_energy(k_s, x), do: 0.5 * k_s * :math.pow(x, 2)

  @doc """
  Calculates the angular frequency of a spring-mass system.

  ## Parameters
    - `k_s`: Spring constant (N/m)
    - `m`: Mass of the object (kg)

  ## Examples
      iex> Physics.Spring.angular_frequency(10, 2.5)
      2.0
  """
  @spec angular_frequency(number, number) :: float
  def angular_frequency(k_s, m), do: :math.sqrt(k_s / m)
end
