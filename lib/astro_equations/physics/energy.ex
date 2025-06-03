defmodule AstroEquations.Physics.Energy do
  @moduledoc """
  A collection of fundamental physics formulas and calculations.

  This module provides functions for calculating work, buoyancy, friction,
  wave properties, and related physics concepts.
  """

  @doc """
  Calculates the work done by a force along a path.

  ## Parameters
    - force_vector: The force vector [Fx, Fy, Fz]
    - displacement_vector: The displacement vector [dx, dy, dz]

  ## Examples
      iex> Physics.work([2, 0, 0], [3, 0, 0])
      6.0

      iex> Physics.work([1, 2, 3], [4, 5, 6])
      32.0
  """
  @spec work(list(number), list(number)) :: float
  def work(force_vector, displacement_vector) do
    Enum.zip(force_vector, displacement_vector)
    |> Enum.map(fn {f, s} -> f * s end)
    |> Enum.sum()
    |> Kernel./(1)
  end
end
