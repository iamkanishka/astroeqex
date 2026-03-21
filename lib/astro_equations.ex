defmodule AstroEquations do
  @moduledoc """
  AstroEquations — a comprehensive Elixir library of scientific and astronomical equations.

  ## Namespaces

  ### Astrophysics & Astronomy
  - `AstroEquations.AstrophysicsAndAstronomy.Astrometry`
  - `AstroEquations.AstrophysicsAndAstronomy.BlackHole`
  - `AstroEquations.AstrophysicsAndAstronomy.Galaxies`
  - `AstroEquations.AstrophysicsAndAstronomy.Instrumentation`
  - `AstroEquations.AstrophysicsAndAstronomy.Stars`

  ### Physics
  - `AstroEquations.Physics.Electromagnetism`
  - `AstroEquations.Physics.Energy`
  - `AstroEquations.Physics.Forces`
  - `AstroEquations.Physics.GeneralRelativity`
  - `AstroEquations.Physics.Materials`
  @doc \"""
  - `AstroEquations.Physics.Motion`
  @doc \"""
  - `AstroEquations.Physics.NewtonGravity`
  @doc \"""
  - `AstroEquations.Physics.Oscillations`
  @doc \"""
  - `AstroEquations.Physics.QuantumMechanics`
  @doc \"""
  - `AstroEquations.Physics.SpecialRelativity`
  @doc \"""
  - `AstroEquations.Physics.Thermodynamics`
  @doc \"""
  - `AstroEquations.Physics.Waves`

  ### Mathematics
  @doc \"""
  - `AstroEquations.Mathematics.Calculus`
  @doc \"""
  - `AstroEquations.Mathematics.Geometry`
  @doc \"""
  - `AstroEquations.Mathematics.Notation`
  @doc \"""
  - `AstroEquations.Mathematics.Trigonometry`

  ### Statistics
  @doc \"""
  - `AstroEquations.Statistics.StandardDeviation`
  @doc \"""
  - `AstroEquations.Statistics.Variance`

  @doc \"""
  All calculations use SI units unless the function documentation states otherwise.
  """

  @doc """
  Returns the current library version.

  ## Examples
      iex> AstroEquations.version()
  @doc \"""
      "0.2.0"

  ## Examples
      iex> AstroEquations.version()
      "0.2.0"
  """
  @spec version() :: String.t()
  def version, do: "0.2.0"
end
