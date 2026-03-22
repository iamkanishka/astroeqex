defmodule AstroEquations do
  @moduledoc """
  **AstroEquations** — a comprehensive Elixir library of scientific and
  astronomical equations.

  696+ pure functions across 23 modules, all in SI units with full
  `@spec`, `@doc`, `@type`, and input guard coverage.

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
  - `AstroEquations.Physics.Motion`
  - `AstroEquations.Physics.NewtonGravity`
  - `AstroEquations.Physics.Oscillations`
  - `AstroEquations.Physics.QuantumMechanics`
  - `AstroEquations.Physics.SpecialRelativity`
  - `AstroEquations.Physics.Thermodynamics`
  - `AstroEquations.Physics.Waves`

  ### Mathematics

  - `AstroEquations.Mathematics.Calculus`
  - `AstroEquations.Mathematics.Geometry`
  - `AstroEquations.Mathematics.Notation`
  - `AstroEquations.Mathematics.Trigonometry`

  ### Statistics

  - `AstroEquations.Statistics.StandardDeviation`
  - `AstroEquations.Statistics.Variance`

  ## Deprecated

  The modules below have spelling errors in their names. They delegate all
  calls to the correct module and will be removed in v1.0.

  - `AstroEquations.AstrophysicsAndAstronomy.Instrumentaion` →
    use `Instrumentation`
  - `AstroEquations.Mathematics.Trignometry` →
    use `Trigonometry`

  ## Design

  All functions are **pure mathematical transformations** with no side effects.
  Physics-critical functions carry input guards (`when is_positive(mass)`,
  `when is_non_negative(v)`) that reject physically nonsensical inputs at the
  call boundary.

  Every module exports three reusable guard predicates:

      import AstroEquations.Physics.Energy,
        only: [is_positive: 1, is_non_negative: 1, is_real: 1]

  ## Units

  All inputs and outputs use SI base units unless the function documentation
  explicitly states otherwise.
  """

  @doc since: "0.1.0"
  @doc """
  Returns the current library version string.

  ## Examples

      iex> AstroEquations.version()
      "0.3.0"
  """
  @spec version :: String.t()
  def version, do: "0.3.0"

  @doc """
  Returns a map of all live modules to their public function counts.

  ## Examples

      iex> AstroEquations.modules() |> map_size()
      23
  """
  @spec modules :: %{module() => non_neg_integer()}
  def modules do
    %{
      AstroEquations.AstrophysicsAndAstronomy.Astrometry => 32,
      AstroEquations.AstrophysicsAndAstronomy.BlackHole => 15,
      AstroEquations.AstrophysicsAndAstronomy.Galaxies => 22,
      AstroEquations.AstrophysicsAndAstronomy.Instrumentation => 30,
      AstroEquations.AstrophysicsAndAstronomy.Stars => 27,
      AstroEquations.Physics.Electromagnetism => 85,
      AstroEquations.Physics.Energy => 18,
      AstroEquations.Physics.Forces => 22,
      AstroEquations.Physics.GeneralRelativity => 25,
      AstroEquations.Physics.Materials => 28,
      AstroEquations.Physics.Motion => 54,
      AstroEquations.Physics.NewtonGravity => 18,
      AstroEquations.Physics.Oscillations => 31,
      AstroEquations.Physics.QuantumMechanics => 37,
      AstroEquations.Physics.SpecialRelativity => 26,
      AstroEquations.Physics.Thermodynamics => 42,
      AstroEquations.Physics.Waves => 39,
      AstroEquations.Mathematics.Calculus => 16,
      AstroEquations.Mathematics.Geometry => 25,
      AstroEquations.Mathematics.Notation => 58,
      AstroEquations.Mathematics.Trigonometry => 14,
      AstroEquations.Statistics.StandardDeviation => 18,
      AstroEquations.Statistics.Variance => 14
    }
  end
end
