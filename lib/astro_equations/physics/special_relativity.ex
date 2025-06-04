defmodule AstroEquations.Physics.SpecialRelativity do
  @moduledoc """
  A module for special relativity calculations including time dilation, length contraction,
  relativistic energy, and more.

  All formulas are based on the principles of special relativity with the speed of light c.
  """

  @doc """
  Calculates the gamma factor (Lorentz factor) for a given velocity.

  ## Parameters
    - v: velocity in meters/second
    - c: speed of light in meters/second (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.gamma_factor(150_000_000)
      1.3416407864998738
  """
  @spec gamma_factor(number, number) :: float
  def gamma_factor(v, c \\ 299_792_458) do
    1 / :math.sqrt(1 - (v / c) ** 2)
  end

  @doc """
  Calculates relativistic time dilation.

  ## Parameters
    - t0: proper time (time in the rest frame) in seconds
    - v: relative velocity in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.time_dilation(1, 150_000_000)
      1.3416407864998738
  """
  @spec time_dilation(number, number, number) :: float
  def time_dilation(t0, v, c \\ 299_792_458) do
    t0 * gamma_factor(v, c)
  end

  @doc """
  Calculates length contraction.

  ## Parameters
    - l0: proper length (length in the rest frame) in meters
    - v: relative velocity in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.length_contraction(1, 150_000_000)
      0.7453559924999299
  """
  @spec length_contraction(number, number, number) :: float
  def length_contraction(l0, v, c \\ 299_792_458) do
    l0 / gamma_factor(v, c)
  end

  @doc """
  Calculates relativistic mass.

  ## Parameters
    - m0: rest mass in kg
    - v: velocity in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.relativistic_mass(1, 150_000_000)
      1.3416407864998738
  """
  @spec relativistic_mass(number, number, number) :: float
  def relativistic_mass(m0, v, c \\ 299_792_458) do
    m0 * gamma_factor(v, c)
  end

  @doc """
  Calculates rest energy (E = mc²).

  ## Parameters
    - m: mass in kg
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.rest_energy(1)
      8.987551787368176e16
  """
  @spec rest_energy(number, number) :: float
  def rest_energy(m, c \\ 299_792_458) do
    m * c ** 2
  end

  @doc """
  Calculates total relativistic energy (E = γmc²).

  ## Parameters
    - m: rest mass in kg
    - v: velocity in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.total_energy(1, 150_000_000)
      1.2057362987360216e17
  """
  @spec total_energy(number, number, number) :: float
  def total_energy(m, v, c \\ 299_792_458) do
    gamma_factor(v, c) * rest_energy(m, c)
  end

  @doc """
  Calculates relativistic kinetic energy (K = (γ-1)mc²).

  ## Parameters
    - m: rest mass in kg
    - v: velocity in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.kinetic_energy(1, 150_000_000)
      3.0698111999220405e16
  """
  @spec kinetic_energy(number, number, number) :: float
  def kinetic_energy(m, v, c \\ 299_792_458) do
    (gamma_factor(v, c) - 1) * rest_energy(m, c)
  end

  @doc """
  Calculates relative velocity in one dimension.

  ## Parameters
    - u: velocity of object in frame S in meters/second
    - v: velocity of frame S' relative to S in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.relative_velocity(200_000_000, 150_000_000)
      111111111.1111111
  """
  @spec relative_velocity(number, number, number) :: float
  def relative_velocity(u, v, c \\ 299_792_458) do
    (u - v) / (1 - (v * u) / c ** 2)
  end

  @doc """
  Calculates relativistic momentum.

  ## Parameters
    - m: rest mass in kg
    - v: velocity in meters/second
    - c: speed of light (default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.relativistic_momentum(1, 150_000_000)
      2.0124611797498107e8
  """
  @spec relativistic_momentum(number, number, number) :: float
  def relativistic_momentum(m, v, c \\ 299_792_458) do
    gamma_factor(v, c) * m * v
  end
end
