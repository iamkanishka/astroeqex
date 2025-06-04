
defmodule SpecialRelativity do
  @moduledoc """
  A module for special relativity calculations including four-vectors, reference frame transformations,
  and proper time calculations.
  """

  alias __MODULE__

  @doc """
  Represents a four-vector in spacetime.

  ## Fields
    - ct: time component multiplied by speed of light
    - x: x spatial component
    - y: y spatial component
    - z: z spatial component
  """
  defstruct [:ct, :x, :y, :z]
end

defmodule AstroEquations.Physics.SpecialRelativity do
  @moduledoc """
  A module for special relativity calculations including time dilation, length contraction,
  relativistic energy, four-vectors, reference frame transformations, and proper time calculations.


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
    (u - v) / (1 - v * u / c ** 2)
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

  @doc """
  Creates a four-vector from components.

  ## Parameters
    - ct: time component (c*t)
    - x: x spatial component
    - y: y spatial component
    - z: z spatial component

  ## Examples
      iex> SpecialRelativity.four_vector(1, 2, 3, 4)
      %SpecialRelativity{ct: 1, x: 2, y: 3, z: 4}
  """
  @spec four_vector(number, number, number, number) :: %SpecialRelativity{}
  def four_vector(ct, x, y, z) do
    %SpecialRelativity{ct: ct, x: x, y: y, z: z}
  end

  @doc """
  Calculates four-velocity given velocity components and proper time.

  ## Parameters
    - vx: x-component of velocity (m/s)
    - vy: y-component of velocity (m/s)
    - vz: z-component of velocity (m/s)
    - c: speed of light (default: 299_792_458 m/s)

  ## Examples
      iex> SpecialRelativity.four_velocity(0.5 * 299792458, 0, 0)
      %SpecialRelativity{ct: 1.1547005383792517, x: 0.5773502691896258, y: 0, z: 0}
  """
  @spec four_velocity(number, number, number, number) :: %SpecialRelativity{}
  def four_velocity(vx, vy, vz, c \\ 299_792_458) do
    γ = gamma_factor(:math.sqrt(vx ** 2 + vy ** 2 + vz ** 2), c)
    four_vector(γ * c, γ * vx, γ * vy, γ * vz)
  end

  @doc """
  Calculates four-momentum given mass and velocity components.

  ## Parameters
    - m: mass (kg)
    - vx: x-component of velocity (m/s)
    - vy: y-component of velocity (m/s)
    - vz: z-component of velocity (m/s)
    - c: speed of light (default: 299_792_458 m/s)

  ## Examples
      iex> SpecialRelativity.four_momentum(1, 0.5 * 299792458, 0, 0)
      %SpecialRelativity{ct: 1.1547005383792517e8, x: 0.5773502691896258e8, y: 0, z: 0}
  """
  @spec four_momentum(number, number, number, number, number) :: %SpecialRelativity{}
  def four_momentum(m, vx, vy, vz, c \\ 299_792_458) do
    %SpecialRelativity{ct: e, x: px, y: py, z: pz} = four_velocity(vx, vy, vz, c)
    four_vector(m * e, m * px, m * py, m * pz)
  end

  @doc """
  Applies a Galilean transformation to spatial coordinates.

  ## Parameters
    - x: x coordinate
    - y: y coordinate
    - z: z coordinate
    - t: time
    - v: relative velocity between frames (m/s)

  ## Examples
      iex> SpecialRelativity.galilean_transform(10, 5, 3, 2, 5)
      {20, 5, 3}
  """
  @spec galilean_transform(number, number, number, number, number) :: {number, number, number}
  def galilean_transform(x, y, z, t, v) do
    {x + v * t, y, z}
  end

  @doc """
  Applies a Lorentz boost to spacetime coordinates.

  ## Parameters
    - ct: time component (c*t)
    - x: x coordinate
    - y: y coordinate
    - z: z coordinate
    - v: relative velocity between frames (m/s)
    - c: speed of light (default: 299_792_458 m/s)

  ## Examples
      iex> SpecialRelativity.lorentz_boost(1, 0, 0, 0, 0.5 * 299792458)
      {1.1547005383792517, -0.5773502691896258, 0, 0}
  """
  @spec lorentz_boost(number, number, number, number, number, number) ::
          {number, number, number, number}
  def lorentz_boost(ct, x, y, z, v, c \\ 299_792_458) do
    γ = gamma_factor(v, c)
    β = v / c

    ct_prime = γ * (ct - β * x)
    x_prime = γ * (x - β * ct)
    {ct_prime, x_prime, y, z}
  end

  @doc """
  Calculates proper time between two events for a given velocity profile.

  ## Parameters
    - t_a: start time
    - t_b: end time
    - velocity_fn: function that returns velocity at time t
    - c: speed of light (default: 299_792_458 m/s)

  ## Examples
      iex> SpecialRelativity.proper_time(0, 1, fn _ -> 0.5 * 299792458 end)
      0.8660254037844386
  """
  @spec proper_time(number, number, (number -> number), number) :: number
  def proper_time(t_a, t_b, velocity_fn, c \\ 299_792_458) do
    # Numerical integration using trapezoidal rule
    steps = 1000
    delta_t = (t_b - t_a) / steps

    Enum.reduce(0..steps, 0, fn i, acc ->
      t = t_a + i * delta_t
      v = velocity_fn.(t)
      γ = gamma_factor(v, c)
      acc + delta_t / γ
    end)
  end

end
