defmodule AstroEquations.Physics.SpecialRelativity do
  @moduledoc """
  Special Relativity calculations.

  Covers:
  - Lorentz factor γ
  - Time dilation and length contraction
  - Relativistic mass, momentum, energy (total, kinetic, rest)
  - Relativistic velocity addition
  - Four-vectors (position, velocity, momentum)
  - Lorentz boost and Galilean transformation
  - Proper time (numerical integration)
  - Relativistic Doppler effect
  - Invariant spacetime interval
  - Rapidity
  - Relativistic aberration
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Speed in metres per second (m/s). Must satisfy 0 ≤ v < c."
  @type speed :: float()

  @typedoc "Lorentz factor γ ≥ 1."
  @type lorentz_factor :: float()

  @typedoc "Mass in kilograms (kg). Must be positive."
  @type mass :: float()

  @typedoc "Energy in joules (J). Must be non-negative."
  @type energy :: float()

  @typedoc "Momentum in kg·m/s."
  @type momentum :: float()

  # m/s
  @speed_of_light 299_792_458

  # ---------------------------------------------------------------------------
  # Lorentz Factor
  # ---------------------------------------------------------------------------

  @doc """
  Lorentz factor: γ = 1 / √(1 - β²)  where β = v/c

  ## Parameters
    - v: Speed (m/s)
    - c: Speed of light (m/s, default: 299_792_458)

  ## Examples
      iex> SpecialRelativity.gamma_factor(150_000_000) |> Float.round(4)
      1.3416
  """

  @spec gamma_factor(number, number) :: float
  def gamma_factor(v, c \\ @speed_of_light) when is_non_negative(v) and is_positive(c) do
    1 / :math.sqrt(:math.pow(1 - v / c, 2))
  end

  @doc "Normalised velocity: β = v/c."

  @spec beta(number, number) :: float
  def beta(v, c \\ @speed_of_light), do: v / c

  @doc "Rapidity (additive under collinear Lorentz boosts): φ = arctanh(β)."

  @spec rapidity(number, number) :: float
  def rapidity(v, c \\ @speed_of_light) do
    b = v / c
    :math.log((1 + b) / (1 - b)) / 2
  end

  @doc "Speed from rapidity: v = c tanh(φ)."

  @spec speed_from_rapidity(number, number) :: float
  def speed_from_rapidity(phi, c \\ @speed_of_light) do
    c * (:math.exp(2 * phi) - 1) / (:math.exp(2 * phi) + 1)
  end

  # ---------------------------------------------------------------------------
  # Kinematic Effects
  # ---------------------------------------------------------------------------

  @doc """
  Time dilation: t = γ t₀

  ## Examples
      iex> SpecialRelativity.time_dilation(1, 150_000_000) |> Float.round(4)
      1.3416
  """

  @spec time_dilation(number, number, number) :: float
  def time_dilation(t0, v, c \\ @speed_of_light) when is_positive(t0) and is_non_negative(v),
    do: t0 * gamma_factor(v, c)

  @doc """
  Length contraction: L = L₀ / γ

  ## Examples
      iex> SpecialRelativity.length_contraction(1, 150_000_000) |> Float.round(4)
      0.7454
  """

  @spec length_contraction(number, number, number) :: float
  def length_contraction(l0, v, c \\ @speed_of_light), do: l0 / gamma_factor(v, c)

  @doc "Relativistic velocity addition: u' = (u − v)/(1 − uv/c²)."

  @spec relative_velocity(number, number, number) :: float
  def relative_velocity(u, v, c \\ @speed_of_light) do
    (u - v) / (1 - v * u / c * c)
  end

  # ---------------------------------------------------------------------------
  # Relativistic Dynamics
  # ---------------------------------------------------------------------------

  @doc "Relativistic (inertial) mass: m_rel = γm₀."

  @spec relativistic_mass(number, number, number) :: float
  def relativistic_mass(m0, v, c \\ @speed_of_light), do: m0 * gamma_factor(v, c)

  @doc "Relativistic momentum: p = γm₀v."

  @spec relativistic_momentum(number, number, number) :: float
  def relativistic_momentum(m, v, c \\ @speed_of_light)
      when is_positive(m) and is_subluminal(v, c) do
    gamma_factor(v, c) * m * v
  end

  @doc "Rest energy of a body: E₀ = m₀c²."

  @spec rest_energy(number, number) :: number()
  def rest_energy(m, c \\ @speed_of_light) when is_positive(m), do: m * c * c

  @doc "Total relativistic energy: E = γm₀c²."

  @spec total_energy(number, number, number) :: float
  def total_energy(m, v, c \\ @speed_of_light) do
    gamma_factor(v, c) * m * c * c
  end

  @doc "Relativistic kinetic energy: K = (γ − 1)m₀c²."

  @spec kinetic_energy(number, number, number) :: float
  def kinetic_energy(m, v, c \\ @speed_of_light) when is_positive(m) and is_non_negative(v) do
    (gamma_factor(v, c) - 1) * m * c * c
  end

  @doc """
  Energy-momentum relation: E = √((pc)² + (mc²)²)

  ## Examples
      iex> SpecialRelativity.energy_momentum_relation(0, 1.0) > 0
      true
  """

  @spec energy_momentum_relation(number, number, number) :: float
  def energy_momentum_relation(p, m, c \\ @speed_of_light) do
    :math.sqrt(:math.pow(p * c, 2) + :math.pow(m * c * c, 2))
  end

  @doc "Speed from relativistic kinetic energy: v = c √(1 − 1/γ²)."

  @spec speed_from_kinetic_energy(number, number, number) :: float
  def speed_from_kinetic_energy(ke, m, c \\ @speed_of_light) do
    gamma_val = ke / (m * c * c) + 1
    c * :math.sqrt(1 - 1 / gamma_val * gamma_val)
  end

  # ---------------------------------------------------------------------------
  # Relativistic Doppler Effect
  # ---------------------------------------------------------------------------

  @doc """
  Relativistic (longitudinal) Doppler factor: f_obs = f₀ √((1 + β) / (1 - β))

  Positive v = source approaching observer.

  ## Examples
      iex> SpecialRelativity.relativistic_doppler(1.0e14, 0) |> Float.round(2)
      1.0e14
  """

  @spec relativistic_doppler(number, number, number) :: float
  def relativistic_doppler(f0, v, c \\ @speed_of_light) when is_positive(f0) do
    b = v / c
    f0 * :math.sqrt((1 + b) / (1 - b))
  end

  @doc "Transverse Doppler effect (motion perpendicular to line of sight): f_obs = f₀/γ."

  @spec transverse_doppler(number, number, number) :: float
  def transverse_doppler(f0, v, c \\ @speed_of_light) when is_positive(f0) do
    f0 / gamma_factor(v, c)
  end

  @doc """
  Relativistic aberration of light: cos θ' = (cos θ - β) / (1 - β cos θ).

  Converts angle θ in the source frame to angle θ' in the observer frame.

  ## Parameters
    - theta: Angle of incoming light in source frame (radians)
    - v:     Observer speed (m/s)
    - c:     Speed of light (m/s)

  ## Returns
    Aberrated angle in radians

  ## Examples
      iex> SpecialRelativity.relativistic_aberration(:math.pi()/2, 0) |> Float.round(4)
      1.5708
  """

  @spec relativistic_aberration(number, number, number) :: float
  def relativistic_aberration(theta, v, c \\ @speed_of_light) do
    b = v / c
    cos_theta_prime = (:math.cos(theta) - b) / (1 - b * :math.cos(theta))
    :math.acos(max(-1.0, min(1.0, cos_theta_prime)))
  end

  # ---------------------------------------------------------------------------
  # Four-Vectors
  # ---------------------------------------------------------------------------

  @doc "Creates a spacetime four-vector map %{ct:, x:, y:, z:}."

  @spec four_vector(number, number, number, number) :: map()
  def four_vector(ct, x, y, z), do: %{ct: ct, x: x, y: y, z: z}

  @doc "Four-velocity: Uμ = γ(c, vₓ, v_y, v_z)."

  @spec four_velocity(number, number, number, number) :: map()
  def four_velocity(vx, vy, vz, c \\ @speed_of_light) do
    v = :math.sqrt(vx * vx + vy * vy + vz * vz)
    gamma = gamma_factor(v, c)
    four_vector(gamma * c, gamma * vx, gamma * vy, gamma * vz)
  end

  @doc "Four-momentum: Pμ = m₀Uμ."

  @spec four_momentum(number, number, number, number, number) :: map()
  def four_momentum(m, vx, vy, vz, c \\ @speed_of_light) do
    %{ct: e, x: px, y: py, z: pz} = four_velocity(vx, vy, vz, c)
    four_vector(m * e, m * px, m * py, m * pz)
  end

  @doc "Minkowski inner product of two four-vectors (signature −+++)."

  @spec four_product(map(), map()) :: number()
  def four_product(%{ct: ct1, x: x1, y: y1, z: z1}, %{ct: ct2, x: x2, y: y2, z: z2}) do
    -ct1 * ct2 + x1 * x2 + y1 * y2 + z1 * z2
  end

  @doc "Invariant mass from energy and momentum: m = √(E² - (pc)²) / c²."

  @spec invariant_mass(number, number, number) :: float
  def invariant_mass(energy, momentum, c \\ @speed_of_light) do
    :math.sqrt(max(energy * energy - :math.pow(momentum * c, 2), 0.0)) / (c * c)
  end

  # ---------------------------------------------------------------------------
  # Lorentz Transformations
  # ---------------------------------------------------------------------------

  @doc """
  Galilean position transformation (non-relativistic limit): x' = x − vt.

  ## Returns
    {x', y', z'} transformed coordinates
  """

  @spec galilean_transform(number, number, number, number, number) :: {number, number, number}
  def galilean_transform(x, y, z, t, v), do: {x - v * t, y, z}

  @doc """
  Lorentz boost along the x-axis.

  ct' = γ(ct - β x),  x' = γ(x - β ct)

  ## Returns
    {ct', x', y', z'}

  ## Examples
      iex> SpecialRelativity.lorentz_boost(0, 0, 0, 0, 0) |> elem(0) |> Float.round(4)
      0.0
  """

  @spec lorentz_boost(number, number, number, number, number, number) ::
          {number, number, number, number}
  def lorentz_boost(ct, x, y, z, v, c \\ @speed_of_light) do
    gamma = gamma_factor(v, c)
    b = v / c
    {gamma * (ct - b * x), gamma * (x - b * ct), y, z}
  end

  @doc """
  Proper time calculated by numerical integration over a velocity profile.

  τ = ∫ dt / γ(v(t))

  ## Parameters
    - t_a, t_b:    Integration bounds (s)
    - velocity_fn: Function v(t) → speed in m/s
    - c:           Speed of light (m/s)
    - steps:       Integration steps (default: 1000)

  ## Examples
      iex> SpecialRelativity.proper_time(0, 1.0, fn _ -> 0 end) |> Float.round(4)
      1.0
  """

  @spec proper_time(number, number, (number -> number), number, pos_integer) :: float
  def proper_time(t_a, t_b, velocity_fn, c \\ @speed_of_light, steps \\ 1_000) do
    dt = (t_b - t_a) / steps

    Enum.reduce(0..steps, 0.0, fn i, acc ->
      t = t_a + i * dt
      acc + dt / gamma_factor(velocity_fn.(t), c)
    end)
  end

  # ---------------------------------------------------------------------------
  # Invariant Interval
  # ---------------------------------------------------------------------------

  @doc """
  Spacetime interval: s² = -c²Δt² + Δx² + Δy² + Δz²

  Returns:
  - s² < 0 → timelike
  - s² = 0 → lightlike (null)
  - s² > 0 → spacelike
  """

  @spec spacetime_interval(number, number, number, number, number) :: number()

  def spacetime_interval(dt, dx, dy, dz, c \\ @speed_of_light) do
    -:math.pow(c * dt, 2) + dx * dx + dy * dy + dz * dz
  end
end
