defmodule AstroEquations.Physics.Energy do
  @moduledoc """
  Fundamental energy calculations in classical and modern physics.

  Covers:
  - Mechanical work (dot product and angle form)
  - Kinetic and rotational kinetic energy
  - Potential energy (gravitational, spring, electric)
  - Power (average, instantaneous, vector form)
  - Work-energy theorem
  - Escape velocity
  - Relativistic rest energy and binding energy
  - Simple harmonic motion energy
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Mass in kilograms (kg). Must be positive."
  @type mass :: float()

  @typedoc "Velocity in metres per second (m/s). Non-negative."
  @type velocity :: float()

  @typedoc "Energy in joules (J)."
  @type energy :: float()

  @typedoc "Angular velocity in radians per second (rad/s)."
  @type angular_velocity :: float()

  # m³ kg⁻¹ s⁻²
  @gravitational_constant 6.674_30e-11
  # m/s
  @speed_of_light 2.997_924_58e8

  # ---------------------------------------------------------------------------
  # Work
  # ---------------------------------------------------------------------------

  @doc """
  Work done by a constant force along a displacement: W = F · d

  ## Parameters
    - force_vector:        [Fx, Fy, Fz]
    - displacement_vector: [dx, dy, dz]

  ## Returns
    Work in joules

  ## Examples
      iex> Energy.work([2, 0, 0], [3, 0, 0])
      6.0
      iex> Energy.work([1, 2, 3], [4, 5, 6])
      32.0
  """

  @spec work([number], [number]) :: float
  def work(force_vector, displacement_vector) do
    Enum.zip(force_vector, displacement_vector)
    |> Enum.map(fn {f, d} -> f * d end)
    |> Enum.sum()
    |> Kernel.*(1.0)
  end

  @doc """
  Work done by a force at angle θ: W = F d cos θ

  ## Examples
      iex> Energy.work_angle(10, 5, 0)
      50.0
  """

  @spec work_angle(number, number, number) :: float
  def work_angle(force, displacement, theta \\ 0.0) do
    force * displacement * :math.cos(theta)
  end

  # ---------------------------------------------------------------------------
  # Kinetic Energy
  # ---------------------------------------------------------------------------

  @doc """
  Translational kinetic energy: KE = ½ m v²

  ## Examples
      iex> Energy.kinetic_energy(4, 5)
      50.0
  """

  @spec kinetic_energy(number, number) :: float
  def kinetic_energy(mass, velocity) when is_positive(mass) and is_non_negative(velocity) do
    0.5 * mass * velocity * velocity
  end

  @doc """
  Kinetic energy from momentum: KE = p² / (2m)

  ## Examples
      iex> Energy.kinetic_energy_from_momentum(10, 2)
      25.0
  """

  @spec kinetic_energy_from_momentum(number, number) :: float
  def kinetic_energy_from_momentum(momentum, mass) when is_positive(mass) do
    :math.pow(momentum, 2) / (2 * mass)
  end

  @doc """
  Rotational kinetic energy: KE_rot = ½ I omega²

  ## Examples
      iex> Energy.rotational_kinetic_energy(2, 3)
      9.0
  """

  @spec rotational_kinetic_energy(number, number) :: float
  def rotational_kinetic_energy(moment_of_inertia, angular_velocity) do
    0.5 * moment_of_inertia * :math.pow(angular_velocity, 2)
  end

  # ---------------------------------------------------------------------------
  # Potential Energy
  # ---------------------------------------------------------------------------

  @doc """
  Gravitational potential energy near Earth's surface: PE = m g h

  ## Examples
      iex> Energy.gravitational_pe(10, 5)
      490.5
  """

  @spec gravitational_pe(number, number, number) :: number()
  def gravitational_pe(mass, height, g \\ 9.81) when is_positive(mass), do: mass * g * height

  @doc """
  General gravitational potential energy: PE = -G m M / r

  ## Examples
      iex> Energy.gravitational_pe_general(5.972e24, 7.348e22, 3.844e8) < 0
      true
  """

  @spec gravitational_pe_general(number, number, number) :: float
  def gravitational_pe_general(m1, m2, r) do
    -@gravitational_constant * m1 * m2 / r
  end

  @doc """
  Elastic (spring) potential energy: PE = ½ k x²

  ## Examples
      iex> Energy.spring_pe(100, 0.1)
      0.5
  """

  @spec spring_pe(number, number) :: float
  def spring_pe(k, x) when is_positive(k), do: 0.5 * k * :math.pow(x, 2)

  @doc """
  Electric potential energy of two charges: U = q₁ q₂ / (4πε₀ r)

  ## Examples
      iex> Energy.electric_pe(1.0e-9, 1.0e-9, 0.1) > 0
      true
  """

  @spec electric_pe(number, number, number, number) :: float
  def electric_pe(q1, q2, r, epsilon_0 \\ 8.854_187_812_8e-12) do
    q1 * q2 / (4 * :math.pi() * epsilon_0 * r)
  end

  # ---------------------------------------------------------------------------
  # Power
  # ---------------------------------------------------------------------------

  @doc """
  Average power: P = W / t

  ## Examples
      iex> Energy.power_from_work(100, 5)
      20.0
  """

  @spec power_from_work(number, number) :: float
  def power_from_work(work, time), do: work / time

  @doc """
  Instantaneous mechanical power: P = F v cos θ

  ## Examples
      iex> Energy.mechanical_power(50, 10)
      500.0
  """

  @spec mechanical_power(number, number, number) :: float
  def mechanical_power(force, velocity, theta \\ 0.0) do
    force * velocity * :math.cos(theta)
  end

  @doc """
  Instantaneous mechanical power as the dot product of force and velocity: P = F·v.

  ## Parameters
    - force_vector:    [Fx, Fy, Fz]
    - velocity_vector: [vx, vy, vz]
  """

  @spec power_vectors([number], [number]) :: float
  def power_vectors(force_vector, velocity_vector) do
    Enum.zip(force_vector, velocity_vector)
    |> Enum.map(fn {f, v} -> f * v end)
    |> Enum.sum()
    |> Kernel.*(1.0)
  end

  # ---------------------------------------------------------------------------
  # Work-Energy Theorem
  # ---------------------------------------------------------------------------

  @doc """
  Change in kinetic energy equals net work: ΔKE = W_net.

  Returns the final kinetic energy after net work is applied.

  ## Examples
      iex> Energy.work_energy_theorem(100, 50)
      150.0
  """

  @spec work_energy_theorem(number, number) :: number()
  def work_energy_theorem(initial_ke, net_work), do: initial_ke + net_work

  @doc """
  Escape velocity from a gravitational body: v_esc = √(2 G M / r)

  ## Examples
      iex> Energy.escape_velocity(5.972e24, 6.371e6) > 0
      true
  """

  @spec escape_velocity(number, number) :: float
  def escape_velocity(mass, radius) when is_positive(mass) and is_positive(radius) do
    :math.sqrt(2 * @gravitational_constant * mass / radius)
  end

  # ---------------------------------------------------------------------------
  # Relativistic Energy
  # ---------------------------------------------------------------------------

  @doc """
  Rest energy: E₀ = m c²

  ## Examples
      iex> Energy.rest_energy(1.0) > 0
      true
  """

  @spec rest_energy(number, number) :: number()
  def rest_energy(mass, c \\ @speed_of_light) when is_positive(mass), do: mass * c * c

  @doc """
  Mass-energy equivalent from binding (mass defect): ΔE = Δm c²

  ## Examples
      iex> Energy.binding_energy(3.565e-29) > 0
      true
  """

  @spec binding_energy(number, number) :: number()
  def binding_energy(delta_mass, c \\ @speed_of_light), do: delta_mass * c * c

  # ---------------------------------------------------------------------------
  # Simple Harmonic Motion Energy
  # ---------------------------------------------------------------------------

  @doc """
  Total mechanical energy in SHM (constant): E = ½ k A²

  ## Examples
      iex> Energy.shm_total_energy(10, 0.1) |> Float.round(4)
      0.05
  """

  @spec shm_total_energy(number, number) :: float
  def shm_total_energy(k, amplitude), do: 0.5 * k * :math.pow(amplitude, 2)

  @doc """
  Instantaneous KE in SHM: KE = ½ m omega² (A² - x²)

  ## Examples
      iex> Energy.shm_kinetic_energy(1.0, 10, 0.1, 0.0) |> Float.round(4)
      0.5
  """

  @spec shm_kinetic_energy(number, number, number, number) :: float

  def shm_kinetic_energy(mass, omega, amplitude, x) when is_positive(mass) do
    0.5 * mass * :math.pow(omega, 2) * (:math.pow(amplitude, 2) - :math.pow(x, 2))
  end
end
