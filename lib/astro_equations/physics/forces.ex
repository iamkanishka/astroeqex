defmodule AstroEquations.Physics.Forces do
  @moduledoc """
  Fundamental force calculations in classical mechanics and fluid dynamics.

  Covers:
  - Newton's laws
  - Buoyancy (Archimedes)
  - Friction (kinetic, static, rolling, viscous drag)
  - Spring force (Hooke's law)
  - Centripetal force and acceleration
  - Gravitational force
  - Pressure and normal force
  - Drag (linear and quadratic)
  - Terminal velocity
  - Tension and impulse / momentum
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Force in newtons (N)."
  @type force :: float()

  @typedoc "Mass in kilograms (kg). Must be positive."
  @type mass :: float()

  @typedoc "Coefficient of friction or drag (dimensionless). Non-negative."
  @type coefficient :: float()

  @gravitational_constant 6.674_30e-11

  # ---------------------------------------------------------------------------
  # Newton's Laws
  # ---------------------------------------------------------------------------

  @doc """
  Newton's Second Law: F = m a

  ## Examples
      iex> Forces.newtons_second_law(5, 2)
      10.0
  """

  @spec newtons_second_law(number, number) :: number()
  def newtons_second_law(mass, acceleration) when is_positive(mass), do: mass * acceleration

  @doc """
  Newton's Law of Universal Gravitation: F = G m₁ m₂ / r²

  ## Parameters
    - m1, m2: Masses (kg)
    - r:      Separation distance (m)
  """

  @spec gravitational_force(number, number, number) :: float
  def gravitational_force(m1, m2, r) when r > 0 do
    @gravitational_constant * m1 * m2 / :math.pow(r, 2)
  end

  @doc """
  Weight of an object: W = m g

  ## Examples
      iex> Forces.weight(70) |> Float.round(2)
      686.7
  """

  @spec weight(number, number) :: number()
  def weight(mass, g \\ 9.81) when is_positive(mass), do: mass * g

  # ---------------------------------------------------------------------------
  # Buoyancy
  # ---------------------------------------------------------------------------

  @doc """
  Buoyant force: F_b = m_fluid g  (displaced mass form)

  ## Examples
      iex> Forces.buoyancy(5)
      49.05
  """

  @spec buoyancy(number, number) :: number()
  def buoyancy(displaced_mass, gravity \\ 9.81), do: displaced_mass * gravity

  @doc """
  Buoyant force from fluid density and displaced volume: F_b = ρ V g

  ## Examples
      iex> Forces.buoyancy_from_density(1000, 0.005)
      49.05
  """

  @spec buoyancy_from_density(number, number, number) :: number()
  def buoyancy_from_density(density, volume, gravity \\ 9.81) when is_positive(density),
    do: density * volume * gravity

  # ---------------------------------------------------------------------------
  # Friction
  # ---------------------------------------------------------------------------

  @doc """
  Kinetic friction force: f_k = μ_k N

  ## Examples
      iex> Forces.kinetic_friction(0.3, 10)
      3.0
  """

  @spec kinetic_friction(number, number) :: number()
  def kinetic_friction(coefficient, normal_force), do: coefficient * normal_force

  @doc """
  Maximum static friction force: f_s ≤ μ_s N

  ## Examples
      iex> Forces.static_friction(0.4, 10)
      4.0
  """

  @spec static_friction(number, number) :: number()
  def static_friction(coefficient, normal_force), do: coefficient * normal_force

  @doc """
  Rolling friction force: f_r = μ_r N

  ## Examples
      iex> Forces.rolling_friction(0.01, 1000)
      10.0
  """

  @spec rolling_friction(number, number) :: number()
  def rolling_friction(coefficient, normal_force), do: coefficient * normal_force

  # ---------------------------------------------------------------------------
  # Spring Force
  # ---------------------------------------------------------------------------

  @doc """
  Hooke's Law: F = -k x

  ## Examples
      iex> Forces.spring_force(10, 0.5)
      -5.0
  """

  @spec spring_force(number, number) :: number()
  def spring_force(k, x), do: -k * x

  # ---------------------------------------------------------------------------
  # Centripetal / Centrifugal
  # ---------------------------------------------------------------------------

  @doc """
  Centripetal force: F_c = m v² / r

  ## Examples
      iex> Forces.centripetal_force(2, 5, 10)
      5.0
  """

  @spec centripetal_force(number, number, number) :: float
  def centripetal_force(mass, velocity, radius) when is_positive(mass) and is_positive(radius) do
    mass * :math.pow(velocity, 2) / radius
  end

  @doc """
  Centripetal acceleration: a_c = v² / r

  ## Examples
      iex> Forces.centripetal_acceleration(5, 10)
      2.5
  """

  @spec centripetal_acceleration(number, number) :: float
  def centripetal_acceleration(velocity, radius) when is_positive(radius),
    do: :math.pow(velocity, 2) / radius

  @doc """
  Centripetal force using angular velocity: F_c = m omega² r

  ## Examples
      iex> Forces.centripetal_force_angular(1000, 2, 50)
      200000.0
  """

  @spec centripetal_force_angular(number, number, number) :: float
  def centripetal_force_angular(mass, angular_velocity, radius) when is_positive(mass) do
    mass * :math.pow(angular_velocity, 2) * radius
  end

  # ---------------------------------------------------------------------------
  # Drag Forces
  # ---------------------------------------------------------------------------

  @doc """
  Linear (Stokes) drag force: F_d = b v.

  Applicable at low Reynolds numbers.

  ## Examples
      iex> Forces.stokes_drag(0.5, 3.0)
      1.5
  """

  @spec stokes_drag(number, number) :: number()
  def stokes_drag(drag_coefficient, velocity), do: drag_coefficient * velocity

  @doc """
  Stokes drag on a sphere: F_d = 6π η r v

  ## Examples
      iex> Forces.stokes_drag_sphere(1.0e-3, 0.001, 1.0) > 0
      true
  """

  @spec stokes_drag_sphere(number, number, number) :: float
  def stokes_drag_sphere(eta, radius, velocity) when is_positive(radius) do
    6 * :math.pi() * eta * radius * velocity
  end

  @doc """
  Quadratic drag force: F_d = ½ C_d ρ A v².

  Applicable at high Reynolds numbers.

  ## Examples
      iex> Forces.quadratic_drag(0.47, 1.225, 0.04, 30) > 0
      true
  """

  @spec quadratic_drag(number, number, number, number) :: float
  def quadratic_drag(drag_coeff, fluid_density, area, velocity) do
    0.5 * drag_coeff * fluid_density * area * :math.pow(velocity, 2)
  end

  @doc """
  Terminal velocity (quadratic drag): v_t = √(2mg / (C_d ρ A))

  ## Examples
      iex> Forces.terminal_velocity(75, 0.47, 1.225, 0.6) > 0
      true
  """

  @spec terminal_velocity(number, number, number, number, number) :: float
  def terminal_velocity(mass, drag_coeff, fluid_density, area, g \\ 9.81)
      when is_positive(mass)
      when is_positive(mass) and is_positive(drag_coeff) and is_positive(fluid_density) and
             is_positive(area) do
    :math.sqrt(2 * mass * g / (drag_coeff * fluid_density * area))
  end

  # ---------------------------------------------------------------------------
  # Pressure & Normal Force
  # ---------------------------------------------------------------------------

  @doc """
  Hydrostatic pressure: P = ρ g h

  ## Examples
      iex> Forces.hydrostatic_pressure(1000, 10) |> Float.round(1)
      98100.0
  """

  @spec hydrostatic_pressure(number, number, number) :: number()
  def hydrostatic_pressure(density, height, g \\ 9.81) when is_positive(density),
    do: density * g * height

  @doc """
  Pressure from force and area: P = F / A

  ## Examples
      iex> Forces.pressure(100, 0.01)
      10000.0
  """

  @spec pressure(number, number) :: float
  def pressure(force, area), do: force / area

  @doc """
  Normal force on an inclined plane: N = m g cos θ

  ## Examples
      iex> Forces.normal_force_incline(10, 0) |> Float.round(2)
      98.1
  """

  @spec normal_force_incline(number, number, number) :: float
  def normal_force_incline(mass, theta, g \\ 9.81) when is_positive(mass) do
    mass * g * :math.cos(theta)
  end

  @doc """
  Force component along an inclined plane: F_// = m g sin θ

  ## Examples
      iex> Forces.incline_force_parallel(10, :math.pi()/6) |> Float.round(2)
      49.05
  """

  @spec incline_force_parallel(number, number, number) :: float
  def incline_force_parallel(mass, theta, g \\ 9.81) when is_positive(mass) do
    mass * g * :math.sin(theta)
  end

  # ---------------------------------------------------------------------------
  # Impulse & Momentum
  # ---------------------------------------------------------------------------

  @doc """
  Impulse: J = F Δt

  ## Examples
      iex> Forces.impulse(10, 5)
      50.0
  """

  @spec impulse(number, number) :: number()
  def impulse(force, time), do: force * time

  @doc """
  Linear momentum: p = m v

  ## Examples
      iex> Forces.momentum(10, 5)
      50.0
  """

  @spec momentum(number, number) :: number()

  def momentum(mass, velocity) when is_positive(mass), do: mass * velocity
end
