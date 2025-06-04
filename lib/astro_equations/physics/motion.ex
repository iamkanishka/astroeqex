defmodule AstroEquations.Physics.Motion do
  @moduledoc """
  This module contains functions related to motion physics concepts including:
  - Velocity
  - Acceleration
  - Newton's Laws
  - Momentum
  - Centripetal Force
  - Kinetic Energy
  - Angular Velocity
  - Angular Acceleration
  - Moment of Inertia (for various shapes)
  - Rotational Kinetic Energy
  - Total Kinetic Energy
  - Angular Momentum
  - Torque
  - Projectile Motion
  - Lagrangian Mechanics
  - Hamiltonian Mechanics
  """

  @doc """
  Calculates velocity given displacement and time.

  ## Parameters
    - displacement: change in position (Δx) in meters
    - time: time interval (Δt) in seconds

  ## Returns
    Velocity in m/s

  ## Examples
      iex> AstroEquations.Physics.Motion.velocity(10, 2)
      5.0
  """
  def velocity(displacement, time) do
    displacement / time
  end

  @doc """
  Calculates acceleration given change in velocity and time.

  ## Parameters
    - delta_v: change in velocity (Δv) in m/s
    - time: time interval (Δt) in seconds

  ## Returns
    Acceleration in m/s²

  ## Examples
      iex> AstroEquations.Physics.Motion.acceleration(20, 5)
      4.0
  """
  def acceleration(delta_v, time) do
    delta_v / time
  end

  @doc """
  Newton's Second Law: calculates force given mass and acceleration.

  ## Parameters
    - mass: in kilograms
    - acceleration: in m/s²

  ## Returns
    Force in Newtons

  ## Examples
      iex> AstroEquations.Physics.Motion.newtons_second_law(5, 2)
      10.0
  """
  def newtons_second_law(mass, acceleration) do
    mass * acceleration
  end

  @doc """
  Calculates momentum given mass and velocity.

  ## Parameters
    - mass: in kilograms
    - velocity: in m/s

  ## Returns
    Momentum in kg·m/s

  ## Examples
      iex> AstroEquations.Physics.Motion.momentum(10, 5)
      50.0
  """
  def momentum(mass, velocity) do
    mass * velocity
  end

  @doc """
  Calculates change in momentum (impulse) given force and time interval.

  ## Parameters
    - force: in Newtons
    - time: time interval in seconds

  ## Returns
    Change in momentum in kg·m/s

  ## Examples
      iex> AstroEquations.Physics.Motion.impulse(10, 5)
      50.0
  """
  def impulse(force, time) do
    force * time
  end

  @doc """
  Calculates centripetal force given mass, velocity, and radius.

  ## Parameters
    - mass: in kilograms
    - velocity: in m/s
    - radius: in meters

  ## Returns
    Centripetal force in Newtons

  ## Examples
      iex> AstroEquations.Physics.Motion.centripetal_force(2, 5, 10)
      5.0
  """
  def centripetal_force(mass, velocity, radius) do
    mass * :math.pow(velocity, 2) / radius
  end

  @doc """
  Calculates kinetic energy given mass and velocity.

  ## Parameters
    - mass: in kilograms
    - velocity: in m/s

  ## Returns
    Kinetic energy in Joules

  ## Examples
      iex> AstroEquations.Physics.Motion.kinetic_energy(4, 5)
      50.0
  """
  def kinetic_energy(mass, velocity) do
    0.5 * mass * :math.pow(velocity, 2)
  end

  @doc """
  Calculates angular velocity given angular displacement and time.

  ## Parameters
    - angular_displacement: change in angle (Δθ) in radians
    - time: time interval (Δt) in seconds

  ## Returns
    Angular velocity in rad/s

  ## Examples
      iex> AstroEquations.Physics.Motion.angular_velocity(:math.pi(), 2)
      1.5707963267948966
  """
  def angular_velocity(angular_displacement, time) do
    angular_displacement / time
  end

  @doc """
  Calculates angular acceleration given change in angular velocity and time.

  ## Parameters
    - delta_omega: change in angular velocity (Δω) in rad/s
    - time: time interval (Δt) in seconds

  ## Returns
    Angular acceleration in rad/s²

  ## Examples
      iex> AstroEquations.Physics.Motion.angular_acceleration(4, 2)
      2.0
  """
  def angular_acceleration(delta_omega, time) do
    delta_omega / time
  end

  @doc """
  Calculates moment of inertia for a point mass.

  ## Parameters
    - mass: in kilograms
    - radius: distance from axis in meters

  ## Returns
    Moment of inertia in kg·m²

  ## Examples
      iex> AstroEquations.Physics.Motion.point_mass_moment_of_inertia(2, 3)
      18.0
  """
  def point_mass_moment_of_inertia(mass, radius) do
    mass * :math.pow(radius, 2)
  end

  @doc """
  Calculates moment of inertia for several point masses.

  ## Parameters
    - masses: list of masses in kilograms
    - radii: list of distances from axis in meters

  ## Returns
    Total moment of inertia in kg·m²

  ## Examples
      iex> AstroEquations.Physics.Motion.multiple_point_masses_moment_of_inertia([1, 2], [3, 4])
      35.0
  """
  def multiple_point_masses_moment_of_inertia(masses, radii) do
    Enum.zip_with(masses, radii, fn m, r -> m * :math.pow(r, 2) end)
    |> Enum.sum()
  end

  @doc """
  Calculates moment of inertia for a thin disk rotating about its center.

  ## Parameters
    - mass: in kilograms
    - radius: in meters

  ## Returns
    Moment of inertia in kg·m²

  ## Examples
      iex> AstroEquations.Physics.Motion.thin_disk_moment_of_inertia(4, 2)
      8.0
  """
  def thin_disk_moment_of_inertia(mass, radius) do
    mass * :math.pow(radius, 2) / 2
  end

  @doc """
  Calculates moment of inertia for a thin loop rotating about its center.

  ## Parameters
    - mass: in kilograms
    - radius: in meters

  ## Returns
    Moment of inertia in kg·m²

  ## Examples
      iex> AstroEquations.Physics.Motion.thin_loop_moment_of_inertia(3, 2)
      12.0
  """
  def thin_loop_moment_of_inertia(mass, radius) do
    mass * :math.pow(radius, 2)
  end

  @doc """
  Calculates moment of inertia for a thin rod rotating about its center.

  ## Parameters
    - mass: in kilograms
    - length: in meters

  ## Returns
    Moment of inertia in kg·m²

  ## Examples
      iex> AstroEquations.Physics.Motion.thin_rod_center_moment_of_inertia(6, 2)
      2.0
  """
  def thin_rod_center_moment_of_inertia(mass, length) do
    mass * :math.pow(length, 2) / 12
  end

  @doc """
  Calculates moment of inertia for a thin rod rotating about its end.

  ## Parameters
    - mass: in kilograms
    - length: in meters

  ## Returns
    Moment of inertia in kg·m²

  ## Examples
      iex> AstroEquations.Physics.Motion.thin_rod_end_moment_of_inertia(6, 2)
      8.0
  """
  def thin_rod_end_moment_of_inertia(mass, length) do
    mass * :math.pow(length, 2) / 3
  end

  @doc """
  Calculates rotational kinetic energy.

  ## Parameters
    - moment_of_inertia: in kg·m²
    - angular_velocity: in rad/s

  ## Returns
    Rotational kinetic energy in Joules

  ## Examples
      iex> AstroEquations.Physics.Motion.rotational_kinetic_energy(2, 3)
      9.0
  """
  def rotational_kinetic_energy(moment_of_inertia, angular_velocity) do
    0.5 * moment_of_inertia * :math.pow(angular_velocity, 2)
  end

  @doc """
  Calculates total kinetic energy (translational + rotational).

  ## Parameters
    - mass: in kilograms
    - velocity: translational velocity in m/s
    - moment_of_inertia: in kg·m²
    - angular_velocity: in rad/s

  ## Returns
    Total kinetic energy in Joules

  ## Examples
      iex> AstroEquations.Physics.Motion.total_kinetic_energy(2, 3, 4, 5)
      59.0
  """
  def total_kinetic_energy(mass, velocity, moment_of_inertia, angular_velocity) do
    translational = 0.5 * mass * :math.pow(velocity, 2)
    rotational = 0.5 * moment_of_inertia * :math.pow(angular_velocity, 2)
    translational + rotational
  end

  @doc """
  Calculates angular momentum.

  ## Parameters
    - moment_of_inertia: in kg·m²
    - angular_velocity: in rad/s

  ## Returns
    Angular momentum in kg·m²/s

  ## Examples
      iex> AstroEquations.Physics.Motion.angular_momentum(3, 4)
      12.0
  """
  def angular_momentum(moment_of_inertia, angular_velocity) do
    moment_of_inertia * angular_velocity
  end

  @doc """
  Calculates torque.

  ## Parameters
    - force: in Newtons
    - lever_arm: distance from axis in meters
    - angle: angle between force and lever arm in radians (defaults to π/2 for perpendicular force)

  ## Returns
    Torque in N·m

  ## Examples
      iex> AstroEquations.Physics.Motion.torque(10, 2)
      20.0
      iex> AstroEquations.Physics.Motion.torque(10, 2, :math.pi()/4)
      14.142135623730951
  """
  def torque(force, lever_arm, angle \\ :math.pi() / 2) do
    force * lever_arm * :math.sin(angle)
  end

  @doc """
  Calculates vertical velocity component using kinematic equation.

  ## Parameters
    - initial_vy: initial vertical velocity (u_y) in m/s
    - acceleration_y: vertical acceleration (a_y) in m/s²
    - delta_y: vertical displacement (Δy) in meters

  ## Returns
    Final vertical velocity in m/s

  ## Examples
      iex> AstroEquations.Physics.Motion.vertical_velocity(10, -9.8, 5)
      8.280123850695517
  """
  def vertical_velocity(initial_vy, acceleration_y, delta_y) do
    :math.sqrt(:math.pow(initial_vy, 2) + 2 * acceleration_y * delta_y)
  end

  @doc """
  Calculates horizontal displacement in projectile motion.

  ## Parameters
    - initial_vx: initial horizontal velocity (u_x) in m/s
    - time: flight time (t) in seconds

  ## Returns
    Horizontal displacement in meters

  ## Examples
      iex> AstroEquations.Physics.Motion.horizontal_displacement(15, 3)
      45.0
  """
  def horizontal_displacement(initial_vx, time) do
    initial_vx * time
  end

  @doc """
  Calculates vertical displacement in projectile motion.

  ## Parameters
    - initial_vy: initial vertical velocity (u_y) in m/s
    - time: flight time (t) in seconds
    - acceleration_y: vertical acceleration (a_y) in m/s²

  ## Returns
    Vertical displacement in meters

  ## Examples
      iex> AstroEquations.Physics.Motion.vertical_displacement(20, 2, -9.8)
      20.4
  """
  def vertical_displacement(initial_vy, time, acceleration_y) do
    initial_vy * time + 0.5 * acceleration_y * :math.pow(time, 2)
  end

  @doc """
  Calculates the Lagrangian (L = T - V).

  ## Parameters
    - kinetic: kinetic energy (T) in Joules
    - potential: potential energy (V) in Joules

  ## Returns
    Lagrangian in Joules

  ## Examples
      iex> AstroEquations.Physics.Motion.lagrangian(50, 20)
      30.0
  """
  def lagrangian(kinetic, potential) do
    kinetic - potential
  end

  @doc """
  Calculates generalized momentum (p_k = ∂L/∂q̇_k).

  Note: This is a placeholder for the partial derivative calculation.
  In practice, you would need to provide the derivative function.

  ## Parameters
    - lagrangian_func: function representing the Lagrangian
    - q_dot: generalized velocity (q̇_k)
    - h: small step for numerical differentiation (default 1.0e-6)

  ## Returns
    Generalized momentum

  ## Examples
      iex> AstroEquations.Physics.Motion.generalized_momentum(fn x -> 0.5 * x ** 2 end, 2)
      2.0000000000575113
  """
  def generalized_momentum(lagrangian_func, q_dot, h \\ 1.0e-6) do
    # Numerical differentiation to approximate ∂L/∂q̇
    (lagrangian_func.(q_dot + h) - lagrangian_func.(q_dot - h)) / (2 * h)
  end

  @doc """
  Calculates the Hamiltonian (H = Σp_l q̇_l - L).

  ## Parameters
    - momenta: list of generalized momenta [p_1, p_2, ...]
    - velocities: list of generalized velocities [q̇_1, q̇_2, ...]
    - lagrangian: Lagrangian value (L)

  ## Returns
    Hamiltonian in Joules

  ## Examples
      iex> AstroEquations.Physics.Motion.hamiltonian([2, 3], [1.5, 2], 10)
      2.0
  """
  def hamiltonian(momenta, velocities, lagrangian) do
    Enum.zip_with(momenta, velocities, fn p, v -> p * v end)
    |> Enum.sum()
    |> Kernel.-(lagrangian)
  end

  @doc """
  Solves Hamilton's equations for simple harmonic motion.

  ## Parameters
    - p: momentum
    - q: position
    - omega: angular frequency

  ## Returns
    Tuple of {dp/dt, dq/dt}

  ## Examples
      iex> AstroEquations.Physics.Motion.hamiltons_equations(2, 3, 1.5)
      {-6.75, 2}
  """
  def hamiltons_equations(p, q, omega) do
    dpdt = -:math.pow(omega, 2) * q
    dqdt = p
    {dpdt, dqdt}
  end
end
