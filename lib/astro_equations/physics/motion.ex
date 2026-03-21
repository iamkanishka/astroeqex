defmodule AstroEquations.Physics.Motion do
  @moduledoc """
  Classical mechanics of motion: kinematics, dynamics, rotation, and analytical mechanics.

  Covers:
  - Kinematics (velocity, acceleration, SUVAT equations)
  - Newton's Second Law, momentum, impulse
  - Centripetal dynamics
  - Kinetic energy (translational and rotational)
  - Moments of inertia (point mass, disk, ring, rod, sphere, hollow sphere)
  - Parallel-axis theorem
  - Angular momentum and torque
  - Projectile motion
  - Circular and simple harmonic motion
  - Lagrangian and Hamiltonian mechanics
  - Orbital mechanics (Kepler's laws)
  """

  @gravitational_constant 6.67430e-11

  # ---------------------------------------------------------------------------
  # Kinematics
  # ---------------------------------------------------------------------------

  @doc "Average velocity: v = Δx/Δt."
  @spec velocity(number, number) :: float
  def velocity(displacement, time), do: displacement / time * 1.0

  @doc "Average acceleration: a = Δv/Δt."
  @spec acceleration(number, number) :: float
  def acceleration(delta_v, time), do: delta_v / time * 1.0

  @doc "SUVAT final velocity: v = u + at."
  @spec final_velocity(number, number, number) :: float
  def final_velocity(u, a, t), do: u + a * t * 1.0

  @doc "SUVAT displacement: s = ut + ½at²."
  @spec displacement_suvat(number, number, number) :: float
  def displacement_suvat(u, t, a), do: u * t + 0.5 * a * t * t * 1.0

  @doc "SUVAT speed from v² = u² + 2as (returns non-negative magnitude)."
  @spec velocity_squared(number, number, number) :: float
  def velocity_squared(u, a, s), do: :math.sqrt(max(:math.pow(u, 2) + 2 * a * s, 0.0))

  @doc "SUVAT displacement: s = (u + v)t/2."
  @spec displacement_average(number, number, number) :: float
  def displacement_average(u, v, t), do: (u + v) * t / 2 * 1.0

  # ---------------------------------------------------------------------------
  # Newton's Laws & Momentum
  # ---------------------------------------------------------------------------

  @doc "Newton's Second Law of Motion: F = ma."
  @spec newtons_second_law(number, number) :: float
  def newtons_second_law(mass, accel), do: mass * accel * 1.0

  @doc "Linear momentum: p = mv."
  @spec momentum(number, number) :: float
  def momentum(mass, velocity), do: mass * velocity * 1.0

  @doc "Impulse: J = FΔt = Δp."
  @spec impulse(number, number) :: float
  def impulse(force, time), do: force * time * 1.0

  @doc "Acceleration from net force: a = F_net/m."
  @spec acceleration_from_force(number, number) :: float
  def acceleration_from_force(net_force, mass), do: net_force / mass * 1.0

  # ---------------------------------------------------------------------------
  # Kinetic Energy
  # ---------------------------------------------------------------------------

  @doc "Translational kinetic energy: KE = ½mv²."
  @spec kinetic_energy(number, number) :: float
  def kinetic_energy(mass, velocity), do: 0.5 * mass * :math.pow(velocity, 2)

  @doc "Rotational kinetic energy: KE_rot = ½Iomega²."
  @spec rotational_kinetic_energy(number, number) :: float
  def rotational_kinetic_energy(moment_of_inertia, angular_velocity) do
    0.5 * moment_of_inertia * :math.pow(angular_velocity, 2)
  end

  @doc "Total kinetic energy (translational + rotational): KE_total = ½mv² + ½Iomega²."
  @spec total_kinetic_energy(number, number, number, number) :: float
  def total_kinetic_energy(mass, velocity, moment_of_inertia, angular_velocity) do
    0.5 * mass * :math.pow(velocity, 2) +
      0.5 * moment_of_inertia * :math.pow(angular_velocity, 2)
  end

  # ---------------------------------------------------------------------------
  # Centripetal Motion
  # ---------------------------------------------------------------------------

  @doc "Centripetal force required for circular motion: F_c = mv²/r."
  @spec centripetal_force(number, number, number) :: float
  def centripetal_force(mass, velocity, radius) do
    mass * :math.pow(velocity, 2) / radius
  end

  @doc "Centripetal acceleration: a_c = v²/r."
  @spec centripetal_acceleration(number, number) :: float
  def centripetal_acceleration(velocity, radius), do: :math.pow(velocity, 2) / radius

  @doc "Circular orbital speed from gravitational balance: v = √(GM/r)."
  @spec orbital_speed(number, number) :: float
  def orbital_speed(central_mass, radius) do
    :math.sqrt(@gravitational_constant * central_mass / radius)
  end

  # ---------------------------------------------------------------------------
  # Rotation
  # ---------------------------------------------------------------------------

  @doc "Average angular velocity: omega = Δθ/Δt."
  @spec angular_velocity(number, number) :: float
  def angular_velocity(angular_displacement, time), do: angular_displacement / time * 1.0

  @doc "Average angular acceleration: α = Δomega/Δt."
  @spec angular_acceleration(number, number) :: float
  def angular_acceleration(delta_omega, time), do: delta_omega / time * 1.0

  @doc "Tangential (linear) velocity of a rotating body: v_t = omegar."
  @spec tangential_velocity(number, number) :: float
  def tangential_velocity(omega, radius), do: omega * radius * 1.0

  @doc "Tangential acceleration: a_t = αr."
  @spec tangential_acceleration(number, number) :: float
  def tangential_acceleration(alpha, radius), do: alpha * radius * 1.0

  @doc "Angular momentum of a rigid body: L = Iomega."
  @spec angular_momentum(number, number) :: float
  def angular_momentum(moment_of_inertia, angular_velocity) do
    moment_of_inertia * angular_velocity * 1.0
  end

  @doc "Torque about a pivot: τ = r × F (magnitude: rF sin θ)."
  @spec torque(number, number, number) :: float
  def torque(force, lever_arm, angle \\ :math.pi() / 2) do
    force * lever_arm * :math.sin(angle)
  end

  @doc "Newton's Second Law for rotation: τ = Iα."
  @spec rotational_newtons_law(number, number) :: float
  def rotational_newtons_law(moment_of_inertia, angular_accel) do
    moment_of_inertia * angular_accel * 1.0
  end

  @doc "Angular impulse: L_impulse = τΔt."
  @spec angular_impulse(number, number) :: float
  def angular_impulse(torque, time), do: torque * time * 1.0

  @doc """
  Coriolis acceleration in a rotating frame: a_Cor = -2 Ω × v.

  Returns the magnitude 2 Ω v sin θ where θ is the angle between Ω and v.

  ## Parameters
    - omega_rot:    Rotation rate Ω of the frame (rad/s)
    - velocity:     Object speed v in the rotating frame (m/s)
    - theta:        Angle between Ω and v (default π/2 — fully perpendicular)

  ## Examples
      iex> Motion.coriolis_acceleration(7.3e-5, 300, :math.pi()/2) |> Float.round(4)
      0.0438
  """
  @spec coriolis_acceleration(number, number, number) :: float
  def coriolis_acceleration(omega_rot, velocity, theta \\ :math.pi() / 2) do
    2 * omega_rot * velocity * :math.sin(theta)
  end

  # ---------------------------------------------------------------------------
  # Moments of Inertia
  # ---------------------------------------------------------------------------

  @doc "Moment of inertia of a point mass: I = mr²."
  @spec point_mass_moment_of_inertia(number, number) :: float
  def point_mass_moment_of_inertia(mass, radius) do
    mass * :math.pow(radius, 2)
  end

  @doc "Moment of inertia of a system of point masses: I = Σ mᵢrᵢ²."
  @spec multiple_point_masses_moment_of_inertia([number], [number]) :: float
  def multiple_point_masses_moment_of_inertia(masses, radii) do
    Enum.zip_with(masses, radii, fn m, r -> m * :math.pow(r, 2) end) |> Enum.sum()
  end

  @doc "Moment of inertia of a solid disk about its central axis: I = ½mR²."
  @spec thin_disk_moment_of_inertia(number, number) :: float
  def thin_disk_moment_of_inertia(mass, radius) do
    mass * :math.pow(radius, 2) / 2
  end

  @doc "Moment of inertia of a thin ring about its central axis: I = mR²."
  @spec thin_loop_moment_of_inertia(number, number) :: float
  def thin_loop_moment_of_inertia(mass, radius) do
    mass * :math.pow(radius, 2)
  end

  @doc "Moment of inertia of a thin rod about its midpoint: I = mL²/12."
  @spec thin_rod_center_moment_of_inertia(number, number) :: float
  def thin_rod_center_moment_of_inertia(mass, length) do
    mass * :math.pow(length, 2) / 12
  end

  @doc "Moment of inertia of a thin rod about one end: I = mL²/3."
  @spec thin_rod_end_moment_of_inertia(number, number) :: float
  def thin_rod_end_moment_of_inertia(mass, length) do
    mass * :math.pow(length, 2) / 3
  end

  @doc "Moment of inertia of a solid sphere about a diameter: I = 2mR²/5."
  @spec solid_sphere_moment_of_inertia(number, number) :: float
  def solid_sphere_moment_of_inertia(mass, radius) do
    2 / 5 * mass * :math.pow(radius, 2)
  end

  @doc "Moment of inertia of a hollow spherical shell about a diameter: I = 2mR²/3."
  @spec hollow_sphere_moment_of_inertia(number, number) :: float
  def hollow_sphere_moment_of_inertia(mass, radius) do
    2 / 3 * mass * :math.pow(radius, 2)
  end

  @doc "Moment of inertia of a hollow cylinder about its axis: I = ½m(R₁² + R₂²)."
  @spec hollow_cylinder_moment_of_inertia(number, number, number) :: float
  def hollow_cylinder_moment_of_inertia(mass, inner_radius, outer_radius) do
    0.5 * mass * (:math.pow(inner_radius, 2) + :math.pow(outer_radius, 2))
  end

  @doc "Parallel-axis theorem: I = I_cm + md², where d is distance from centre of mass."
  @spec parallel_axis(number, number, number) :: float
  def parallel_axis(i_cm, mass, d) do
    i_cm + mass * :math.pow(d, 2)
  end

  # ---------------------------------------------------------------------------
  # Projectile Motion
  # ---------------------------------------------------------------------------

  @doc "Projectile horizontal displacement: x = v_x t."
  @spec horizontal_displacement(number, number) :: float
  def horizontal_displacement(initial_vx, time), do: initial_vx * time * 1.0

  @doc "Projectile vertical displacement: y = v_y₀t + ½at²."
  @spec vertical_displacement(number, number, number) :: float
  def vertical_displacement(initial_vy, time, acceleration_y) do
    initial_vy * time + 0.5 * acceleration_y * :math.pow(time, 2)
  end

  @doc "Final vertical speed from v_y² = v_y₀² + 2aΔy."
  @spec vertical_velocity(number, number, number) :: float
  def vertical_velocity(initial_vy, acceleration_y, delta_y) do
    :math.sqrt(:math.pow(initial_vy, 2) + 2 * acceleration_y * delta_y)
  end

  @doc """
  Time of flight for projectile launched and landing at equal height: t = 2v₀ sin θ / g.

  ## Parameters
    - v0:    Launch speed (m/s)
    - theta: Launch angle (radians)
    - g:     g (m/s², default: 9.81)

  ## Examples
      iex> Motion.time_of_flight(20, :math.pi()/4) |> Float.round(3)
      2.887
  """
  @spec time_of_flight(number, number, number) :: float
  def time_of_flight(v0, theta, g \\ 9.81) do
    2 * v0 * :math.sin(theta) / g
  end

  @doc "Horizontal range of a projectile launched and landing at equal height: R = v₀² sin(2θ)/g."
  @spec projectile_range(number, number, number) :: float
  def projectile_range(v0, theta, g \\ 9.81) do
    :math.pow(v0, 2) * :math.sin(2 * theta) / g
  end

  @doc "Maximum height of a projectile: H = v₀² sin²θ/(2g)."
  @spec projectile_max_height(number, number, number) :: float
  def projectile_max_height(v0, theta, g \\ 9.81) do
    :math.pow(v0 * :math.sin(theta), 2) / (2 * g)
  end

  # ---------------------------------------------------------------------------
  # Simple Harmonic Motion
  # ---------------------------------------------------------------------------

  @doc "Simple harmonic displacement: x(t) = A cos(omegat + φ)."
  @spec shm_displacement(number, number, number, number) :: float
  def shm_displacement(amplitude, omega, t, phi \\ 0.0) do
    amplitude * :math.cos(omega * t + phi)
  end

  @doc "Simple harmonic velocity: v(t) = −Aomega sin(omegat + φ)."
  @spec shm_velocity(number, number, number, number) :: float
  def shm_velocity(amplitude, omega, t, phi \\ 0.0) do
    -amplitude * omega * :math.sin(omega * t + phi)
  end

  @doc "Simple harmonic acceleration: a(t) = -Aomega² cos(omegat + φ)."
  @spec shm_acceleration(number, number, number, number) :: float
  def shm_acceleration(amplitude, omega, t, phi \\ 0.0) do
    -amplitude * :math.pow(omega, 2) * :math.cos(omega * t + phi)
  end

  @doc "Angular frequency of a spring-mass oscillator: omega = √(k/m)."
  @spec spring_angular_frequency(number, number) :: float
  def spring_angular_frequency(k, mass), do: :math.sqrt(k / mass)

  @doc "Small-angle pendulum angular frequency: omega = √(g/L)."
  @spec pendulum_angular_frequency(number, number) :: float
  def pendulum_angular_frequency(length, g \\ 9.81), do: :math.sqrt(g / length)

  @doc "Oscillation period from angular frequency: T = 2π/omega."
  @spec period_from_omega(number) :: float
  def period_from_omega(omega), do: 2 * :math.pi() / omega

  # ---------------------------------------------------------------------------
  # Analytical Mechanics
  # ---------------------------------------------------------------------------

  @doc "Lagrangian of a mechanical system: L = T − V (kinetic minus potential energy)."
  @spec lagrangian(number, number) :: float
  def lagrangian(kinetic, potential), do: kinetic - potential

  @doc "Generalised canonical momentum via numerical differentiation of the Lagrangian: p = ∂L/∂q̇."
  @spec generalized_momentum((number -> number), number, number) :: float
  def generalized_momentum(lagrangian_func, q_dot, h \\ 1.0e-6) do
    (lagrangian_func.(q_dot + h) - lagrangian_func.(q_dot - h)) / (2 * h)
  end

  @doc "Hamiltonian (total energy in canonical coordinates): H = Σ pᵢq̇ᵢ − L."
  @spec hamiltonian([number], [number], number) :: float
  def hamiltonian(momenta, velocities, lagrangian) do
    Enum.zip_with(momenta, velocities, fn p, v -> p * v end)
    |> Enum.sum()
    |> Kernel.-(lagrangian)
    |> Kernel.*(1.0)
  end

  @doc "Hamilton's equations of motion for a 1-D harmonic oscillator: {dp/dt, dq/dt}."
  @spec hamiltons_equations(number, number, number) :: {float, float}
  def hamiltons_equations(p, q, omega) do
    {-:math.pow(omega, 2) * q, p * 1.0}
  end

  @doc "Poisson bracket {f, g}(q, p) computed by numerical partial differentiation."
  @spec poisson_bracket(
          (number, number -> number),
          (number, number -> number),
          number,
          number,
          number
        ) :: float
  def poisson_bracket(f, g, q, p, h \\ 1.0e-6) do
    df_dq = (f.(q + h, p) - f.(q - h, p)) / (2 * h)
    df_dp = (f.(q, p + h) - f.(q, p - h)) / (2 * h)
    dg_dq = (g.(q + h, p) - g.(q - h, p)) / (2 * h)
    dg_dp = (g.(q, p + h) - g.(q, p - h)) / (2 * h)
    df_dq * dg_dp - df_dp * dg_dq
  end

  # ---------------------------------------------------------------------------
  # Kepler's Laws
  # ---------------------------------------------------------------------------

  @doc """
  Kepler's Third Law: T = 2π √(a³ / (G M))

  ## Parameters
    - a:            Semi-major axis (m)
    - central_mass: M (kg)

  ## Examples
      iex> Motion.keplers_third_law(1.496e11, 1.989e30) > 0
      true
  """
  @spec keplers_third_law(number, number) :: float
  def keplers_third_law(a, central_mass) do
    2 * :math.pi() * :math.sqrt(:math.pow(a, 3) / (@gravitational_constant * central_mass))
  end

  @doc """
  Vis-viva equation: v = √(G M (2/r - 1/a))

  ## Parameters
    - central_mass: M (kg)
    - r:            Current orbital radius (m)
    - a:            Semi-major axis (m)

  ## Examples
      iex> Motion.vis_viva(5.972e24, 6.771e6, 6.771e6) > 0
      true
  """
  @spec vis_viva(number, number, number) :: float
  def vis_viva(central_mass, r, a) do
    :math.sqrt(@gravitational_constant * central_mass * (2 / r - 1 / a))
  end
end
