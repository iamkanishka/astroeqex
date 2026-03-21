defmodule AstroEquations.Physics.NewtonGravity do
  @moduledoc """
  Newtonian gravitational mechanics.

  Covers:
  - Gravitational force, field, potential, and potential energy
  - Near-surface approximation (mgh)
  - Orbital mechanics (Kepler's Third Law, orbital speed, escape velocity)
  - Tidal forces and Roche limit
  - Hill sphere
  - Shell theorem
  - Surface gravity
  - Vis-viva equation
  - Orbital energy and angular momentum
  - Two-body reduced mass
  - Lagrange L1 point estimate
  """

  # m³ kg⁻¹ s⁻²
  @gravitational_constant 6.67430e-11

  # ---------------------------------------------------------------------------
  # Core Gravitational Equations
  # ---------------------------------------------------------------------------

  @doc """
  Newton's Law of Universal Gravitation: F = G m₁ m₂ / r²

  ## Parameters
    - m1, m2: Masses (kg)
    - r:      Separation distance (m)

  ## Examples
      iex> NewtonGravity.force(5.972e24, 7.348e22, 3.844e8)
      1.9817963747937598e20
  """
  @spec force(number, number, number) :: float
  def force(m1, m2, r) when is_number(m1) and is_number(m2) and is_number(r) and r > 0 do
    @gravitational_constant * m1 * m2 / :math.pow(r, 2)
  end

  def force(_, _, _), do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Gravitational field strength (acceleration): g = G M / r²

  ## Examples
      iex> NewtonGravity.field(5.972e24, 6.371e6) |> Float.round(2)
      9.82
  """
  @spec field(number, number) :: float
  def field(mass, r) when r > 0, do: @gravitational_constant * mass / :math.pow(r, 2)
  def field(_, _), do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Gravitational potential: Φ = -G M / r

  ## Examples
      iex> NewtonGravity.potential(5.972e24, 6.371e6) |> Float.round(0)
      -62565146.0
  """
  @spec potential(number, number) :: float
  def potential(mass, r) when r > 0, do: -@gravitational_constant * mass / r
  def potential(_, _), do: raise(ArgumentError, "Distance must be positive")

  @doc """
  Gravitational potential energy between two masses: U = -G m₁ m₂ / r

  ## Examples
      iex> NewtonGravity.potential_energy(5.972e24, 7.348e22, 3.844e8) < 0
      true
  """
  @spec potential_energy(number, number, number) :: float
  def potential_energy(m1, m2, r) when r > 0 do
    -@gravitational_constant * m1 * m2 / r
  end

  def potential_energy(_, _, _), do: raise(ArgumentError, "Distance must be positive")

  @doc "Near-surface gravitational potential energy (uniform field): U = mgh."
  @spec approximate_potential_energy(number, number, number) :: float
  def approximate_potential_energy(m, g, h), do: m * g * h * 1.0

  # ---------------------------------------------------------------------------
  # Orbital Mechanics
  # ---------------------------------------------------------------------------

  @doc """
  Kepler's Third Law (two-body): T = 2π √(a³ / (G (m₁ + m₂)))

  ## Examples
      iex> NewtonGravity.orbital_period(3.844e8, 5.972e24, 7.348e22) |> Float.round(0)
      2360449.0
  """
  @spec orbital_period(number, number, number) :: float
  def orbital_period(a, m1, m2) when a > 0 do
    :math.sqrt(
      4 * :math.pow(:math.pi(), 2) * :math.pow(a, 3) /
        (@gravitational_constant * (m1 + m2))
    )
  end

  def orbital_period(_, _, _), do: raise(ArgumentError, "Semi-major axis must be positive")

  @doc "Circular orbit speed: v_c = √(GM/r)."
  @spec circular_orbit_speed(number, number) :: float
  def circular_orbit_speed(mass, r) when r > 0 do
    :math.sqrt(@gravitational_constant * mass / r)
  end

  @doc """
  Escape velocity: v_esc = √(2 G M / r)

  ## Examples
      iex> NewtonGravity.escape_velocity(5.972e24, 6.371e6) |> Float.round(0)
      11186.0
  """
  @spec escape_velocity(number, number) :: float
  def escape_velocity(mass, r) when r > 0 do
    :math.sqrt(2 * @gravitational_constant * mass / r)
  end

  @doc """
  Vis-viva orbital speed: v = √(G M (2/r - 1/a))

  ## Parameters
    - central_mass: M (kg)
    - r:            Current distance (m)
    - a:            Semi-major axis (m)

  ## Examples
      iex> NewtonGravity.vis_viva(5.972e24, 6.771e6, 6.771e6) > 0
      true
  """
  @spec vis_viva(number, number, number) :: float
  def vis_viva(central_mass, r, a) do
    :math.sqrt(@gravitational_constant * central_mass * (2 / r - 1 / a))
  end

  @doc """
  Specific orbital energy: ε = -G M / (2 a)

  ## Examples
      iex> NewtonGravity.specific_orbital_energy(5.972e24, 6.771e6) < 0
      true
  """
  @spec specific_orbital_energy(number, number) :: float
  def specific_orbital_energy(central_mass, a) do
    -@gravitational_constant * central_mass / (2 * a)
  end

  @doc """
  Specific angular momentum of a circular orbit: h = √(G M r)

  ## Examples
      iex> NewtonGravity.specific_angular_momentum(5.972e24, 6.771e6) > 0
      true
  """
  @spec specific_angular_momentum(number, number) :: float
  def specific_angular_momentum(central_mass, r) do
    :math.sqrt(@gravitational_constant * central_mass * r)
  end

  @doc """
  Semi-major axis from period (Kepler III inversion): a = (G M T² / 4π²)^(1/3)

  ## Examples
      iex> NewtonGravity.semi_major_axis_from_period(5.972e24, 5400) > 0
      true
  """
  @spec semi_major_axis_from_period(number, number) :: float
  def semi_major_axis_from_period(central_mass, period) do
    :math.pow(
      @gravitational_constant * central_mass * :math.pow(period, 2) /
        (4 * :math.pow(:math.pi(), 2)),
      1 / 3
    )
  end

  @doc """
  Two-body reduced mass: μ = m₁ m₂ / (m₁ + m₂)

  Used to convert a two-body problem into an equivalent one-body problem.

  ## Examples
      iex> NewtonGravity.reduced_mass(1.0, 1.0) |> Float.round(4)
      0.5
  """
  @spec reduced_mass(number, number) :: float
  def reduced_mass(m1, m2), do: m1 * m2 / (m1 + m2)

  @doc """
  Approximate distance of the L1 Lagrange point from the secondary body:
  r_L1 ≈ a (m₂ / (3 m₁))^(1/3)

  Valid when m₂ ≪ m₁.

  ## Parameters
    - a:   Orbital semi-major axis (m)
    - m1:  Primary mass (kg)
    - m2:  Secondary mass (kg)

  ## Examples
      iex> NewtonGravity.lagrange_l1(1.496e11, 1.989e30, 5.972e24) > 0
      true
  """
  @spec lagrange_l1(number, number, number) :: float
  def lagrange_l1(a, m1, m2) do
    a * :math.pow(m2 / (3 * m1), 1 / 3)
  end

  # ---------------------------------------------------------------------------
  # Surface Gravity
  # ---------------------------------------------------------------------------

  @doc """
  Surface gravity of a body: g = G M / R²

  ## Examples
      iex> NewtonGravity.surface_gravity(5.972e24, 6.371e6) |> Float.round(2)
      9.82
  """
  @spec surface_gravity(number, number) :: float
  def surface_gravity(mass, radius) do
    @gravitational_constant * mass / :math.pow(radius, 2)
  end

  # ---------------------------------------------------------------------------
  # Tidal Forces
  # ---------------------------------------------------------------------------

  @doc """
  Tidal acceleration across an object of size d at distance r from mass M:
  Δa ≈ 2 G M d / r³

  ## Examples
      iex> NewtonGravity.tidal_acceleration(7.348e22, 3.844e8, 6.371e6) > 0
      true
  """
  @spec tidal_acceleration(number, number, number) :: float
  def tidal_acceleration(mass, r, d) do
    2 * @gravitational_constant * mass * d / :math.pow(r, 3)
  end

  @doc """
  Roche limit (fluid body): d_R = a (2 M_primary / M_secondary)^(1/3)

  ## Examples
      iex> NewtonGravity.roche_limit(1.737e6, 5.972e24, 7.348e22) > 0
      true
  """
  @spec roche_limit(number, number, number) :: float
  def roche_limit(a, m_primary, m_secondary) do
    a * :math.pow(2 * m_primary / m_secondary, 1 / 3)
  end

  @doc """
  Hill sphere radius: r_H ≈ a (m / (3 M))^(1/3)

  ## Examples
      iex> NewtonGravity.hill_sphere(1.496e11, 5.972e24, 1.989e30) > 0
      true
  """
  @spec hill_sphere(number, number, number) :: float
  def hill_sphere(a, m, big_m) do
    a * :math.pow(m / (3 * big_m), 1 / 3)
  end
end
