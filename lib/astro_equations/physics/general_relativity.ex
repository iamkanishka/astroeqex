defmodule AstroEquations.Physics.GeneralRelativity do
  @moduledoc """
  General Relativity: spacetime metrics, geodesics, curvature, and gravitational phenomena.

  Covers:
  - Spacetime metrics (Minkowski, Schwarzschild, Kerr equatorial, Rindler, FLRW)
  - Tensor algebra (index raising/lowering, inner product, inverse metric)
  - Christoffel symbols (Schwarzschild)
  - Geodesic equations (radial free-fall)
  - Gravitational time dilation and redshift
  - Gravitational wave strain
  - Orbital precession (Mercury / post-Newtonian)
  - Light deflection by a mass
  - Cosmology (Friedmann equation, critical density, Hubble parameter)
  - Schwarzschild radius, Hawking temperature, photon sphere

  Natural units c = G = 1 are used unless physical constants are explicit parameters.
  """

  @gravitational_constant 6.67430e-11
  @speed_of_light 2.99792458e8
  # Reduced Planck constant (J·s)
  @hbar 1.054571817e-34
  # Boltzmann constant (J/K)
  @boltzmann 1.380649e-23

  # ---------------------------------------------------------------------------
  # Metrics
  # ---------------------------------------------------------------------------

  @doc """
  Returns the Minkowski metric tensor η (signature -+++)

  ## Examples
      iex> GeneralRelativity.minkowski_metric()
      [[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
  """
  @spec minkowski_metric() :: [[integer]]
  def minkowski_metric do
    [[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
  end

  @doc """
  Minkowski spacetime interval: ds² = -c²dt² + dx² + dy² + dz²

  ## Examples
      iex> GeneralRelativity.minkowski_interval(1, 0, 0, 0)
      -1
  """
  @spec minkowski_interval(number, number, number, number, number) :: float
  def minkowski_interval(dt, dx, dy, dz, c \\ 1) do
    -(c * c * dt * dt) + dx * dx + dy * dy + dz * dz
  end

  @doc """
  Schwarzschild metric tensor g_μν for a spherically symmetric mass.

  ## Parameters
    - m:     Mass (in units where G = c = 1)
    - r:     Radial coordinate
    - theta: Polar angle
    - c, g:  Speed of light and G (default: 1 each)

  ## Examples
      iex> GeneralRelativity.schwarzschild_metric(1, 3, :math.pi()/2) |> length()
      4
  """
  @spec schwarzschild_metric(number, number, number, number, number) :: [[number]]
  def schwarzschild_metric(m, r, theta, c \\ 1, g \\ 1) do
    rs = 2 * g * m / (c * c)
    factor = 1 - rs / r

    [
      [-factor, 0, 0, 0],
      [0, 1 / factor, 0, 0],
      [0, 0, r * r, 0],
      [0, 0, 0, r * r * :math.sin(theta) ** 2]
    ]
  end

  @doc """
  Schwarzschild spacetime interval ds².

  ## Examples
      iex> GeneralRelativity.schwarzschild_interval(1, 3, :math.pi()/2, 1, 0, 0, 0) |> Float.round(4)
      -0.3333
  """
  @spec schwarzschild_interval(
          number,
          number,
          number,
          number,
          number,
          number,
          number,
          number,
          number
        ) :: float
  def schwarzschild_interval(m, r, theta, dt, dr, dtheta, dphi, c \\ 1, g \\ 1) do
    rs = 2 * g * m / (c * c)
    factor = 1 - rs / r

    -factor * c * c * dt * dt +
      1 / factor * dr * dr +
      r * r * dtheta * dtheta +
      r * r * :math.sin(theta) ** 2 * dphi * dphi
  end

  @doc """
  Rindler spacetime interval (uniformly accelerating frame):
  ds² = -(1 + g x'/c²)² c² dt'² + dx'²

  ## Examples
      iex> GeneralRelativity.rindler_interval(9.8, 1, 1, 0) < 0
      true
  """
  @spec rindler_interval(number, number, number, number, number) :: float
  def rindler_interval(g, x_prime, dt_prime, dx_prime, c \\ 1) do
    -(1 + g * x_prime / (c * c)) ** 2 * c * c * dt_prime * dt_prime +
      dx_prime * dx_prime
  end

  @doc """
  Flat FLRW spacetime interval: ds² = -c²dt² + a(t)²(dx² + dy² + dz²)

  ## Examples
      iex> GeneralRelativity.flrw_interval(1, 1, 0, 1)
      -1
  """
  @spec flrw_interval(number, number, number, number) :: float
  def flrw_interval(dt, a, dr, c \\ 1) do
    -(c * c * dt * dt) + a * a * dr * dr
  end

  # ---------------------------------------------------------------------------
  # Tensor Algebra
  # ---------------------------------------------------------------------------

  @doc "Lowers an index of a 4-vector using the metric: V_μ = g_μν V^ν."
  @spec lower_index([number], [[number]], [atom], non_neg_integer) :: [number]
  def lower_index(tensor, metric, _index_positions, index_to_lower) do
    tensor
    |> Enum.with_index()
    |> Enum.map(fn {component, i} ->
      if i == index_to_lower do
        Enum.reduce(0..3, 0, fn j, acc ->
          acc + Enum.at(Enum.at(metric, i), j) * component
        end)
      else
        component
      end
    end)
  end

  @doc """
  Raises an index using the inverse metric: V^μ = g^μν V_ν.

  For diagonal metrics the inverse is simply 1/g_μμ on the diagonal.
  """
  @spec raise_index([number], [[number]], [atom], non_neg_integer) :: [number]
  def raise_index(tensor, metric, _index_positions, index_to_raise) do
    inv = inverse_metric(metric)

    tensor
    |> Enum.with_index()
    |> Enum.map(fn {component, i} ->
      if i == index_to_raise do
        Enum.reduce(0..3, 0, fn j, acc ->
          acc + Enum.at(Enum.at(inv, i), j) * component
        end)
      else
        component
      end
    end)
  end

  @doc "Transforms a contravariant vector under a Jacobian: V'^i = (∂x'^i/∂x^j) V^j."
  @spec transform_tensor([number], [[number]], [[number]]) :: [number]
  def transform_tensor(tensor, jacobian, _inverse_jacobian) when is_list(tensor) do
    if not is_list(List.first(tensor)) do
      Enum.map(0..3, fn i ->
        Enum.reduce(0..3, 0, fn j, acc ->
          acc + Enum.at(tensor, j) * Enum.at(Enum.at(jacobian, j), i)
        end)
      end)
    else
      tensor
    end
  end

  @doc "Minkowski (or general metric) inner product of two four-vectors: a·b = g_μν aμ bν."
  @spec four_vector_product([number], [number], [[number]]) :: float
  def four_vector_product(a, b, metric) do
    Enum.reduce(0..3, 0, fn i, acc1 ->
      Enum.reduce(0..3, acc1, fn j, acc2 ->
        acc2 + Enum.at(Enum.at(metric, i), j) * Enum.at(a, i) * Enum.at(b, j)
      end)
    end)
  end

  @doc """
  Inverse of a diagonal metric tensor: g^μμ = 1/g_μμ.

  For non-diagonal metrics use a proper matrix inversion library.
  """
  @spec inverse_metric([[number]]) :: [[number]]
  def inverse_metric(metric) do
    Enum.map(metric, fn row ->
      Enum.map(row, fn
        0 -> 0
        x -> 1 / x
      end)
    end)
  end

  # ---------------------------------------------------------------------------
  # Gravitational Effects
  # ---------------------------------------------------------------------------

  @doc """
  Gravitational time dilation: dt_local / dt_inf = √(1 - r_s / r)

  ## Returns
    Time dilation factor (approaches 0 at event horizon, 1 at infinity)

  ## Examples
      iex> GeneralRelativity.gravitational_time_dilation(1.989e30, 1.0e15) |> Float.round(6)
      1.0
  """
  @spec gravitational_time_dilation(number, number) :: float
  def gravitational_time_dilation(mass, radius) do
    rs = 2 * @gravitational_constant * mass / @speed_of_light ** 2
    :math.sqrt(max(1 - rs / radius, 0.0))
  end

  @doc """
  Gravitational redshift: z = 1/√(1 - r_s/r) - 1

  ## Examples
      iex> GeneralRelativity.gravitational_redshift(1.989e30, 6.957e8) > 0
      true
  """
  @spec gravitational_redshift(number, number) :: float
  def gravitational_redshift(mass, r_emit) do
    rs = 2 * @gravitational_constant * mass / @speed_of_light ** 2
    1.0 / :math.sqrt(max(1 - rs / r_emit, 1.0e-15)) - 1
  end

  @doc """
  Schwarzschild radius: r_s = 2 G M / c²

  The radius at which the escape velocity equals the speed of light.

  ## Examples
      iex> GeneralRelativity.schwarzschild_radius(1.989e30) > 0
      true
  """
  @spec schwarzschild_radius(number) :: float
  def schwarzschild_radius(mass) do
    2 * @gravitational_constant * mass / @speed_of_light ** 2
  end

  @doc """
  Hawking temperature of a Schwarzschild black hole: T_H = ħ c³ / (8π G M k_B)

  ## Parameters
    - mass: Black hole mass (kg)

  ## Returns
    Hawking temperature in Kelvin

  ## Examples
      iex> GeneralRelativity.hawking_temperature(1.989e30) > 0
      true
  """
  @spec hawking_temperature(number) :: float
  def hawking_temperature(mass) do
    @hbar * @speed_of_light ** 3 /
      (8 * :math.pi() * @gravitational_constant * mass * @boltzmann)
  end

  @doc """
  Photon sphere radius for a Schwarzschild black hole: r_ph = 3 G M / c² = 3r_s/2

  ## Examples
      iex> GeneralRelativity.photon_sphere_radius(1.989e30) > 0
      true
  """
  @spec photon_sphere_radius(number) :: float
  def photon_sphere_radius(mass) do
    3 * @gravitational_constant * mass / @speed_of_light ** 2
  end

  @doc """
  Gravitational wave strain amplitude (quadrupole approximation):
  h ~ (2G / c⁴) × (d²Q/dt²) / r

  ## Examples
      iex> GeneralRelativity.gravitational_wave_strain(1.0e47, 1.0e25) > 0
      true
  """
  @spec gravitational_wave_strain(number, number) :: float
  def gravitational_wave_strain(d2Q_dt2, distance) do
    2 * @gravitational_constant / @speed_of_light ** 4 * d2Q_dt2 / distance
  end

  @doc """
  Orbital precession per orbit (Schwarzschild geodesic / post-Newtonian):
  Δφ = 6π G M / (a (1 - e²) c²)

  ## Examples
      iex> GeneralRelativity.orbital_precession(1.989e30, 5.79e10, 0.206) > 0
      true
  """
  @spec orbital_precession(number, number, number) :: float
  def orbital_precession(mass, a, e) do
    6 * :math.pi() * @gravitational_constant * mass /
      (a * (1 - e * e) * @speed_of_light ** 2)
  end

  @doc """
  Light deflection by a point mass (GR prediction): α = 4 G M / (b c²)

  ## Examples
      iex> GeneralRelativity.light_deflection(1.989e30, 6.957e8) > 0
      true
  """
  @spec light_deflection(number, number) :: float
  def light_deflection(mass, impact_param) do
    4 * @gravitational_constant * mass / (impact_param * @speed_of_light ** 2)
  end

  @doc """
  Radial free-fall coordinate velocity in Schwarzschild spacetime:
  dr/dt = -(1 - r_s/r) √(r_s/r)  (from rest at infinity)

  ## Returns
    dr/dt in m/s (negative = inward)

  ## Examples
      iex> GeneralRelativity.radial_freefall_velocity(1.989e30, 1.0e9) < 0
      true
  """
  @spec radial_freefall_velocity(number, number) :: float
  def radial_freefall_velocity(mass, r) do
    rs = 2 * @gravitational_constant * mass / @speed_of_light ** 2

    if r <= rs,
      do: 0.0,
      else: -(1 - rs / r) * @speed_of_light * :math.sqrt(rs / r)
  end

  # ---------------------------------------------------------------------------
  # Cosmology (Friedmann Equations)
  # ---------------------------------------------------------------------------

  @doc """
  Friedmann equation (flat universe): H = √(8πG ρ / 3)

  Returns the Hubble parameter from energy density.

  ## Examples
      iex> GeneralRelativity.hubble_from_density(9.47e-27) > 0
      true
  """
  @spec hubble_from_density(number) :: float
  def hubble_from_density(rho) do
    :math.sqrt(8 * :math.pi() * @gravitational_constant * rho / 3)
  end

  @doc """
  Critical density for a flat universe: ρ_c = 3 H² / (8π G)

  ## Examples
      iex> GeneralRelativity.critical_density(2.27e-18) > 0
      true
  """
  @spec critical_density(number) :: float
  def critical_density(h) do
    3 * h * h / (8 * :math.pi() * @gravitational_constant)
  end

  @doc """
  Density parameter: Ω = ρ / ρ_c

  ## Examples
      iex> GeneralRelativity.density_parameter(9.47e-27, 9.47e-27) |> Float.round(4)
      1.0
  """
  @spec density_parameter(number, number) :: float
  def density_parameter(rho, rho_c), do: rho / rho_c

  @doc """
  Scale factor in a matter-dominated flat universe: a(t) = a₀ (t/t₀)^(2/3)

  ## Examples
      iex> GeneralRelativity.scale_factor_matter(1, 1)
      1.0
  """
  @spec scale_factor_matter(number, number, number) :: float
  def scale_factor_matter(t, t0, a0 \\ 1.0) do
    a0 * :math.pow(t / t0, 2 / 3)
  end

  @doc """
  Lookback time approximation (flat matter-dominated universe):
  t_lb ≈ (2/3) (1/H₀) [1 - 1/(1+z)^(3/2)]

  ## Examples
      iex> GeneralRelativity.lookback_time(0, 2.27e-18)
      0.0
  """
  @spec lookback_time(number, number) :: float
  def lookback_time(z, h0) do
    2.0 / 3.0 / h0 * (1 - 1 / :math.pow(1 + z, 1.5))
  end
end
