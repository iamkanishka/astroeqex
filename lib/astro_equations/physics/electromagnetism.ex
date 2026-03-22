defmodule AstroEquations.Physics.Electromagnetism do
  @moduledoc """
  Fundamental electromagnetic equations including:
  - Maxwell's equations (integral and differential forms)
  - Lorentz force
  - Electric field, potential, and energy
  - Dipole fields and moments
  - Charge and current densities
  - Circuit theory (Ohm's law, series/parallel, RC/RL/LC)
  - Capacitors and inductors
  - Magnetic fields (Biot-Savart, solenoid, toroid)
  - Materials (permittivity, susceptibility, polarisation, magnetisation)
  - Electromagnetic waves (speed, impedance, Poynting vector)
  - Skin depth, plasma frequency, Hall effect

  All calculations use SI units.
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Electric charge in coulombs (C)."
  @type charge :: float()

  @typedoc "Electric field magnitude in V/m."
  @type electric_field :: float()

  @typedoc "Magnetic field magnitude in tesla (T). Non-negative."
  @type magnetic_field :: float()

  @typedoc "Distance in metres (m). Must be positive."
  @type distance :: float()

  @typedoc "Electric potential in volts (V)."
  @type potential :: float()

  @typedoc "Current in amperes (A)."
  @type current :: float()

  @typedoc "Permittivity in F/m."
  @type permittivity :: float()

  @typedoc "Permeability in H/m."
  @type permeability :: float()

  # Vacuum permittivity ε₀ (F/m)
  @epsilon_0 8.854_187_812_8e-12
  # Vacuum permeability μ₀ (N/A²)
  @mu_0 1.256_637_062_12e-6
  # Elementary charge e (C)
  @elementary_charge 1.602_176_634e-19
  # Speed of light c (m/s)
  @speed_of_light 2.997_924_58e8

  # ---------------------------------------------------------------------------
  # Maxwell's Equations
  # ---------------------------------------------------------------------------

  @doc """
  Electric flux through a closed surface (Gauss's Law): Φ_E = Q_enc / ε₀

  ## Parameters
    - q_enc:     Enclosed charge (C)
    - epsilon_0: Permittivity of free space (default: ε₀)

  ## Examples
      iex> Electromagnetism.gauss_law(1.0) > 0
      true
  """

  @spec gauss_law(number, number) :: float
  def gauss_law(q_enc, epsilon_0 \\ @epsilon_0), do: q_enc / epsilon_0

  @doc "Differential form of Gauss's Law: ∇·E = ρ/ε₀."

  @spec gauss_law_differential(number, number) :: float
  def gauss_law_differential(rho, epsilon_0 \\ @epsilon_0), do: rho / epsilon_0

  @doc "Gauss's Law for magnetism: ∇·B = 0, reflecting the absence of magnetic monopoles. Always returns 0."

  @spec gauss_law_magnetism() :: 0
  def gauss_law_magnetism, do: 0

  @doc "Faraday's Law of electromagnetic induction (differential form): ∇×E = -∂B/∂t."

  @spec faraday_law(number) :: number()
  def faraday_law(dB_dt), do: -dB_dt

  @doc """
  Faraday induction EMF: ε = -dΦ_B/dt

  ## Parameters
    - d_phi_b_dt: Rate of change of magnetic flux (Wb/s)

  ## Examples
      iex> Electromagnetism.faraday_emf(0.05)
      -0.05
  """

  @spec faraday_emf(number) :: number()
  def faraday_emf(d_phi_b_dt), do: -d_phi_b_dt

  @doc """
  Ampère-Maxwell law: ∇×B = μ₀(J + ε₀ ∂E/∂t)

  ## Parameters
    - current_density: J (A/m²)
    - dE_dt:           ∂E/∂t (V/m·s)

  ## Examples
      iex> Electromagnetism.ampere_law(1.0e6, 0.0) > 0
      true
  """

  @spec ampere_law(number, number, number, number) :: number()
  def ampere_law(current_density, dE_dt, mu_0 \\ @mu_0, epsilon_0 \\ @epsilon_0) do
    mu_0 * (current_density + epsilon_0 * dE_dt)
  end

  @doc """
  Displacement current: I_d = ε₀ dΦ_E/dt

  The displacement current term added by Maxwell; ensures charge conservation
  and predicts electromagnetic waves.

  ## Parameters
    - d_phi_e_dt: Rate of change of electric flux Φ_E (V·m/s)
    - epsilon_0:  Permittivity of free space (default: ε₀)

  ## Examples
      iex> Electromagnetism.displacement_current(1.0e6) > 0
      true
  """

  @spec displacement_current(number, number) :: number()
  def displacement_current(d_phi_e_dt, epsilon_0 \\ @epsilon_0) do
    epsilon_0 * d_phi_e_dt
  end

  @doc """
  Electric flux through a flat surface: Φ_E = E A cos θ

  ## Parameters
    - e_field: Electric field magnitude (V/m)
    - area:    Surface area (m²)
    - theta:   Angle between field and surface normal (radians, default: 0)

  ## Examples
      iex> Electromagnetism.electric_flux(100.0, 0.5, 0)
      50.0
  """

  @spec electric_flux(number, number, number) :: float
  def electric_flux(e_field, area, theta \\ 0.0) do
    e_field * area * :math.cos(theta)
  end

  # ---------------------------------------------------------------------------
  # Lorentz Force
  # ---------------------------------------------------------------------------

  @doc """
  Lorentz force on a point charge: F = q(E + v × B)

  ## Parameters
    - q:        Charge (C)
    - e_field:  Electric field {Ex, Ey, Ez} (V/m)
    - velocity: Charge velocity {vx, vy, vz} (m/s)
    - b_field:  Magnetic field {Bx, By, Bz} (T)

  ## Returns
    Force vector {Fx, Fy, Fz} (N)

  ## Examples
      iex> e = {1.0, 0.0, 0.0}; v = {0.0, 0.0, 0.0}; b = {0.0, 0.0, 0.0}
      ...> Electromagnetism.lorentz_force_point(1.0, e, v, b) |> elem(0) |> Float.round(2)
      1.0
  """

  @spec lorentz_force_point(
          number,
          {number, number, number},
          {number, number, number},
          {number, number, number}
        ) :: {float, float, float}
  def lorentz_force_point(q, e_field, velocity, b_field) do
    {ex, ey, ez} = e_field
    {vx, vy, vz} = velocity
    {bx, by, bz} = b_field

    {q * (ex + vy * bz - vz * by), q * (ey + vz * bx - vx * bz), q * (ez + vx * by - vy * bx)}
  end

  @doc """
  Cyclotron radius for a charged particle in a magnetic field: r = m v_⊥ / (|q| B)

  ## Parameters
    - mass:    Particle mass (kg)
    - v_perp:  Speed perpendicular to B (m/s)
    - charge:  Particle charge magnitude (C)
    - b_field: Magnetic field magnitude (T)

  ## Examples
      iex> Electromagnetism.cyclotron_radius(9.109e-31, 1.0e6, 1.602e-19, 0.01) > 0
      true
  """

  @spec cyclotron_radius(number, number, number, number) :: float
  def cyclotron_radius(mass, v_perp, charge, b_field) when is_positive(mass) do
    mass * v_perp / (abs(charge) * b_field)
  end

  @doc """
  Cyclotron (gyro) frequency: omega_c = |q| B / m

  ## Parameters
    - charge:  Charge magnitude (C)
    - b_field: Magnetic field magnitude (T)
    - mass:    Particle mass (kg)

  ## Examples
      iex> Electromagnetism.cyclotron_frequency(1.602e-19, 1.0, 1.673e-27) > 0
      true
  """

  @spec cyclotron_frequency(number, number, number) :: float
  def cyclotron_frequency(charge, b_field, mass) when is_positive(mass),
    do: abs(charge) * b_field / mass

  # ---------------------------------------------------------------------------
  # Electric Field
  # ---------------------------------------------------------------------------

  @doc "Coulomb's constant: k_e = 1/(4πε₀) ≈ 8.988×10⁹ N·m²/C²."

  @spec coulombs_constant(number) :: float
  def coulombs_constant(epsilon_0 \\ @epsilon_0) do
    1.0 / (4 * :math.pi() * epsilon_0)
  end

  @doc """
  Electric field from a point charge: E = q / (4πε₀ r²)

  ## Examples
      iex> Electromagnetism.electric_field_point(1.0e-9, 1.0) |> Float.round(2)
      8.99
  """

  @spec electric_field_point(number, number, number) :: float
  def electric_field_point(q, r, epsilon_0 \\ @epsilon_0) do
    q / (4 * :math.pi() * epsilon_0 * :math.pow(r, 2))
  end

  @doc """
  Electric field on the axis of a uniformly charged ring:
  E = q z / (4πε₀ (z² + R²)^(3/2))

  ## Parameters
    - q:   Total ring charge (C)
    - z:   Axial distance from ring centre (m)
    - r:   Ring radius (m)

  ## Examples
      iex> Electromagnetism.electric_field_ring(1.0e-9, 0.1, 0.1) > 0
      true
  """

  @spec electric_field_ring(number, number, number, number) :: float
  def electric_field_ring(q, z, r, epsilon_0 \\ @epsilon_0) do
    q * z / (4 * :math.pi() * epsilon_0 * :math.pow(z * z + r * r, 1.5))
  end

  @doc "Electric field on the axis of an electric dipole: E = 2p/(4πε₀r³)."

  @spec dipole_field_parallel(number, number, number) :: float
  def dipole_field_parallel(p, r, epsilon_0 \\ @epsilon_0) do
    2 * p / (4 * :math.pi() * epsilon_0 * :math.pow(r, 3))
  end

  @doc "Electric field perpendicular to an electric dipole axis: E = p/(4πε₀r³)."

  @spec dipole_field_perpendicular(number, number, number) :: float
  def dipole_field_perpendicular(p, r, epsilon_0 \\ @epsilon_0) do
    p / (4 * :math.pi() * epsilon_0 * :math.pow(r, 3))
  end

  @doc "Electric dipole moment vector: p = q d (charge × displacement vector)."

  @spec dipole_moment(number, {number, number, number}) :: {float, float, float}
  def dipole_moment(q, {dx, dy, dz}), do: {q * dx, q * dy, q * dz}

  # ---------------------------------------------------------------------------
  # Electric Potential & Energy
  # ---------------------------------------------------------------------------

  @doc "Electric potential from a point charge: V = q/(4πε₀r)."

  @spec electric_potential_point(number, number, number) :: float
  def electric_potential_point(q, r, epsilon_0 \\ @epsilon_0) do
    q / (4 * :math.pi() * epsilon_0 * r)
  end

  @doc "Potential difference V(a) - V(b) between two radii from a point charge."

  @spec potential_difference_point(number, number, number, number) :: float
  def potential_difference_point(q, a, b, epsilon_0 \\ @epsilon_0) do
    1 / (4 * :math.pi() * epsilon_0) * q * (1 / b - 1 / a)
  end

  @doc "Electrostatic potential energy of two point charges: U = q₁q₂/(4πε₀r)."

  @spec potential_energy(number, number, number, number) :: float
  def potential_energy(q1, q2, r, epsilon_0 \\ @epsilon_0) do
    q1 * q2 / (4 * :math.pi() * epsilon_0 * r)
  end

  @doc "Energy density of an electric field: u_E = ½ε₀E²."

  @spec field_energy_density(number, number) :: float
  def field_energy_density(e_field, epsilon_0 \\ @epsilon_0) do
    0.5 * epsilon_0 * e_field * e_field
  end

  @doc "Total electric field energy over a volume: U = ½ε₀E² V."

  @spec field_energy(number, number, number) :: float
  def field_energy(e_field, volume, epsilon_0 \\ @epsilon_0) do
    0.5 * epsilon_0 * :math.pow(e_field, 2) * volume
  end

  @doc """
  Total electromagnetic energy density: u = ½ε₀E² + B²/(2μ₀)

  ## Parameters
    - e_field: Electric field magnitude (V/m)
    - b_field: Magnetic field magnitude (T)

  ## Examples
      iex> Electromagnetism.em_energy_density(1000.0, 3.33e-6) > 0
      true
  """

  @spec em_energy_density(number, number, number, number) :: float
  def em_energy_density(e_field, b_field, epsilon_0 \\ @epsilon_0, mu_0 \\ @mu_0) do
    0.5 * epsilon_0 * e_field * e_field + b_field * b_field / (2 * mu_0)
  end

  # ---------------------------------------------------------------------------
  # Charge & Current Densities
  # ---------------------------------------------------------------------------

  @doc "Surface charge density: σ = Q/A."

  @spec surface_charge_density(number, number) :: float
  def surface_charge_density(q, area), do: q / area

  @doc "Linear charge density along a wire: λ = Q/L."

  @spec linear_charge_density(number, number) :: float
  def linear_charge_density(q, length) when is_positive(length), do: q / length

  @doc "Volume current density: J = I/A."

  @spec current_density(number, number) :: float
  def current_density(current, area), do: current / area

  @doc "Electron drift velocity: v_d = μE, where μ is carrier mobility."

  @spec drift_velocity(number, number) :: number()
  def drift_velocity(mobility, e_field), do: mobility * e_field

  @doc "Electric current from charge carrier properties: I = nAqv_d."

  @spec current_from_properties(number, number, number, number, number) :: number()
  def current_from_properties(n, area, charge \\ @elementary_charge, mobility, e_field) do
    n * area * charge * mobility * e_field
  end

  @doc "Continuity equation (charge conservation): ∂ρ/∂t = -∇·J."

  @spec continuity(number) :: number()
  def continuity(div_j), do: -div_j

  # ---------------------------------------------------------------------------
  # Circuit Theory
  # ---------------------------------------------------------------------------

  @doc "Ohm's Law: V = IR."

  @spec ohms_law(number, number) :: number()
  def ohms_law(current, resistance), do: current * resistance

  @doc "Electrical power delivered to a component: P = IV."

  @spec electrical_power(number, number) :: number()
  def electrical_power(current, voltage), do: current * voltage

  @doc "Electrical power dissipated in a resistor: P = I²R."

  @spec electrical_power_from_resistance(number, number) :: float
  def electrical_power_from_resistance(current, resistance) do
    :math.pow(current, 2) * resistance
  end

  @doc "Electrical power from voltage and resistance: P = V²/R."

  @spec electrical_power_from_voltage(number, number) :: float
  def electrical_power_from_voltage(voltage, resistance) do
    :math.pow(voltage, 2) / resistance
  end

  @doc "Resistance of a conductor from its resistivity: R = ρL/A."

  @spec resistance(number, number, number) :: float
  def resistance(resistivity, length, area) when is_positive(length),
    do: resistivity * length / area

  @doc "Total resistance of resistors in series: R = Σ Rᵢ."

  @spec series_resistance([number]) :: number()
  def series_resistance(resistances), do: Enum.sum(resistances)

  @doc "Total resistance of resistors in parallel: 1/R = Σ 1/Rᵢ."

  @spec parallel_resistance([number]) :: float
  def parallel_resistance(resistances) do
    1 / Enum.reduce(resistances, 0, fn r, acc -> acc + 1 / r end)
  end

  @doc "Total capacitance of capacitors in series: 1/C = Σ 1/Cᵢ."

  @spec series_capacitance([number]) :: float
  def series_capacitance(capacitances) do
    1 / Enum.reduce(capacitances, 0, fn c, acc -> acc + 1 / c end)
  end

  @doc "Total capacitance of capacitors in parallel: C = Σ Cᵢ."

  @spec parallel_capacitance([number]) :: number()
  def parallel_capacitance(capacitances), do: Enum.sum(capacitances)

  @doc "Time constant of an RC circuit: τ = RC."

  @spec rc_time_constant(number, number) :: number()
  def rc_time_constant(resistance, capacitance), do: resistance * capacitance

  @doc "Time constant of an RL circuit: τ = L/R."

  @spec rl_time_constant(number, number) :: float
  def rl_time_constant(inductance, resistance), do: inductance / resistance

  @doc "Natural (resonant) angular frequency of an LC circuit: omega₀ = 1/√(LC)."

  @spec lc_resonant_frequency(number, number) :: float
  def lc_resonant_frequency(inductance, capacitance) do
    1.0 / :math.sqrt(inductance * capacitance)
  end

  @doc "Impedance magnitude of a capacitor: |Z_C| = 1/(ωC)."

  @spec capacitor_impedance(number, number) :: float
  def capacitor_impedance(omega, capacitance), do: 1.0 / (omega * capacitance)

  @doc "Impedance magnitude of an inductor: |Z_L| = ωL."

  @spec inductor_impedance(number, number) :: number()
  def inductor_impedance(omega, inductance), do: omega * inductance

  @doc """
  RLC series resonance Q-factor: Q = ω₀L/R = (1/R)√(L/C)

  ## Parameters
    - inductance:  L (H)
    - capacitance: C (F)
    - resistance:  R (Ω)

  ## Examples
      iex> Electromagnetism.rlc_quality_factor(1.0e-3, 10.0e-6, 0.1) > 0
      true
  """

  @spec rlc_quality_factor(number, number, number) :: float
  def rlc_quality_factor(inductance, capacitance, resistance) do
    1.0 / resistance * :math.sqrt(inductance / capacitance)
  end

  # ---------------------------------------------------------------------------
  # Capacitors
  # ---------------------------------------------------------------------------

  @doc "Capacitance from stored charge and applied voltage: C = Q/V."

  @spec capacitance(number, number) :: float
  def capacitance(charge, voltage), do: charge / voltage

  @doc "Parallel-plate capacitance: C = ε A/d."

  @spec capacitance_from_geometry(number, number, number) :: float
  def capacitance_from_geometry(epsilon, area, distance) when is_positive(distance),
    do: epsilon * area / distance

  @doc "Energy stored in a charged capacitor: U = ½CV²."

  @spec capacitor_energy(number, number) :: float
  def capacitor_energy(capacitance, voltage) do
    0.5 * capacitance * :math.pow(voltage, 2)
  end

  @doc "Electric field between the plates of a parallel-plate capacitor: E = Q/(ε₀A)."

  @spec capacitor_field(number, number, number) :: float
  def capacitor_field(charge, area, epsilon_0 \\ @epsilon_0), do: charge / (epsilon_0 * area)

  # ---------------------------------------------------------------------------
  # Magnetic Fields
  # ---------------------------------------------------------------------------

  @doc """
  Biot-Savart law — differential magnetic field dB from a current element:
  dB = (μ₀/4π) I (dl × r̂) / r²

  ## Parameters
    - current:  Current I (A)
    - dl:       Current element vector {dlx, dly, dlz} (m)
    - r_vector: Unit displacement vector {rx, ry, rz}
    - r_mag:    Distance magnitude (m)
  """

  @spec biot_savart(number, {number, number, number}, {number, number, number}, number) ::
          {float, float, float}
  def biot_savart(current, dl, r_vector, r_mag) do
    factor = @mu_0 / (4 * :math.pi()) * current / :math.pow(r_mag, 2)
    {cx, cy, cz} = cross3(dl, r_vector)
    {factor * cx, factor * cy, factor * cz}
  end

  @doc "Magnetic field produced by a moving point charge."

  @spec moving_charge_field(number, {number, number, number}, {number, number, number}, number) ::
          {float, float, float}
  def moving_charge_field(q, velocity, r_vector, r_mag) do
    factor = @mu_0 / (4 * :math.pi()) * q / :math.pow(r_mag, 2)
    {cx, cy, cz} = cross3(velocity, r_vector)
    {factor * cx, factor * cy, factor * cz}
  end

  @doc "Magnetic field surrounding an infinite straight current-carrying wire: B = μ₀I/(2πr)."

  @spec wire_magnetic_field(number, number, number) :: float
  def wire_magnetic_field(current, distance, mu_0 \\ @mu_0) when is_positive(distance) do
    mu_0 * current / (2 * :math.pi() * distance)
  end

  @doc """
  Magnetic field inside an ideal solenoid: B = μ₀ n I

  ## Parameters
    - n:       Turns per metre (m⁻¹)
    - current: Current (A)

  ## Examples
      iex> Electromagnetism.solenoid_field(1000, 2.0) > 0
      true
  """

  @spec solenoid_field(number, number, number) :: number()
  def solenoid_field(n, current, mu_0 \\ @mu_0), do: mu_0 * n * current

  @doc """
  Magnetic field inside a toroid: B = μ₀ N I / (2π r)

  ## Examples
      iex> Electromagnetism.toroid_field(100, 1.0, 0.1) > 0
      true
  """

  @spec toroid_field(number, number, number, number) :: float
  def toroid_field(n_turns, current, r, mu_0 \\ @mu_0) do
    mu_0 * n_turns * current / (2 * :math.pi() * r)
  end

  @doc "Force per unit length between two parallel current-carrying wires."

  @spec parallel_wire_force(number, number, number, number) :: float
  def parallel_wire_force(i1, i2, distance, mu_0 \\ @mu_0) when is_positive(distance) do
    mu_0 * i1 * i2 / (2 * :math.pi() * distance)
  end

  @doc "Magnetic flux through a surface: Φ = BA cos θ."

  @spec magnetic_flux(number, number, number) :: float
  def magnetic_flux(b_field, area, theta \\ 0.0) do
    b_field * area * :math.cos(theta)
  end

  # ---------------------------------------------------------------------------
  # Inductors
  # ---------------------------------------------------------------------------

  @doc "Self-inductance of an ideal solenoid: L = μ₀n²V."

  @spec solenoid_inductance(number, number, number) :: number()
  def solenoid_inductance(n, volume, mu_0 \\ @mu_0), do: mu_0 * n * n * volume

  @doc """
  Self-inductance from flux linkage: L = NΦ / I

  ## Parameters
    - n_turns: Number of turns N
    - flux:    Magnetic flux Φ per turn (Wb)
    - current: Current I (A)

  ## Examples
      iex> Electromagnetism.self_inductance_from_flux(100, 5.0e-4, 2.0) > 0
      true
  """

  @spec self_inductance_from_flux(number, number, number) :: float
  def self_inductance_from_flux(n_turns, flux, current) do
    n_turns * flux / current
  end

  @doc "EMF induced in an inductor: ε = -L dI/dt."

  @spec inductor_emf(number, number) :: number()
  def inductor_emf(inductance, dI_dt), do: -inductance * dI_dt

  @doc "Energy stored in an inductor carrying current I: U = ½LI²."

  @spec inductor_energy(number, number) :: float
  def inductor_energy(inductance, current) do
    0.5 * inductance * :math.pow(current, 2)
  end

  @doc "Energy density of a magnetic field: u = B²/(2μ₀)."

  @spec magnetic_energy_density(number, number) :: float
  def magnetic_energy_density(b_field, mu_0 \\ @mu_0) do
    b_field * b_field / (2 * mu_0)
  end

  @doc "EMF induced in a secondary coil by mutual inductance: ε₂ = -M dI₁/dt."

  @spec mutual_inductance_emf(number, number) :: number()
  def mutual_inductance_emf(m, dI1_dt), do: -m * dI1_dt

  @doc "Magnetic vector potential (far-field approximation): A ≈ (μ₀/4π) J/r."

  @spec vector_potential({number, number, number}, number, number) :: {float, float, float}
  def vector_potential({jx, jy, jz}, distance, mu_0 \\ @mu_0) when is_positive(distance) do
    factor = mu_0 / (4 * :math.pi()) / distance
    {factor * jx, factor * jy, factor * jz}
  end

  # ---------------------------------------------------------------------------
  # Materials
  # ---------------------------------------------------------------------------

  @doc "Relative permittivity (dielectric constant) of a material: εᵣ = ε/ε₀."

  @spec relative_permittivity(number, number) :: float
  def relative_permittivity(permittivity, vacuum_permittivity \\ @epsilon_0) do
    permittivity / vacuum_permittivity
  end

  @doc "Electric susceptibility: χₑ = εᵣ - 1."

  @spec electric_susceptibility(number) :: number()
  def electric_susceptibility(relative_permittivity), do: relative_permittivity - 1

  @doc "Absolute permittivity of a medium: ε = εᵣε₀."

  @spec absolute_permittivity(number, number) :: number()
  def absolute_permittivity(relative_permittivity, vacuum_permittivity \\ @epsilon_0) do
    relative_permittivity * vacuum_permittivity
  end

  @doc "Electric polarisation vector: P = ε₀χₑE (or P = n·p for discrete dipoles)."

  @spec polarization(number, [number], number | nil, [number] | nil) :: [float]
  def polarization(chi_e, electric_field, n \\ nil, dipole_moment \\ nil) do
    if n && dipole_moment do
      Enum.map(dipole_moment, fn p -> n * p end)
    else
      Enum.map(electric_field, fn e -> @epsilon_0 * chi_e * e end)
    end
  end

  @doc "Surface bound charge density: σ_b = P·n̂."

  @spec surface_bound_charge([number], [number]) :: number()
  def surface_bound_charge(polarization, normal_vector) do
    Enum.zip(polarization, normal_vector)
    |> Enum.map(fn {p, n} -> p * n end)
    |> Enum.sum()
  end

  @doc "Volume bound charge density (finite-difference approximation of −∇·P)."

  @spec volume_bound_charge([[number]], number) :: float
  def volume_bound_charge([p1, p2], delta_x) do
    delta_p = Enum.zip(p1, p2) |> Enum.map(fn {a, b} -> b - a end)
    -Enum.sum(delta_p) / delta_x
  end

  @doc "Total bound charge from surface and volume contributions: Q_b = σ_b·A + ρ_b·V."

  @spec total_bound_charge(number, number) :: number()
  def total_bound_charge(surface_charge, volume_charge), do: surface_charge + volume_charge

  @doc "Electric displacement field: D = ε₀E + P."

  @spec electric_displacement(number, [number], [number]) :: [float]
  def electric_displacement(permittivity, electric_field, polarization \\ [0, 0, 0]) do
    d = Enum.map(electric_field, fn e -> permittivity * e end)
    Enum.zip(d, polarization) |> Enum.map(fn {di, pi} -> di + pi end)
  end

  @doc "Magnetic field strength: H = B/μ₀ − M."

  @spec magnetic_field_strength([number], [number]) :: [float]
  def magnetic_field_strength(magnetic_flux, magnetization \\ [0, 0, 0]) do
    h = Enum.map(magnetic_flux, fn b -> b / @mu_0 end)
    Enum.zip(h, magnetization) |> Enum.map(fn {hi, mi} -> hi - mi end)
  end

  @doc "Magnetic susceptibility: χ_m = μᵣ − 1."

  @spec magnetic_susceptibility(number) :: number()
  def magnetic_susceptibility(relative_permeability), do: relative_permeability - 1

  @doc "Relative permeability of a material: μᵣ = μ/μ₀."

  @spec relative_permeability(number, number) :: float
  def relative_permeability(permeability, mu_0 \\ @mu_0), do: permeability / mu_0

  @doc "Magnetic dipole moment vector: m = I·A (current times area vector)."

  @spec magnetic_dipole_moment(number, [number]) :: [float]
  def magnetic_dipole_moment(current, area_vector) do
    Enum.map(area_vector, fn a -> current * a end)
  end

  @doc "Bound surface current density: K_b = M × n̂."

  @spec bound_surface_current([number], [number]) :: [float]
  def bound_surface_current([mx, my, mz], [nx, ny, nz]) do
    [my * nz - mz * ny, mz * nx - mx * nz, mx * ny - my * nx]
  end

  @doc "Bound charge from polarisation (integral form): Q_b = -∮P·dA."

  @spec gauss_law_polarization([number], [number]) :: number()
  def gauss_law_polarization(polarization_vector, surface_area) do
    dot =
      Enum.zip(polarization_vector, surface_area)
      |> Enum.map(fn {p, da} -> p * da end)
      |> Enum.sum()

    -dot
  end

  @doc "Free charge from displacement field (integral form): Q_f = ∮D·dA."

  @spec gauss_law_displacement([number], [number]) :: number()
  def gauss_law_displacement(displacement_vector, surface_area) do
    Enum.zip(displacement_vector, surface_area)
    |> Enum.map(fn {d, da} -> d * da end)
    |> Enum.sum()
  end

  # ---------------------------------------------------------------------------
  # Electromagnetic Waves
  # ---------------------------------------------------------------------------

  @doc """
  Speed of light in a medium: v = 1 / √(με)

  ## Examples
      iex> Electromagnetism.wave_speed(1.257e-6, 8.854e-12) |> round()
      299792030
  """

  @spec wave_speed(number, number) :: float
  def wave_speed(mu, epsilon), do: 1.0 / :math.sqrt(mu * epsilon)

  @doc """
  Characteristic impedance of free space: Z₀ = √(μ₀/ε₀) ≈ 377 Ω

  ## Examples
      iex> Electromagnetism.free_space_impedance() |> Float.round(2)
      376.73
  """

  @spec free_space_impedance(number, number) :: float
  def free_space_impedance(mu_0 \\ @mu_0, epsilon_0 \\ @epsilon_0) do
    :math.sqrt(mu_0 / epsilon_0)
  end

  @doc """
  Poynting vector magnitude: S = E × H = E B / μ₀

  ## Examples
      iex> Electromagnetism.poynting_magnitude(1000, 3.0e-4) > 0
      true
  """

  @spec poynting_magnitude(number, number, number) :: float
  def poynting_magnitude(e_field, b_field, mu_0 \\ @mu_0) do
    e_field * b_field / mu_0
  end

  @doc """
  Radiation pressure for perfect absorption: P_rad = S / c

  ## Examples
      iex> Electromagnetism.radiation_pressure(1000) > 0
      true
  """

  @spec radiation_pressure(number, number) :: float
  def radiation_pressure(intensity, c \\ @speed_of_light), do: intensity / c

  @doc """
  Skin depth in a conductor: δ = √(2 / (μ σ ω))

  ## Examples
      iex> Electromagnetism.skin_depth(1.257e-6, 5.8e7, 6_283_185.0) > 0
      true
  """

  @spec skin_depth(number, number, number) :: float
  def skin_depth(mu, sigma, omega), do: :math.sqrt(2.0 / (mu * sigma * omega))

  @doc """
  Plasma frequency: ω_p = √(n e² / (ε₀ m_e))

  ## Examples
      iex> Electromagnetism.plasma_frequency(1.0e18) > 0
      true
  """

  @spec plasma_frequency(number, number, number) :: float
  def plasma_frequency(n, m_e \\ 9.10938e-31, epsilon_0 \\ @epsilon_0) do
    :math.sqrt(n * :math.pow(@elementary_charge, 2) / (epsilon_0 * m_e))
  end

  @doc """
  Brewster's angle for polarisation by reflection: θ_B = arctan(n₂/n₁)

  ## Examples
      iex> Electromagnetism.brewsters_angle(1.0, 1.5) > 0
      true
  """

  @spec brewsters_angle(number, number) :: float
  def brewsters_angle(n1, n2), do: :math.atan2(n2, n1)

  @doc """
  Hall voltage across a conductor in a magnetic field: V_H = IB / (n q t)

  Used in astronomical plasma diagnostics and semiconductor characterisation.

  ## Parameters
    - current:        I (A)
    - b_field:        Magnetic field B (T)
    - carrier_density: n (m⁻³)
    - charge:         Carrier charge q (C, default: elementary charge)
    - thickness:      Sample thickness t in field direction (m)

  ## Examples
      iex> Electromagnetism.hall_voltage(1.0, 0.5, 1.0e28, 1.602e-19, 1.0e-3) > 0
      true
  """

  @spec hall_voltage(number, number, number, number, number) :: float

  def hall_voltage(current, b_field, carrier_density, charge \\ @elementary_charge, thickness) do
    current * b_field / (carrier_density * charge * thickness)
  end

  # 3D cross product: a × b
  defp cross3({ax, ay, az}, {bx, by, bz}) do
    {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}
  end
end
