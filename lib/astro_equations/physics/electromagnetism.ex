defmodule AstroEquations.Physics.Electromagnetism do
  @moduledoc """
  A module containing implementations of fundamental electromagnetic equations including:
  - Maxwell's equations
  - Lorentz force
  - Electric field calculations
  - Dipole fields and moments
  - Electric potential and energy
  - Charge and current densities
  - Circuit theory
  - Capacitor
  - Magnetic fields
  - Inductors

  All calculations use SI units.
  """

  # Vacuum permittivity (F/m)
  @epsilon_0 8.8541878128e-12
  # Vacuum permeability (N/A^2)
  @mu_0 1.25663706212e-6
  # Elementary charge (C)
  @elementary_charge 1.602176634e-19
  # Pi constant
  @pi :math.pi()

  @doc """
  Calculates the electric flux through a closed surface using Gauss's Law.

  ## Parameters
    - q_enc: enclosed charge (Coulombs)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Electric flux (N·m²/C)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.gauss_law(1.0)
      1.1294090673382828e11
  """
  def gauss_law(q_enc, epsilon_0 \\ @epsilon_0) do
    q_enc / epsilon_0
  end

  @doc """
  Calculates the divergence of the electric field at a point with given charge density.

  ## Parameters
    - rho: charge density (C/m³)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Divergence of E field (V/m²)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.gauss_law_differential(1.0)
      1.1294090673382828e11
  """
  def gauss_law_differential(rho, epsilon_0 \\ @epsilon_0) do
    rho / epsilon_0
  end

  @doc """
  Calculates the magnetic flux through a closed surface (always zero according to Gauss's Law for Magnetism).

  ## Returns
    - Always returns 0 (no magnetic monopoles)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.gauss_law_magnetism()
      0
  """
  def gauss_law_magnetism, do: 0

  @doc """
  Calculates the curl of the electric field from a changing magnetic field (Faraday's Law of Induction).

  ## Parameters
    - dB_dt: time derivative of magnetic field (T/s)

  ## Returns
    - Curl of E field (V/m²)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.faraday_law(1.0)
      -1.0
  """
  def faraday_law(dB_dt), do: -dB_dt

  @doc """
  Calculates the curl of the magnetic field from current and changing electric field (Ampère's circuital law with Maxwell's correction).

  ## Parameters
    - current_density: current density J (A/m²)
    - dE_dt: time derivative of electric field (V/m·s)
    - mu_0: permeability of free space (defaults to @mu_0)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Curl of B field (T/m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.ampere_law(1.0, 1.0)
      1.25663706212e-6
  """
  def ampere_law(current_density, dE_dt, mu_0 \\ @mu_0, epsilon_0 \\ @epsilon_0) do
    mu_0 * (current_density + epsilon_0 * dE_dt)
  end

  @doc """
  Calculates the Lorentz force on a point charge.

  ## Parameters
    - q: charge (Coulombs)
    - e_field: electric field vector (V/m)
    - velocity: velocity vector of charge (m/s)
    - b_field: magnetic field vector (Tesla)

  ## Returns
    - Force vector (Newtons)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.lorentz_force_point(1.0, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0})
      {1.0, 0.0, -1.0}
  """
  def lorentz_force_point(q, e_field, velocity, b_field) do
    {ex, ey, ez} = e_field
    {vx, vy, vz} = velocity
    {bx, by, bz} = b_field

    # Calculate v × B
    cross_x = vy * bz - vz * by
    cross_y = vz * bx - vx * bz
    cross_z = vx * by - vy * bx

    # Calculate F = q(E + v × B)
    {
      q * (ex + cross_x),
      q * (ey + cross_y),
      q * (ez + cross_z)
    }
  end

  @doc """
  Calculates the electric field from a point charge.

  ## Parameters
    - q: charge (Coulombs)
    - r: distance from charge (meters)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Electric field magnitude (V/m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.electric_field_point(1.0, 1.0)
      8.987551787368176e9
  """
  def electric_field_point(q, r, epsilon_0 \\ @epsilon_0) do
    1 / (4 * :math.pi() * epsilon_0) * (q / :math.pow(r, 2))
  end

  @doc """
  Calculates the electric field magnitude parallel to the dipole axis at a distance r.

  ## Parameters
    - p: dipole moment (C·m)
    - r: distance from dipole (m)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Electric field magnitude parallel to dipole axis (V/m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.dipole_field_parallel(1e-12, 0.1)
      1.7975103574736352e-8
  """
  def dipole_field_parallel(p, r, epsilon_0 \\ @epsilon_0) do
    2 * p / (4 * :math.pi() * epsilon_0 * :math.pow(r, 3))
  end

  @doc """
  Calculates the electric field magnitude perpendicular to the dipole axis at a distance r.

  ## Parameters
    - p: dipole moment (C·m)
    - r: distance from dipole (m)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Electric field magnitude perpendicular to dipole axis (V/m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.dipole_field_perpendicular(1e-12, 0.1)
      8.987551787368176e-9
  """
  def dipole_field_perpendicular(p, r, epsilon_0 \\ @epsilon_0) do
    p / (4 * :math.pi() * epsilon_0 * :math.pow(r, 3))
  end

  @doc """
  Calculates the dipole moment from charge and separation distance.

  ## Parameters
    - q: charge magnitude (C)
    - d: separation distance vector (m)

  ## Returns
    - Dipole moment vector (C·m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.dipole_moment(1e-6, {1e-3, 0, 0})
      {1.0e-9, 0.0, 0.0}
  """
  def dipole_moment(q, {dx, dy, dz}) do
    {q * dx, q * dy, q * dz}
  end

  @doc """
  Calculates the electric potential from a point charge.

  ## Parameters
    - q: charge (C)
    - r: distance from charge (m)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Electric potential (V)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.electric_potential_point(1e-9, 0.1)
      89.87551787368176
  """
  def electric_potential_point(q, r, epsilon_0 \\ @epsilon_0) do
    1 / (4 * :math.pi() * epsilon_0) * (q / r)
  end

  @doc """
  Calculates the electric potential difference between two points in a field.

  ## Parameters
    - q: charge creating the field (C)
    - a: first distance from charge (m)
    - b: second distance from charge (m)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Potential difference (V)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.potential_difference_point(1e-9, 0.1, 0.2)
      44.93775893684088
  """
  def potential_difference_point(q, a, b, epsilon_0 \\ @epsilon_0) do
    1 / (4 * :math.pi() * epsilon_0) * q * (1 / b - 1 / a)
  end

  @doc """
  Calculates the electric potential energy between two charges.

  ## Parameters
    - q: first charge (C)
    - Q: second charge (C)
    - r: separation distance (m)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Potential energy (J)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.potential_energy(1e-9, 1e-9, 0.1)
      8.987551787368176e-8
  """
  def potential_energy(first_charge, second_charge, r, epsilon_0 \\ @epsilon_0) do
    1 / (4 * :math.pi() * epsilon_0) * (first_charge * second_charge / r)
  end

  @doc """
  Calculates the energy stored in an electrostatic field distribution.

  ## Parameters
    - e_field: electric field magnitude (V/m)
    - volume: volume of space (m³)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Stored energy (J)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.field_energy(1000, 0.001)
      4.4270939064e-6
  """
  def field_energy(e_field, volume, epsilon_0 \\ @epsilon_0) do
    0.5 * epsilon_0 * :math.pow(e_field, 2) * volume
  end

  @doc """
  Calculates surface charge density.

  ## Parameters
    - q: charge (C)
    - area: surface area (m²)

  ## Returns
    - Surface charge density (C/m²)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.surface_charge_density(1e-6, 0.01)
      1.0e-4
  """
  def surface_charge_density(q, area), do: q / area

  @doc """
  Calculates linear charge density.

  ## Parameters
    - q: charge (C)
    - length: length (m)

  ## Returns
    - Linear charge density (C/m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.linear_charge_density(1e-6, 0.1)
      1.0e-5
  """
  def linear_charge_density(q, length), do: q / length

  @doc """
  Calculates volume current density from current and perpendicular area.

  ## Parameters
    - current: electric current (A)
    - area: perpendicular cross-sectional area (m²)

  ## Returns
    - Current density magnitude (A/m²)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.current_density(1.0, 0.0001)
      10000.0
  """
  def current_density(current, area), do: current / area

  @doc """
  Calculates electron drift velocity in a conductor.

  ## Parameters
    - mobility: electron mobility (m²/(V·s))
    - e_field: electric field (V/m)

  ## Returns
    - Drift velocity (m/s)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.drift_velocity(0.001, 100)
      0.1
  """
  def drift_velocity(mobility, e_field), do: mobility * e_field

  @doc """
  Calculates current from charge carrier properties.

  ## Parameters
    - n: charge carrier density (1/m³)
    - area: cross-sectional area (m²)
    - charge: carrier charge (C, defaults to elementary charge)
    - mobility: electron mobility (m²/(V·s))
    - e_field: electric field (V/m)

  ## Returns
    - Current (A)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.current_from_properties(1e28, 1e-6, 1.6e-19, 0.001, 100)
      1.6
  """
  def current_from_properties(n, area, charge \\ @elementary_charge, mobility, e_field) do
    n * area * charge * mobility * e_field
  end

  @doc """
  Calculates electrical power from current and voltage.

  ## Parameters
    - current: electric current (A)
    - voltage: potential difference (V)

  ## Returns
    - Power (W)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.electrical_power(2, 10)
      20
  """
  def electrical_power(current, voltage), do: current * voltage

  @doc """
  Calculates electrical power from current and resistance.

  ## Parameters
    - current: electric current (A)
    - resistance: electrical resistance (Ω)

  ## Returns
    - Power (W)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.electrical_power_from_resistance(2, 5)
      20
  """
  def electrical_power_from_resistance(current, resistance),
    do: :math.pow(current, 2) * resistance

  @doc """
  Calculates voltage from current and resistance (Ohm's Law).

  ## Parameters
    - current: electric current (A)
    - resistance: electrical resistance (Ω)

  ## Returns
    - Voltage (V)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.ohms_law(2, 5)
      10
  """
  def ohms_law(current, resistance), do: current * resistance

  @doc """
  Calculates resistance from resistivity.

  ## Parameters
    - resistivity: material resistivity (Ω·m)
    - length: conductor length (m)
    - area: cross-sectional area (m²)

  ## Returns
    - Resistance (Ω)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.resistance(1.68e-8, 1, 1e-6) # Copper wire 1m long, 1mm² cross-section
      0.0168
  """
  def resistance(resistivity, length, area), do: resistivity * length / area

  @doc """
  Calculates series resistance for multiple resistors.

  ## Parameters
    - resistances: list of resistances (Ω)

  ## Returns
    - Total series resistance (Ω)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.series_resistance([10, 20, 30])
      60
  """
  def series_resistance(resistances), do: Enum.sum(resistances)

  @doc """
  Calculates parallel resistance for multiple resistors.

  ## Parameters
    - resistances: list of resistances (Ω)

  ## Returns
    - Total parallel resistance (Ω)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.parallel_resistance([10, 10])
      5.0
  """
  def parallel_resistance(resistances) do
    sum_reciprocals = Enum.reduce(resistances, 0, fn r, acc -> acc + 1 / r end)
    1 / sum_reciprocals
  end

  @doc """
  Calculates capacitance from charge and voltage.

  ## Parameters
    - charge: stored charge (C)
    - voltage: potential difference (V)

  ## Returns
    - Capacitance (F)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.capacitance(1e-6, 10)
      1.0e-7
  """
  def capacitance(charge, voltage), do: charge / voltage

  @doc """
  Calculates capacitance from physical properties.

  ## Parameters
    - epsilon: permittivity of dielectric (F/m)
    - area: plate area (m²)
    - distance: separation distance (m)

  ## Returns
    - Capacitance (F)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.capacitance_from_geometry(8.854e-12, 0.01, 1e-3)
      8.854e-11
  """
  def capacitance_from_geometry(epsilon, area, distance), do: epsilon * area / distance

  @doc """
  Calculates energy stored in a capacitor.

  ## Parameters
    - capacitance: capacitance (F)
    - voltage: potential difference (V)

  ## Returns
    - Stored energy (J)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.capacitor_energy(1e-6, 10)
      5.0e-5
  """
  def capacitor_energy(capacitance, voltage), do: 0.5 * capacitance * :math.pow(voltage, 2)

  @doc """
  Calculates electric field in a parallel plate capacitor.

  ## Parameters
    - charge: plate charge (C)
    - area: plate area (m²)
    - epsilon_0: permittivity of free space (defaults to @epsilon_0)

  ## Returns
    - Electric field (V/m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.capacitor_field(1e-6, 0.01)
      1.129409067518456e7
  """
  def capacitor_field(charge, area, epsilon_0 \\ @epsilon_0), do: charge / (epsilon_0 * area)

  @doc """
  Calculates the magnetic field from a current element (Biot-Savart Law).

  ## Parameters
    - current: current (A)
    - dl: length element vector (m)
    - r_vector: displacement vector from current to point (m)
    - r_mag: magnitude of displacement (m)

  ## Returns
    - Magnetic field vector (T)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.biot_savart(1.0, {0.0, 0.0, 1e-3}, {1.0, 0.0, 0.0}, 1.0)
      {0.0, 1.0e-10, 0.0}
  """
  def biot_savart(current, dl, r_vector, r_mag) do
    {dlx, dly, dlz} = dl
    {rx, ry, rz} = r_vector

    # Calculate dl × r̂
    cross_x = dly * rz - dlz * ry
    cross_y = dlz * rx - dlx * rz
    cross_z = dlx * ry - dly * rx

    # Calculate dB = (μ₀/4π) * I * (dl × r̂)/r²
    factor = @mu_0 / (4 * @pi) * current / :math.pow(r_mag, 2)

    {
      factor * cross_x,
      factor * cross_y,
      factor * cross_z
    }
  end

  @doc """
  Calculates the magnetic field from a moving charge (Biot-Savart for point charge).

  ## Parameters
    - q: charge (C)
    - velocity: velocity vector (m/s)
    - r_vector: displacement vector from charge to point (m)
    - r_mag: magnitude of displacement (m)

  ## Returns
    - Magnetic field vector (T)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.moving_charge_field(1.6e-19, {0.0, 1.0e6, 0.0}, {1.0, 0.0, 0.0}, 1.0)
      {0.0, 0.0, -1.6e-19 / :math.pow(1.0, 2) * 1.25663706212e-6 / (4 * :math.pi())}
  """
  def moving_charge_field(q, velocity, r_vector, r_mag) do
    {vx, vy, vz} = velocity
    {rx, ry, rz} = r_vector

    # Calculate v × r̂
    cross_x = vy * rz - vz * ry
    cross_y = vz * rx - vx * rz
    cross_z = vx * ry - vy * rx

    # Calculate B = (μ₀/4π) * q * (v × r̂)/r²
    factor = @mu_0 / (4 * @pi) * q / :math.pow(r_mag, 2)

    {
      factor * cross_x,
      factor * cross_y,
      factor * cross_z
    }
  end

  @doc """
  Calculates the magnetic field around a long straight wire.

  ## Parameters
    - current: current in wire (A)
    - distance: perpendicular distance from wire (m)
    - mu_0: permeability of free space (defaults to @mu_0)

  ## Returns
    - Magnetic field magnitude (T)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.wire_magnetic_field(1.0, 0.1)
      2.0e-6
  """
  def wire_magnetic_field(current, distance, mu_0 \\ @mu_0) do
    mu_0 / (2 * @pi) * current / distance
  end

  @doc """
  Calculates the EMF induced in an inductor.

  ## Parameters
    - inductance: inductance (H)
    - dI_dt: rate of change of current (A/s)

  ## Returns
    - Induced EMF (V)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.inductor_emf(0.1, 10.0)
      -1.0
  """
  def inductor_emf(inductance, dI_dt), do: -inductance * dI_dt

  @doc """
  Calculates the energy stored in an inductor.

  ## Parameters
    - inductance: inductance (H)
    - current: current (A)

  ## Returns
    - Stored energy (J)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.inductor_energy(0.1, 2.0)
      0.2
  """
  def inductor_energy(inductance, current), do: 0.5 * inductance * :math.pow(current, 2)

  @doc """
  Calculates the magnetic vector potential from current density.

  ## Parameters
    - current_density: current density vector (A/m²)
    - distance: distance from current element (m)
    - mu_0: permeability of free space (defaults to @mu_0)

  ## Returns
    - Vector potential (T·m)

  ## Examples
      iex> AstroEquations.Physics.Electromagnetism.vector_potential({1.0, 0.0, 0.0}, 1.0)
      {1.25663706212e-7, 0.0, 0.0}
  """
  def vector_potential({jx, jy, jz}, distance, mu_0 \\ @mu_0) do
    factor = mu_0 / (4 * @pi) / distance

    {
      factor * jx,
      factor * jy,
      factor * jz
    }
  end
end
