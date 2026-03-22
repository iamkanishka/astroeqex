defmodule AstroEquations.Physics.Thermodynamics do
  @moduledoc """
  Thermodynamics and statistical mechanics.

  Covers:
  - Ideal gas law (all variable forms)
  - Heat energy, heat capacity, specific heat
  - First Law of Thermodynamics
  - Work done by gas (isothermal, adiabatic, isobaric)
  - Second Law — entropy and Clausius inequality
  - Carnot efficiency and COP
  - Boltzmann entropy (Ω microstates)
  - Thermodynamic potentials (enthalpy, Helmholtz, Gibbs)
  - Equipartition theorem
  - Maxwell-Boltzmann speed distribution
  - Black-body radiation (Planck, Wien, Stefan-Boltzmann)
  - Photon energy
  - Stellar luminosity (Stefan-Boltzmann)
  - Specific heat of gases (Cp, Cv, γ relation)
  - Newton's law of cooling
  - Van der Waals equation of state
  - Virial temperature (gravitational systems)
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Thermodynamic temperature in kelvin (K). Must be positive."
  @type temperature :: float()

  @typedoc "Pressure in pascals (Pa). Must be positive."
  @type pressure :: float()

  @typedoc "Volume in cubic metres (m³). Must be positive."
  @type volume :: float()

  @typedoc "Energy in joules (J)."
  @type energy :: float()

  @typedoc "Entropy in joules per kelvin (J/K)."
  @type entropy :: float()

  # J/K
  @boltzmann 1.380_649e-23
  # J·s
  @planck 6.626_070_15e-34
  # m/s
  @speed_light 299_792_458
  # m·K
  @wien_b 2.897_771_955e-3
  # W/(m²·K⁴)
  @stefan_sigma 5.670_374_419e-8
  # kg
  @proton_mass 1.672_621_923_69e-27
  # m³ kg⁻¹ s⁻²
  @gravitational_constant 6.674_30e-11

  # ---------------------------------------------------------------------------
  # Ideal Gas Law
  # ---------------------------------------------------------------------------

  @doc """
  Solves the particle-level ideal gas law pV = Nk_BT for the missing variable
  (pass nil for the unknown).

  ## Examples
      iex> Thermodynamics.ideal_gas_law(p: 101_325, v: 0.0224, n: 1, t: nil)
      %{t: _}
  """

  @spec ideal_gas_law(Keyword.t()) :: map()
  def ideal_gas_law(p: nil, v: v, n: n, k_b: k_b, t: t), do: %{p: n * k_b * t / v}
  def ideal_gas_law(p: p, v: nil, n: n, k_b: k_b, t: t), do: %{v: n * k_b * t / p}
  def ideal_gas_law(p: p, v: v, n: nil, k_b: k_b, t: t), do: %{n: p * v / (k_b * t)}
  def ideal_gas_law(p: p, v: v, n: n, k_b: k_b, t: nil), do: %{t: p * v / (n * k_b)}
  def ideal_gas_law(p: p, v: v, n: n, k_b: nil, t: t), do: %{k_b: p * v / (n * t)}

  def ideal_gas_law(p: p, v: v, n: n, t: t) do
    ideal_gas_law(p: p, v: v, n: n, k_b: @boltzmann, t: t)
  end

  @doc """
  Molar ideal gas law residual: |pV - nRT| (should be near zero for a real gas obeying this law).

  ## Parameters
    - p:     Pressure (Pa)
    - v:     Volume (m³)
    - n_mol: Amount (mol)
    - t:     Temperature (K)
    - r:     Gas constant (J/(mol·K), default: 8.31446)

  ## Returns
    %{check: float} — absolute deviation from ideal behaviour
  """

  @spec molar_ideal_gas_law(number | nil, number | nil, number | nil, number | nil, number) ::
          map()
  def molar_ideal_gas_law(p, v, n_mol, t, r \\ 8.31446) do
    %{check: abs(p * v - n_mol * r * t)}
  end

  @doc "Mean translational kinetic energy per molecule in an ideal gas: ⟨KE⟩ = 3k_BT/2."

  @spec mean_kinetic_energy(number, number) :: float
  def mean_kinetic_energy(temperature, k_b \\ @boltzmann) when is_positive(temperature) do
    1.5 * k_b * temperature
  end

  @doc """
  Equipartition theorem: mean energy per degree of freedom = k_BT/2.

  Total mean energy for f degrees of freedom: E = (f/2) k_B T.

  ## Parameters
    - degrees_of_freedom: f (e.g. 3 translational, 5 diatomic, 6 for full rotation+vibration)
    - temperature:        T (K)
    - k_b:                Boltzmann constant (default: k_B)

  ## Examples
      iex> Thermodynamics.equipartition_energy(3, 300) > 0
      true
  """

  @spec equipartition_energy(number, number, number) :: float
  def equipartition_energy(degrees_of_freedom, temperature, k_b \\ @boltzmann)
      when is_positive(temperature) do
    degrees_of_freedom / 2 * k_b * temperature
  end

  @doc "Root-mean-square speed of molecules in an ideal gas: v_rms = √(3k_BT/m)."

  @spec rms_speed(number, number, number) :: float
  def rms_speed(temperature, mass, k_b \\ @boltzmann)
      when is_positive(temperature) and is_positive(mass) do
    :math.sqrt(3 * k_b * temperature / mass)
  end

  @doc "Most probable speed in the Maxwell-Boltzmann distribution: v_p = √(2k_BT/m)."

  @spec most_probable_speed(number, number, number) :: float
  def most_probable_speed(temperature, mass, k_b \\ @boltzmann) when is_positive(temperature) do
    :math.sqrt(2 * k_b * temperature / mass)
  end

  @doc "Mean (average) molecular speed: ⟨v⟩ = √(8k_BT/(πm))."

  @spec mean_speed(number, number, number) :: float
  def mean_speed(temperature, mass, k_b \\ @boltzmann) when is_positive(temperature) do
    :math.sqrt(8 * k_b * temperature / (:math.pi() * mass))
  end

  # ---------------------------------------------------------------------------
  # Heat and Work
  # ---------------------------------------------------------------------------

  @doc "Heat transferred to raise a mass by ΔT: Q = mcΔT."

  @spec heat_energy(number, number, number) :: number()
  def heat_energy(m, c, delta_t), do: m * c * delta_t

  @doc "Heat capacity from the slope of the Q-T curve: C = dQ/dT."

  @spec heat_capacity(number, number) :: float
  def heat_capacity(dq, dt), do: dq / dt

  @doc "Specific heat capacity: c = C/m."

  @spec specific_heat_capacity(number, number) :: float
  def specific_heat_capacity(c_heat, m), do: c_heat / m

  @doc """
  First Law of Thermodynamics: ΔU = Q - W

  ## Examples
      iex> Thermodynamics.first_law(500, 200)
      300.0
  """

  @spec first_law(number, number) :: number()
  def first_law(q, w), do: q - w

  @doc "Work done by a gas at constant pressure: W = PΔV."

  @spec isobaric_work(number, number) :: number()
  def isobaric_work(pressure, delta_v), do: pressure * delta_v

  @doc "Work done by an ideal gas in an isothermal expansion: W = Nk_BT ln(V₂/V₁)."

  @spec isothermal_work(number, number, number, number, number) :: float
  def isothermal_work(n, temperature, v1, v2, k_b \\ @boltzmann) when is_positive(temperature) do
    n * k_b * temperature * :math.log(v2 / v1)
  end

  @doc "Work done by a gas in a reversible adiabatic process: W = (P₁V₁ - P₂V₂)/(γ-1)."

  @spec adiabatic_work(number, number, number, number, number) :: float
  def adiabatic_work(p1, v1, p2, v2, gamma) do
    (p1 * v1 - p2 * v2) / (gamma - 1)
  end

  @doc "Pressure after adiabatic compression/expansion: P₂ = P₁(V₁/V₂)^γ."

  @spec adiabatic_pressure(number, number, number, number) :: float
  def adiabatic_pressure(p1, v1, v2, gamma) do
    p1 * :math.pow(v1 / v2, gamma)
  end

  @doc "Temperature after adiabatic compression/expansion: T₂ = T₁(V₁/V₂)^(γ-1)."

  @spec adiabatic_temperature(number, number, number, number) :: float
  def adiabatic_temperature(t1, v1, v2, gamma) do
    t1 * :math.pow(v1 / v2, gamma - 1)
  end

  # ---------------------------------------------------------------------------
  # Thermodynamic Potentials
  # ---------------------------------------------------------------------------

  @doc """
  Enthalpy: H = U + PV

  ## Parameters
    - internal_energy: U (J)
    - pressure:        P (Pa)
    - volume:          V (m³)

  ## Examples
      iex> Thermodynamics.enthalpy(1000.0, 101_325.0, 0.001) |> Float.round(2)
      1101.325
  """

  @spec enthalpy(number, number, number) :: number()
  def enthalpy(internal_energy, pressure, volume) do
    internal_energy + pressure * volume
  end

  @doc """
  Helmholtz free energy: F = U - TS

  The maximum work extractable from a system at constant temperature.

  ## Parameters
    - internal_energy: U (J)
    - temperature:     T (K)
    - entropy:         S (J/K)

  ## Examples
      iex> Thermodynamics.helmholtz_free_energy(1000.0, 300.0, 2.0) |> Float.round(1)
      400.0
  """

  @spec helmholtz_free_energy(number, number, number) :: number()
  def helmholtz_free_energy(internal_energy, temperature, entropy)
      when is_positive(temperature) do
    internal_energy - temperature * entropy
  end

  @doc """
  Gibbs free energy: G = H - TS = U + PV - TS

  Determines spontaneity of a process at constant T and P; G < 0 is spontaneous.

  ## Parameters
    - internal_energy: U (J)
    - pressure:        P (Pa)
    - volume:          V (m³)
    - temperature:     T (K)
    - entropy:         S (J/K)

  ## Examples
      iex> Thermodynamics.gibbs_free_energy(1000.0, 101_325.0, 0.001, 300.0, 2.0) < 1000.0
      true
  """

  @spec gibbs_free_energy(number, number, number, number, number) :: number()
  def gibbs_free_energy(internal_energy, pressure, volume, temperature, entropy)
      when is_positive(temperature) do
    internal_energy + pressure * volume - temperature * entropy
  end

  # ---------------------------------------------------------------------------
  # Entropy & Second Law
  # ---------------------------------------------------------------------------

  @doc "Boltzmann entropy in terms of microstates: S = k_B ln Ω."

  @spec entropy(number, number) :: float
  def entropy(omega, k_b \\ @boltzmann), do: k_b * :math.log(omega)

  @doc "Clausius entropy change for a reversible process: ΔS = Q_rev/T."

  @spec entropy_change(number, number) :: float
  def entropy_change(q_rev, temperature) when is_positive(temperature), do: q_rev / temperature

  @doc """
  Number of microstates for an Einstein solid: Ω = (q+N-1)! / (q! (N-1)!)

  Note: uses integer arithmetic — this will overflow for large q or N.
  For large values, use the Stirling approximation instead.
  """

  @spec microstates(non_neg_integer, pos_integer) :: pos_integer
  def microstates(q, n) do
    div(factorial(q + n - 1), factorial(q) * factorial(n - 1))
  end

  @doc """
  Stirling approximation for ln Ω of an Einstein solid (large-N safe):
  ln Ω ≈ (q+N) ln(q+N) - q ln(q) - N ln(N)

  ## Parameters
    - q: Number of energy quanta
    - n: Number of oscillators

  ## Examples
      iex> Thermodynamics.microstates_stirling(100, 100) > 0
      true
  """

  @spec microstates_stirling(number, number) :: float
  def microstates_stirling(q, n) do
    qn = q + n
    qn * :math.log(qn) - q * :math.log(max(q, 1)) - n * :math.log(n)
  end

  # ---------------------------------------------------------------------------
  # Heat Engines & Carnot
  # ---------------------------------------------------------------------------

  @doc """
  Carnot efficiency: η = 1 - T_cold / T_hot

  ## Examples
      iex> Thermodynamics.carnot_efficiency(500, 300) |> Float.round(4)
      0.4
  """

  @spec carnot_efficiency(number, number) :: float
  def carnot_efficiency(t_hot, t_cold) when is_positive(t_hot) and is_positive(t_cold),
    do: 1 - t_cold / t_hot

  @doc "Coefficient of performance of a Carnot refrigerator: COP = T_cold/(T_hot − T_cold)."

  @spec cop_refrigerator(number, number) :: float
  def cop_refrigerator(t_hot, t_cold), do: t_cold / (t_hot - t_cold)

  @doc "Coefficient of performance of a Carnot heat pump: COP = T_hot/(T_hot − T_cold)."

  @spec cop_heat_pump(number, number) :: float
  def cop_heat_pump(t_hot, t_cold), do: t_hot / (t_hot - t_cold)

  # ---------------------------------------------------------------------------
  # Specific Heats (Ideal Gases)
  # ---------------------------------------------------------------------------

  @doc "Mayer's relation between molar heat capacities: Cₚ = Cᵥ + R."

  @spec mayers_relation(number, number) :: number()
  def mayers_relation(cv, r \\ 8.31446), do: cv + r

  @doc "Ratio of molar heat capacities: γ = Cₚ/Cᵥ."

  @spec heat_capacity_ratio(number, number) :: float
  def heat_capacity_ratio(cp, cv), do: cp / cv

  @doc "Molar isochoric heat capacity for a monatomic ideal gas: Cᵥ = 3R/2."

  @spec cv_monatomic(number) :: float
  def cv_monatomic(r \\ 8.31446), do: 1.5 * r

  @doc "Molar isochoric heat capacity for a diatomic ideal gas (room temperature): Cᵥ = 5R/2."

  @spec cv_diatomic(number) :: float
  def cv_diatomic(r \\ 8.31446), do: 2.5 * r

  # ---------------------------------------------------------------------------
  # Newton's Law of Cooling
  # ---------------------------------------------------------------------------

  @doc """
  Newton's law of cooling: T(t) = T_env + (T₀ - T_env) exp(-k t)

  ## Examples
      iex> Thermodynamics.newton_cooling(100, 20, 0.1, 0) |> Float.round(4)
      100.0
  """

  @spec newton_cooling(number, number, number, number) :: float
  def newton_cooling(t0, t_env, k, t) do
    t_env + (t0 - t_env) * :math.exp(-k * t)
  end

  # ---------------------------------------------------------------------------
  # Black-Body Radiation
  # ---------------------------------------------------------------------------

  @doc "Photon energy: E = hf."

  @spec photon_energy(number, number) :: number()
  def photon_energy(f, h \\ @planck), do: h * f

  @doc "Wien's displacement law for the peak emission wavelength: λ_max = b/T."

  @spec wiens_displacement(number, number) :: float
  def wiens_displacement(t, b \\ @wien_b) when is_positive(t), do: b / t

  @doc "Radiant flux density from a blackbody surface: F = σT⁴."

  @spec stefan_boltzmann(number, number) :: float
  def stefan_boltzmann(t, sigma \\ @stefan_sigma) when is_positive(t), do: sigma * :math.pow(t, 4)

  @doc """
  Total power radiated by a sphere: P = 4πR²σT⁴.

  In stellar physics this gives the luminosity: L = 4πR²σT_eff⁴.

  ## Examples
      iex> Thermodynamics.stefan_boltzmann_total(6.957e8, 5778.0) > 0
      true
  """

  @spec stefan_boltzmann_total(number, number, number) :: float
  def stefan_boltzmann_total(radius, t, sigma \\ @stefan_sigma) when is_positive(radius) do
    4 * :math.pi() * radius * radius * sigma * :math.pow(t, 4)
  end

  @doc """
  Effective temperature of a star from luminosity and radius:
  T_eff = (L / (4πR²σ))^(1/4)

  ## Parameters
    - luminosity: L (W)
    - radius:     R (m)

  ## Examples
      iex> Thermodynamics.effective_temperature(3.828e26, 6.957e8) |> round()
      5778
  """

  @spec effective_temperature(number, number, number) :: float
  def effective_temperature(luminosity, radius, sigma \\ @stefan_sigma)
      when is_positive(radius) do
    :math.pow(luminosity / (4 * :math.pi() * radius * radius * sigma), 0.25)
  end

  @doc "Spectral radiance (per unit wavelength) from the Planck function."

  @spec planck_wavelength(number, number, number, number, number) :: float
  def planck_wavelength(lambda, t, h \\ @planck, c \\ @speed_light, k_b \\ @boltzmann)
      when is_positive(lambda) and is_positive(t) do
    2 * h * c * c / :math.pow(lambda, 5) * (:math.exp(h * c / (lambda * k_b * t)) - 1)
  end

  @doc "Spectral radiance (per unit frequency) from the Planck function."

  @spec planck_frequency(number, number, number, number, number) :: float
  def planck_frequency(nu, t, h \\ @planck, c \\ @speed_light, k_b \\ @boltzmann) do
    2 * h * :math.pow(nu, 3) /
      (c * c * (:math.exp(h * nu / (k_b * t)) - 1))
  end

  @doc "Classical Rayleigh-Jeans long-wavelength approximation to the Planck function."

  @spec rayleigh_jeans(number, number, number, number) :: float
  def rayleigh_jeans(lambda, t, k_b \\ @boltzmann, c \\ @speed_light) do
    2 * c * k_b * t / :math.pow(lambda, 4)
  end

  # ---------------------------------------------------------------------------
  # Maxwell-Boltzmann Speed Distribution
  # ---------------------------------------------------------------------------

  @doc """
  Maxwell-Boltzmann speed probability density:
  f(v) = 4π (m/2πkT)^(3/2) v² exp(-mv²/2kT)

  ## Examples
      iex> Thermodynamics.maxwell_boltzmann(500, 4.65e-26, 300) > 0
      true
  """

  @spec maxwell_boltzmann(number, number, number, number) :: float
  def maxwell_boltzmann(v, m, t, k_b \\ @boltzmann) do
    a = m / (2 * k_b * t)

    4 * :math.pi() * :math.pow(a / :math.pi(), 1.5) *
      v * v * :math.exp(-a * v * v)
  end

  # ---------------------------------------------------------------------------
  # Van der Waals Equation of State
  # ---------------------------------------------------------------------------

  @doc """
  Van der Waals pressure: P = Nk_BT/(V - Nb) - a N²/V²

  A correction to the ideal gas law accounting for molecular volume (b) and
  intermolecular attraction (a).

  ## Parameters
    - n:           Number of particles N
    - temperature: T (K)
    - volume:      V (m³)
    - a:           Attraction parameter (N·m⁴ per particle pair)
    - b:           Excluded volume per particle (m³)
    - k_b:         Boltzmann constant

  ## Examples
      iex> Thermodynamics.van_der_waals_pressure(6.022e23, 300, 0.0224, 2.253e-49, 3.92e-29) > 0
      true
  """

  @spec van_der_waals_pressure(number, number, number, number, number, number) :: float
  def van_der_waals_pressure(n, temperature, volume, a, b, k_b \\ @boltzmann)
      when is_positive(temperature) do
    n * k_b * temperature / (volume - n * b) - a * n * n / volume * volume
  end

  # ---------------------------------------------------------------------------
  # Gravitational / Astrophysical Thermodynamics
  # ---------------------------------------------------------------------------

  @doc """
  Virial temperature of a self-gravitating system:
  T_vir = μ m_H G M / (2 k_B R)

  Estimates the characteristic temperature at which a gravitationally bound gas
  cloud is in virial equilibrium.

  ## Parameters
    - mean_mol_weight: μ (dimensionless, mean molecular weight, e.g. 0.6 for ionised H/He)
    - mass:            Total mass M (kg)
    - radius:          Characteristic radius R (m)
    - k_b:             Boltzmann constant

  ## Examples
      iex> Thermodynamics.virial_temperature(0.6, 1.989e30, 6.957e8) > 0
      true
  """

  @spec virial_temperature(number, number, number, number) :: float
  def virial_temperature(mean_mol_weight, mass, radius, k_b \\ @boltzmann)
      when is_positive(mass) do
    @gravitational_constant * mean_mol_weight * @proton_mass * mass / (2 * k_b * radius)
  end

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  defp factorial(0), do: 1
  defp factorial(n) when n > 0, do: n * factorial(n - 1)
end
