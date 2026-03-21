defmodule AstroEquations.Physics.QuantumMechanics do
  @moduledoc """
  Fundamental quantum mechanics.

  Covers:
  - Heisenberg uncertainty principle
  - de Broglie wavelength
  - Energy levels (hydrogen atom, infinite square well, harmonic oscillator)
  - Photon energy and momentum
  - Born rule and probability
  - Expectation values (position space and matrix form)
  - Pauli matrices and spin operators
  - Ladder operators (harmonic oscillator, atomic two-level)
  - Density matrix and purity
  - Commutators and Hilbert-Schmidt norm
  - Basic scattering (transmission/reflection coefficients)
  - Compton wavelength and Bohr magneton

  Note: Matrix operations use plain Elixir lists (nested lists for matrices).
  For production use, consider a linear-algebra dependency.
  """

  # Reduced Planck constant (J·s)
  @hbar 1.054571817e-34
  # Planck constant (J·s)
  @h 6.62607015e-34
  # Electron mass (kg)
  @m_electron 9.10938370e-31
  # Speed of light (m/s)
  @speed_of_light 2.99792458e8
  # Bohr radius (m)
  @a0 5.29177210903e-11
  # Fine-structure constant
  @alpha 7.2973525693e-3
  # Boltzmann constant (J/K)
  @boltzmann 1.380649e-23
  # Elementary charge (C)
  @elementary_charge 1.602176634e-19

  # ---------------------------------------------------------------------------
  # Fundamental Relations
  # ---------------------------------------------------------------------------

  @doc """
  Heisenberg uncertainty principle check: Δx Δp ≥ ħ/2

  ## Returns
    true if the product satisfies the principle

  ## Examples
      iex> QuantumMechanics.uncertainty_principle?(1.0, 0.6)
      true
  """
  @spec uncertainty_principle?(number, number) :: boolean
  def uncertainty_principle?(delta_x, delta_p) do
    delta_x * delta_p >= @hbar / 2
  end

  @doc "Minimum position–momentum uncertainty product: Δx Δp_min = ħ/2."
  @spec min_uncertainty_product() :: float
  def min_uncertainty_product, do: @hbar / 2

  @doc "Checks the energy-time uncertainty relation: ΔE Δt ≥ ħ/2."
  @spec energy_time_uncertainty?(number, number) :: boolean
  def energy_time_uncertainty?(delta_e, delta_t) do
    delta_e * delta_t >= @hbar / 2
  end

  @doc """
  de Broglie wavelength: λ = h / p

  ## Parameters
    - momentum: Particle momentum (kg·m/s)

  ## Examples
      iex> QuantumMechanics.de_broglie_wavelength(1.0e-24) > 0
      true
  """
  @spec de_broglie_wavelength(number) :: float
  def de_broglie_wavelength(momentum), do: @h / momentum

  @doc """
  de Broglie wavelength from kinetic energy: λ = h / √(2 m KE)

  ## Parameters
    - kinetic_energy: KE (J)
    - mass:           Particle mass (kg, default: electron mass)

  ## Examples
      iex> QuantumMechanics.de_broglie_wavelength_ke(1.602e-19) > 0
      true
  """
  @spec de_broglie_wavelength_ke(number, number) :: float
  def de_broglie_wavelength_ke(kinetic_energy, mass \\ @m_electron) do
    @h / :math.sqrt(2 * mass * kinetic_energy)
  end

  @doc """
  Photon energy: E = h f = ħ omega

  ## Parameters
    - frequency: f (Hz)

  ## Examples
      iex> QuantumMechanics.photon_energy(6.0e14) > 0
      true
  """
  @spec photon_energy(number) :: float
  def photon_energy(frequency), do: @h * frequency

  @doc "Photon momentum from de Broglie: p = h/λ."
  @spec photon_momentum(number) :: float
  def photon_momentum(wavelength), do: @h / wavelength

  @doc """
  Compton wavelength: λ_C = h / (m c)

  The length scale at which relativistic and quantum effects become comparable.

  ## Parameters
    - mass: Particle mass (kg, default: electron mass)

  ## Returns
    Compton wavelength (m)

  ## Examples
      iex> QuantumMechanics.compton_wavelength() > 0
      true
  """
  @spec compton_wavelength(number) :: float
  def compton_wavelength(mass \\ @m_electron) do
    @h / (mass * @speed_of_light)
  end

  @doc """
  Bohr magneton: μ_B = eħ / (2 m_e)

  The natural unit of electron magnetic dipole moment.

  ## Returns
    Bohr magneton in J/T

  ## Examples
      iex> QuantumMechanics.bohr_magneton() > 0
      true
  """
  @spec bohr_magneton() :: float
  def bohr_magneton do
    @elementary_charge * @hbar / (2 * @m_electron)
  end

  # ---------------------------------------------------------------------------
  # Hydrogen Atom Energy Levels
  # ---------------------------------------------------------------------------

  @doc """
  Hydrogen atom energy levels: E_n = -13.6 eV / n²

  ## Parameters
    - n: Principal quantum number (positive integer)

  ## Returns
    Energy in electron-volts

  ## Examples
      iex> QuantumMechanics.hydrogen_energy_level(1) |> Float.round(2)
      -13.6
  """
  @spec hydrogen_energy_level(number) :: float
  def hydrogen_energy_level(n), do: -13.6 / (n * n)

  @doc """
  Hydrogen atom energy levels in joules.

  E_n = -m_e e⁴ / (8 ε₀² h² n²)

  ## Examples
      iex> QuantumMechanics.hydrogen_energy_joules(1) < 0
      true
  """
  @spec hydrogen_energy_joules(number) :: float
  def hydrogen_energy_joules(n) do
    e = 1.602176634e-19
    eps0 = 8.8541878128e-12
    -@m_electron * e ** 4 / (8 * eps0 ** 2 * @h ** 2 * n ** 2)
  end

  @doc "Returns the Bohr radius: a₀ = 4πε₀ħ²/(mₑe²) ≈ 5.292×10⁻¹¹ m."
  @spec bohr_radius() :: float
  def bohr_radius, do: @a0

  @doc """
  Rydberg formula for spectral lines: 1/λ = R_∞ (1/n₁² - 1/n₂²)

  ## Parameters
    - n1:    Lower principal quantum number
    - n2:    Upper principal quantum number
    - r_inf: Rydberg constant (m⁻¹, default: 1.0973731568×10⁷)

  ## Returns
    Wavelength in metres
  """
  @spec rydberg_wavelength(number, number, number) :: float
  def rydberg_wavelength(n1, n2, r_inf \\ 1.0973731568e7) when n2 > n1 do
    1.0 / (r_inf * (1 / (n1 * n1) - 1 / (n2 * n2)))
  end

  # ---------------------------------------------------------------------------
  # Energy Levels — Idealised Systems
  # ---------------------------------------------------------------------------

  @doc """
  Particle-in-a-box (infinite square well) energy levels:
  E_n = n² π² ħ² / (2 m L²)

  ## Parameters
    - n:    Quantum number
    - mass: Particle mass (kg)
    - l:    Box length (m)

  ## Examples
      iex> QuantumMechanics.infinite_well_energy(1, 9.109e-31, 1.0e-9) > 0
      true
  """
  @spec infinite_well_energy(number, number, number) :: float
  def infinite_well_energy(n, mass, l) do
    n * n * :math.pi() ** 2 * @hbar ** 2 / (2 * mass * l * l)
  end

  @doc """
  Quantum harmonic oscillator energy levels: E_n = ħ omega (n + ½)

  ## Parameters
    - n:     Quantum number (0, 1, 2, …)
    - omega: Angular frequency (rad/s)

  ## Examples
      iex> QuantumMechanics.harmonic_oscillator_energy(0, 1.0e14) > 0
      true
  """
  @spec harmonic_oscillator_energy(number, number) :: float
  def harmonic_oscillator_energy(n, omega), do: @hbar * omega * (n + 0.5)

  # ---------------------------------------------------------------------------
  # Born Rule & Probability
  # ---------------------------------------------------------------------------

  @doc """
  Born rule probability: P = |⟨ψ|φ⟩|²

  ## Examples
      iex> QuantumMechanics.born_rule(1.0) |> Float.round(4)
      1.0
  """
  @spec born_rule(number) :: float
  def born_rule(inner_product_magnitude), do: inner_product_magnitude ** 2

  # ---------------------------------------------------------------------------
  # Expectation Values (Matrix Form)
  # ---------------------------------------------------------------------------

  @doc """
  Expectation value ⟨A⟩ = ⟨ψ|A|ψ⟩ for real state vector and Hermitian operator.

  ## Examples
      iex> QuantumMechanics.expectation_braket([[1,0],[0,-1]], [1.0, 0.0]) |> Float.round(4)
      1.0
  """
  @spec expectation_braket([[number]], [number]) :: float
  def expectation_braket(operator, state) do
    a_psi = matrix_multiply(operator, state)
    inner_product(state, a_psi)
  end

  @doc "Quantum mechanical variance of an observable: Var(A) = ⟨A²⟩ − ⟨A⟩²."
  @spec variance([[number]], [number]) :: float
  def variance(operator, state) do
    a2 = matrix_multiply(operator, operator)
    expectation_braket(a2, state) - expectation_braket(operator, state) ** 2
  end

  @doc "Quantum mechanical standard deviation: δA = √Var(A)."
  @spec standard_deviation([[number]], [number]) :: float
  def standard_deviation(operator, state) do
    :math.sqrt(abs(variance(operator, state)))
  end

  @doc "Expectation value ⟨x⟩ via numerical integration in position space."
  @spec expectation_position((number -> number), (number -> number), Enumerable.t()) :: float
  def expectation_position(operator, wavefunction, x_values) do
    n = Enum.count(x_values) - 1
    x_list = Enum.to_list(x_values)
    dx = if n > 0, do: (List.last(x_list) - List.first(x_list)) / n, else: 1.0

    Enum.reduce(x_list, 0.0, fn x, acc ->
      acc + operator.(x) * wavefunction.(x) ** 2 * dx
    end)
  end

  # ---------------------------------------------------------------------------
  # Pauli Matrices
  # ---------------------------------------------------------------------------

  @doc "Pauli X matrix (bit-flip operator): σₓ = [[0,1],[1,0]]."
  @spec pauli_x() :: [[number]]
  def pauli_x, do: [[0, 1], [1, 0]]

  @doc "Pauli Z matrix (phase-flip operator): σ_z = [[1,0],[0,−1]]."
  @spec pauli_z() :: [[number]]
  def pauli_z, do: [[1, 0], [0, -1]]

  @doc "2×2 identity matrix."
  @spec identity() :: [[number]]
  def identity, do: [[1, 0], [0, 1]]

  @doc "Atomic two-level raising operator: σ₊ = |e⟩⟨g|."
  @spec atomic_raise() :: [[number]]
  def atomic_raise, do: [[0, 1], [0, 0]]

  @doc "Atomic two-level lowering operator: σ₋ = |g⟩⟨e|."
  @spec atomic_lower() :: [[number]]
  def atomic_lower, do: [[0, 0], [1, 0]]

  @doc "Atomic population difference operator: σ_z = |e⟩⟨e| − |g⟩⟨g|."
  @spec atomic_sigma_z() :: [[number]]
  def atomic_sigma_z, do: [[1, 0], [0, -1]]

  @doc "Applies the atomic raising operator σ₊ to a two-level state vector."
  @spec apply_raise([number]) :: [number]
  def apply_raise([g, _e]), do: [0, g]

  @doc "Applies the atomic lowering operator σ₋ to a two-level state vector."
  @spec apply_lower([number]) :: [number]
  def apply_lower([_g, e]), do: [e, 0]

  # ---------------------------------------------------------------------------
  # Ladder Operators (Harmonic Oscillator)
  # ---------------------------------------------------------------------------

  @doc "Applies the bosonic annihilation (lowering) operator: â|n⟩ = √n |n-1⟩."
  @spec annihilate([number], non_neg_integer) :: [float]
  def annihilate(state, n) do
    Enum.with_index(state, fn val, i ->
      if i == n - 1, do: val * :math.sqrt(n), else: 0.0
    end)
  end

  @doc "Applies the bosonic creation (raising) operator: â†|n⟩ = √(n+1) |n+1⟩."
  @spec create([number], non_neg_integer) :: [float]
  def create(state, n) do
    Enum.with_index(state, fn val, i ->
      if i == n + 1, do: val * :math.sqrt(n + 1), else: 0.0
    end)
  end

  # ---------------------------------------------------------------------------
  # Density Matrix
  # ---------------------------------------------------------------------------

  @doc """
  Builds a density matrix ρ = Σᵢ pᵢ |ψᵢ⟩⟨ψᵢ|

  ## Examples
      iex> QuantumMechanics.density_matrix([[1, 0]], [1.0]) |> hd() |> hd() |> Float.round(4)
      1.0
  """
  @spec density_matrix([[number]], [float]) :: [[float]]
  def density_matrix(states, probs) do
    Enum.zip_with(states, probs, fn state, p ->
      matrix_scale(outer_product(state, state), p)
    end)
    |> Enum.reduce(fn m, acc -> matrix_add(m, acc) end)
  end

  @doc "Purity of a quantum state: Tr(ρ²), equal to 1 for a pure state."
  @spec purity([[number]]) :: float
  def purity(rho) do
    rho2 = matrix_multiply(rho, rho)
    matrix_trace(rho2)
  end

  @doc "Trace of a square matrix: Tr(A) = Σᵢ Aᵢᵢ."
  @spec matrix_trace([[number]]) :: float
  def matrix_trace(matrix) do
    Enum.with_index(matrix)
    |> Enum.reduce(0, fn {row, i}, acc -> acc + Enum.at(row, i, 0) end)
  end

  @doc "Hilbert-Schmidt (Frobenius) norm of a matrix: ||A||_HS = √(Tr(A†A))."
  @spec hilbert_schmidt_norm([[number]]) :: float
  def hilbert_schmidt_norm(matrix) do
    matrix
    |> Enum.flat_map(& &1)
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
    |> :math.sqrt()
  end

  # ---------------------------------------------------------------------------
  # Scattering (Step Potential)
  # ---------------------------------------------------------------------------

  @doc """
  Transmission coefficient through a rectangular barrier (E > V₀): T = 4k₁k₂/(k₁+k₂)²

  ## Examples
      iex> QuantumMechanics.transmission_coefficient(5.0e9, 5.0e9) |> Float.round(4)
      1.0
  """
  @spec transmission_coefficient(number, number) :: float
  def transmission_coefficient(k1, k2) do
    4 * k1 * k2 / :math.pow(k1 + k2, 2)
  end

  @doc "Quantum mechanical reflection coefficient at a potential step: R = ((k₁−k₂)/(k₁+k₂))²."
  @spec reflection_coefficient(number, number) :: float
  def reflection_coefficient(k1, k2) do
    :math.pow((k1 - k2) / (k1 + k2), 2)
  end

  @doc """
  Wave vector from energy and potential: k = √(2m(E-V)) / ħ

  ## Examples
      iex> QuantumMechanics.wave_vector(9.109e-31, 5.0e-19, 0.0) > 0
      true
  """
  @spec wave_vector(number, number, number) :: float
  def wave_vector(mass, energy, potential) do
    :math.sqrt(max(2 * mass * (energy - potential), 0.0)) / @hbar
  end

  @doc """
  Returns the **fine-structure constant (α)**.

  The fine-structure constant is a fundamental dimensionless
  constant that characterizes the strength of the electromagnetic interaction.

  Value:

      α ≈ 7.2973525693 × 10⁻³

  ## Examples

      iex> AstroEquations.Physics.QuantumMechanics.fine_structure_constant()
      0.0072973525693
  """
  @spec fine_structure_constant() :: float
  def fine_structure_constant do
    @alpha
  end

  @doc """
  Calculates **thermal energy** using the Boltzmann relation.

  Formula:

      E = k_B * T

  where:

  - `E` = thermal energy (Joules)
  - `k_B` = Boltzmann constant
  - `T` = temperature in Kelvin

  ## Parameters

  - `temperature` — temperature in **Kelvin**

  ## Returns

  - Thermal energy in **Joules**

  ## Examples

      iex> AstroEquations.Physics.QuantumMechanics.thermal_energy(300)
      4.141947e-21
  """
  @spec thermal_energy(number) :: float
  def thermal_energy(temperature) do
    @boltzmann * temperature
  end

  # ---------------------------------------------------------------------------
  # Private helpers
  # ---------------------------------------------------------------------------

  defp matrix_multiply(a, b) when is_list(b) and b != [] and not is_list(hd(b)) do
    Enum.map(a, fn row ->
      Enum.zip_with(row, b, fn x, y -> x * y end) |> Enum.sum()
    end)
  end

  defp matrix_multiply(a, b) do
    n = length(b)
    b_t = Enum.map(0..(n - 1), fn j -> Enum.map(b, &Enum.at(&1, j)) end)

    Enum.map(a, fn row ->
      Enum.map(b_t, fn col ->
        Enum.zip_with(row, col, fn x, y -> x * y end) |> Enum.sum()
      end)
    end)
  end

  defp inner_product(v1, v2) do
    Enum.zip_with(v1, v2, fn a, b -> a * b end) |> Enum.sum()
  end

  defp outer_product(ket, bra) do
    Enum.map(ket, fn k -> Enum.map(bra, fn b -> k * b end) end)
  end

  defp matrix_add(a, b) do
    Enum.zip_with(a, b, fn ra, rb ->
      Enum.zip_with(ra, rb, fn x, y -> x + y end)
    end)
  end

  defp matrix_scale(matrix, s) do
    Enum.map(matrix, fn row -> Enum.map(row, fn x -> x * s end) end)
  end
end
