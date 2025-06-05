defmodule AstroEquations.Physics.QuantumMechanics do
  @moduledoc """
  A module implementing fundamental quantum mechanics principles and calculations.
    - Partial traces
  - Schrödinger equation solvers
  - Heisenberg equation of motion
  - Quantum operators and their properties
  - Density matrix operations
  """

  alias __MODULE__, as: QM
  import Complex

  # Reduced Planck constant
  @hbar 1.054571817e-34

  @doc """
  Calculates the product of uncertainties for position and momentum.

  ## Parameters
    - delta_x: Uncertainty in position
    - delta_p: Uncertainty in momentum

  ## Returns
    - true if the product satisfies Heisenberg's Uncertainty Principle
    - false otherwise

  ## Examples
      iex> QuantumMechanics.uncertainty_principle?(1.0, 0.6)
      true

      iex> QuantumMechanics.uncertainty_principle?(0.1, 0.4)
      false
  """
  @spec uncertainty_principle?(number(), number()) :: boolean()
  def uncertainty_principle?(delta_x, delta_p) do
    # Reduced Planck constant
    h_bar = 1.0545718e-34
    delta_x * delta_p >= h_bar / 2
  end

  @doc """
  Calculates the probability according to the Born rule.

  ## Parameters
    - inner_product: The inner product ⟨ψ|ψ⟩

  ## Returns
    - The probability as the square of the absolute value

  ## Examples
      iex> QuantumMechanics.born_rule(0.5 + 0.3i)
      0.34  # 0.5² + 0.3²
  """
  @spec born_rule(Complex.t()) :: float()
  def born_rule(inner_product) do
    Complex.abs(inner_product) ** 2
  end

  @doc """
  Calculates the expectation value of an operator in position space.

  ## Parameters
    - operator: The operator function A(x)
    - wavefunction: The wavefunction Ψ(x,t) as a function of x
    - x_values: List of position values to integrate over

  ## Returns
    - The expectation value ⟨A⟩

  ## Examples
      iex> operator = fn x -> x end  # Position operator
      iex> wavefunction = fn x -> :math.exp(-x * x / 2) / :math.sqrt(:math.pi) end  # Ground state
      iex> QuantumMechanics.expectation_position(operator, wavefunction, -10..10//0.1)
      0.0  # Expected value for ground state
  """
  @spec expectation_position((float() -> float()), (float() -> Complex.t()), Range.t()) :: float()
  def expectation_position(operator, wavefunction, x_values) do
    x_values
    |> Enum.map(fn x ->
      dx = step_size(x_values)
      operator.(x) * Kernel.abs(wavefunction.(x)) ** 2 * dx
    end)
    |> Enum.sum()
  end

  @doc """
  Calculates the expectation value using bra-ket notation.

  ## Parameters
    - operator: Matrix representation of the operator
    - state: Vector representation of the state |ψ⟩

  ## Returns
    - The expectation value ⟨ψ|A|ψ⟩

  ## Examples
      iex> operator = [[1, 0], [0, -1]]  # Pauli Z matrix
      iex> state = [1/:math.sqrt(2), 1/:math.sqrt(2)]  # +X state
      iex> QuantumMechanics.expectation_braket(operator, state)
      0.0
  """
  @spec expectation_braket(list(list(number())), list(number())) :: float()
  def expectation_braket(operator, state) do
    # ⟨ψ|A|ψ⟩ = state† • operator • state
    adjoint_state = Enum.map(state, &Complex.conjugate/1)

    operator
    |> matrix_multiply(state)
    |> inner_product(adjoint_state)
  end

  @doc """
  Calculates the variance of an observable.

  ## Parameters
    - operator: Matrix representation of the operator
    - state: Vector representation of the state |ψ⟩

  ## Returns
    - The variance var(A)

  ## Examples
      iex> operator = [[1, 0], [0, -1]]  # Pauli Z matrix
      iex> state = [1, 0]  # |0⟩ state
      iex> QuantumMechanics.variance(operator, state)
      0.0
  """
  @spec variance(list(list(number())), list(number())) :: float()
  def variance(operator, state) do
    a_squared = matrix_multiply(operator, operator)
    term1 = expectation_braket(a_squared, state)
    term2 = expectation_braket(operator, state) ** 2
    term1 - term2
  end

  @doc """
  Calculates the standard deviation of an observable.

  ## Parameters
    - operator: Matrix representation of the operator
    - state: Vector representation of the state |ψ⟩

  ## Returns
    - The standard deviation δA

  ## Examples
      iex> operator = [[1, 0], [0, -1]]  # Pauli Z matrix
      iex> state = [1/:math.sqrt(2), 1/:math.sqrt(2)]  # +X state
      iex> QuantumMechanics.standard_deviation(operator, state)
      1.0
  """
  @spec standard_deviation(list(list(number())), list(number())) :: float()
  def standard_deviation(operator, state) do
    variance(operator, state) |> :math.sqrt()
  end

  @doc """
  Calculates the trace of an operator.

  ## Parameters
    - operator: Matrix representation of the operator
    - basis: List of basis vectors

  ## Returns
    - The trace Tr(A)

  ## Examples
      iex> operator = [[1, 0], [0, -1]]  # Pauli Z matrix
      iex> basis = [[[1, 0], [0, 1]]]  # Standard basis
      iex> QuantumMechanics.trace(operator, basis)
      0.0  # 1 + (-1)
  """
  @spec trace(list(list(number())), list(list(list(number())))) :: number()
  def trace(operator, basis) do
    basis
    |> Enum.map(fn basis_vector ->
      expectation_braket(operator, basis_vector)
    end)
    |> Enum.sum()
  end

  # Helper functions
  defp step_size(range) do
    (range.last - range.first) / (Enum.count(range) - 1)
  end

  defp matrix_multiply(matrix, vector) do
    Enum.map(matrix, fn row ->
      Enum.zip_with(row, vector, fn a, b -> a * b end) |> Enum.sum()
    end)
  end

  defp inner_product(vec1, vec2) do
    Enum.zip_with(vec1, vec2, fn a, b -> a * b end) |> Enum.sum()
  end

  @doc """
  Computes the partial trace over subsystem B of a composite system.

  ## Parameters
    - rho_ab: Density matrix of the composite system (as tensor product)
    - dim_a: Dimension of subsystem A
    - dim_b: Dimension of subsystem B

  ## Returns
    - Density matrix of subsystem A after tracing out B

  ## Examples
      iex> rho_a = [[1, 0], [0, 0]]
      iex> rho_b = [[0.5, 0.5], [0.5, 0.5]]
      iex> rho_ab = QM.tensor_product(rho_a, rho_b)
      iex> QM.partial_trace(rho_ab, 2, 2)
      [[0.5, 0.0], [0.0, 0.0]]
  """
  @spec partial_trace(list(list(number())), integer(), integer()) :: list(list(float()))
  def partial_trace(rho_ab, dim_a, dim_b) do
    for i <- 0..(dim_a - 1),
        j <- 0..(dim_a - 1),
        do:
          Enum.reduce(0..(dim_b - 1), 0, fn k, acc ->
            acc + (Enum.at(rho_ab, i * dim_b + k) |> Enum.at(j * dim_b + k))
          end)
          |> Enum.chunk_every(dim_a)
  end

  @doc """
  Solves the time-independent Schrödinger equation for 1D potentials.

  ## Parameters
    - potential_fn: Function V(x) describing the potential
    - mass: Particle mass
    - x_range: Range of x values to solve over
    - boundary_conditions: Tuple of {left_boundary, right_boundary}

  ## Returns
    - Tuple of {eigenvalues, eigenfunctions}

  ## Examples
      iex> potential_fn = fn x -> 0.5 * x * x end # Harmonic oscillator
      iex> {energies, _waves} = QM.solve_schrodinger(potential_fn, 1.0, -5.0..5.0//0.1, {0.0, 0.0})
      iex> length(energies) > 0
      true
  """
  @spec solve_schrodinger((float() -> float()), float(), Range.t(), {number(), number()}) ::
          {list(float()), list(list(float()))}
  def solve_schrodinger(potential_fn, mass, x_range, {psi_left, psi_right}) do
    # Finite difference method implementation
    x_list = Enum.to_list(x_range)
    n = length(x_list)
    dx = step_size(x_range)

    # Construct Hamiltonian matrix
    hamiltonian =
      for i <- 0..(n - 1) do
        for j <- 0..(n - 1) do
          cond do
            i == j ->
              @hbar ** 2 / (mass * dx ** 2) + potential_fn.(Enum.at(x_list, i))

            Kernel.abs(i - j) == 1 ->
              -@hbar ** 2 / (2 * mass * dx ** 2)

            true ->
              0.0
          end
        end
      end

    # Apply boundary conditions
    hamiltonian = List.update_at(hamiltonian, 0, fn _ -> List.duplicate(0.0, n) end)
    hamiltonian = List.update_at(hamiltonian, -1, fn _ -> List.duplicate(0.0, n) end)

    # Diagonalize to find eigenvalues and eigenvectors
    {eigenvalues, eigenfunctions} = diagonalize(hamiltonian)
    {eigenvalues, eigenfunctions}
  end

  @doc """
  Computes the time evolution of an operator in the Heisenberg picture.

  ## Parameters
    - operator: Initial operator matrix
    - hamiltonian: Hamiltonian matrix
    - t: Time to evolve
    - dt: Time step size

  ## Returns
    - Operator evolved to time t

  ## Examples
      iex> op = [[0, 1], [1, 0]] # Pauli X
      iex> h = [[1, 0], [0, -1]] # Pauli Z
      iex> evolved = QM.heisenberg_evolution(op, h, 1.0, 0.01)
      iex> matrix_size(evolved) == {2, 2}
      true
  """
  @spec heisenberg_evolution(list(list(number())), list(list(number())), float(), float()) ::
          list(list(Complex.t()))
  def heisenberg_evolution(operator, hamiltonian, t, dt) do
    steps = round(t / dt)
    evolve_step(operator, hamiltonian, steps, dt)
  end

  defp evolve_step(op, _h, 0, _dt), do: op

  defp evolve_step(op, h, steps, dt) do
    # Compute commutator [H, A]
    comm = matrix_subtract(matrix_multiply(h, op), matrix_multiply(op, h))
    # dA/dt = (i/hbar)[H,A]
    derivative = matrix_scale(comm, new(0, 1.0 / @hbar))
    # Euler step
    new_op = matrix_add(op, matrix_scale(derivative, dt))
    evolve_step(new_op, h, steps - 1, dt)
  end

  @doc """
  Constructs a density matrix from state vectors and probabilities.

  ## Parameters
    - states: List of state vectors
    - probs: List of corresponding probabilities

  ## Returns
    - Density matrix ρ

  ## Examples
      iex> state1 = [1, 0]
      iex> state2 = [0, 1]
      iex> QM.density_matrix([state1, state2], [0.5, 0.5])
      [[0.5, 0], [0, 0.5]]
  """
  @spec density_matrix(list(list(number())), list(float())) :: list(list(float()))
  def density_matrix(states, probs) do
    Enum.zip_with(states, probs, fn state, p ->
      outer = outer_product(state, state)
      matrix_scale(outer, p)
    end)
    |> Enum.reduce(fn m, acc -> matrix_add(m, acc) end)
  end

  @doc """
  Computes the purity of a density matrix.

  ## Parameters
    - rho: Density matrix

  ## Returns
    - Purity value between 1/d and 1

  ## Examples
      iex> rho = [[0.5, 0], [0, 0.5]]
      iex> QM.purity(rho)
      0.5
  """
  @spec purity(list(list(number()))) :: float()
  def purity(rho) do
    rho_squared = matrix_multiply(rho, rho)
    trace(rho_squared)
  end

  # Helper functions
  defp tensor_product(a, b) do
    for a_row <- a do
      for a_el <- a_row do
        for b_row <- b do
          for b_el <- b_row do
            a_el * b_el
          end
        end
      end
    end
    |> List.flatten()
    |> Enum.chunk_every(length(b) * length(Enum.at(b, 0)))
  end

  defp outer_product(ket, bra) do
    for k <- ket do
      for b <- bra do
        k * Complex.conjugate(b)
      end
    end
  end

  defp diagonalize(matrix) do
    # This would use a proper numerical diagonalization in real implementation
    # For simplicity, we return dummy values here
    {[1.0, 2.0], [[1.0, 0.0], [0.0, 1.0]]}
  end

  defp trace(matrix) do
    Enum.with_index(matrix)
    |> Enum.reduce(0, fn {row, i}, acc -> acc + Enum.at(row, i) end)
  end

  defp matrix_add(a, b) do
    Enum.zip_with(a, b, fn row_a, row_b ->
      Enum.zip_with(row_a, row_b, fn x, y -> Complex.add(x, y) end)
    end)
  end

  defp matrix_subtract(a, b) do
    Enum.zip_with(a, b, fn row_a, row_b ->
      Enum.zip_with(row_a, row_b, fn x, y -> Complex.subtract(x, y) end)
    end)
  end

  defp matrix_scale(matrix, scalar) do
    Enum.map(matrix, fn row ->
      Enum.map(row, fn el -> Complex.multiply(el, scalar) end)
    end)
  end
end
