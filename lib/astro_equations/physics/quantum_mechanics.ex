defmodule AstroEquations.Physics.QuantumMechanics do
  @moduledoc """
  A module implementing fundamental quantum mechanics principles and calculations.
  """

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
    h_bar = 1.0545718e-34  # Reduced Planck constant
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
    abs(inner_product) ** 2
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
      operator.(x) * abs(wavefunction.(x)) ** 2 * dx
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
end
