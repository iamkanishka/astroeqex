defmodule AstroEquations.Physics.QuantumMechanicsTest do
  use ExUnit.Case
  import Complex
  doctest AstroEquations.Physics.QuantumMechanics

  test "uncertainty principle" do
    assert AstroEquations.Physics.QuantumMechanics.uncertainty_principle?(1.0, 1.0545718e-34 / 2) ==
             true

    assert AstroEquations.Physics.QuantumMechanics.uncertainty_principle?(1.0, 1.0545718e-34 / 3) ==
             false
  end

  test "born rule" do
    assert_in_delta AstroEquations.Physics.QuantumMechanics.born_rule(new(0.5, 0.3)),
                    0.34,
                    1.0e-10

    assert_in_delta AstroEquations.Physics.QuantumMechanics.born_rule(new(1, 0)), 1.0, 1.0e-10
  end

  test "expectation position" do
    operator = fn x -> x end
    wavefunction = fn x -> new(:math.exp(-x * x / 2) / :math.sqrt(:math.pi())) end

    assert_in_delta AstroEquations.Physics.QuantumMechanics.expectation_position(
                      operator,
                      wavefunction,
                      -10..10//0.1
                    ),
                    0.0,
                    1.0e-10
  end

  test "expectation braket" do
    operator = [[1, 0], [0, -1]]
    state = [1 / :math.sqrt(2), 1 / :math.sqrt(2)]

    assert_in_delta AstroEquations.Physics.QuantumMechanics.expectation_braket(operator, state),
                    0.0,
                    1.0e-10
  end

  test "variance" do
    operator = [[1, 0], [0, -1]]
    state = [1, 0]

    assert_in_delta AstroEquations.Physics.QuantumMechanics.variance(operator, state),
                    0.0,
                    1.0e-10
  end

  test "standard deviation" do
    operator = [[1, 0], [0, -1]]
    state = [1 / :math.sqrt(2), 1 / :math.sqrt(2)]

    assert_in_delta AstroEquations.Physics.QuantumMechanics.standard_deviation(operator, state),
                    1.0,
                    1.0e-10
  end

  test "trace" do
    operator = [[1, 0], [0, -1]]
    basis = [[1, 0], [0, 1]]
    assert_in_delta AstroEquations.Physics.QuantumMechanics.trace(operator, [basis]), 0.0, 1.0e-10
  end


    test "partial trace" do
    rho_a = [[1, 0], [0, 0]]
    rho_b = [[0.5, 0.5], [0.5, 0.5]]
    rho_ab = AstroEquations.Physics.QuantumMechanics.tensor_product(rho_a, rho_b)
    assert AstroEquations.Physics.QuantumMechanics.partial_trace(rho_ab, 2, 2) == [[0.5, 0.0], [0.0, 0.0]]
  end

  test "density matrix construction" do
    state1 = [1, 0]
    state2 = [0, 1]
    assert AstroEquations.Physics.QuantumMechanics.density_matrix([state1, state2], [0.5, 0.5]) == [[0.5, 0], [0, 0.5]]
  end

  test "purity calculation" do
    pure_state = [[1, 0], [0, 0]]
    assert AstroEquations.Physics.QuantumMechanics.purity(pure_state) == 1.0

    mixed_state = [[0.5, 0], [0, 0.5]]
    assert AstroEquations.Physics.QuantumMechanics.purity(mixed_state) == 0.5
  end

  test "heisenberg evolution" do
    # Pauli X
    op = [[0, 1], [1, 0]]
    # Pauli Z
    h = [[1, 0], [0, -1]]
    evolved = AstroEquations.Physics.QuantumMechanics.heisenberg_evolution(op, h, 1.0, 0.01)
    assert length(evolved) == 2
    assert length(hd(evolved)) == 2
  end

  test "schrodinger solver" do
    potential_fn = fn x -> 0.5 * x * x end
    x_range = Enum.to_list(Stream.iterate(-5.0, &(&1 + 0.1)) |> Enum.take_while(&(&1 <= 5.0)))
    {energies, _waves} = AstroEquations.Physics.QuantumMechanics.solve_schrodinger(potential_fn, 1.0, x_range, {0.0, 0.0})
    assert is_list(energies)
    assert length(energies) > 0
  end

end
