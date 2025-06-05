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
end
