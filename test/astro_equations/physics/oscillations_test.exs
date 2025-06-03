defmodule AstroEquations.Physics.OscillationsTest do
  @moduledoc """

  # This module provides unit tests for the AstroEquations.Physics.Oscillations module.
  # It tests the force, potential energy, and angular frequency calculations for springs.
  # The tests cover various scenarios to ensure the correctness of the calculations.
  # The module uses ExUnit for testing and includes doctests for documentation.
  """

  use ExUnit.Case
  doctest AstroEquations.Physics.Oscillations

  describe "force/2" do
    test "calculates spring force correctly" do
      assert AstroEquations.Physics.Oscillations.force(10, 0.5) == -5.0
      assert AstroEquations.Physics.Oscillations.force(5, -0.2) == 1.0
    end
  end

  describe "potential_energy/2" do
    test "calculates potential energy correctly" do
      assert AstroEquations.Physics.Oscillations.potential_energy(10, 0.5) == 1.25
      assert AstroEquations.Physics.Oscillations.potential_energy(20, 0.1) == 0.1
    end
  end

  describe "angular_frequency/2" do
    test "calculates angular frequency correctly" do
      assert AstroEquations.Physics.Oscillations.angular_frequency(10, 2.5) == 2.0
      assert AstroEquations.Physics.Oscillations.angular_frequency(9, 1) == 3.0
    end
  end
end
