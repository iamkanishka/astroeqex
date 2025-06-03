defmodule AstroEquations.Physics.MaterialsTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.Materials

  describe "density/2" do
    test "calculates density correctly" do
      assert Physics.Materials.density(10, 2) == 5.0
      assert Physics.Materials.density(5.5, 2.2) == 2.5
    end

    test "raises error when volume is zero" do
      assert_raise ArgumentError, fn ->
        Physics.Materials.density(10, 0)
      end
    end
  end

  describe "continuous_density/3" do
    test "calculates continuous density for linear function" do
      linear_mass = fn v -> 2.5 * v end
      assert Physics.Materials.continuous_density(linear_mass, 2.0) == 2.5
      assert Physics.Materials.continuous_density(linear_mass, 5.0) == 2.5
    end

    test "calculates continuous density for quadratic function" do
      quadratic_mass = fn v -> 0.5 * :math.pow(v, 2) end
      assert_in_delta Physics.Materials.continuous_density(quadratic_mass, 2.0), 2.0, 1.0e-6
      assert_in_delta Physics.Materials.continuous_density(quadratic_mass, 4.0), 4.0, 1.0e-6
    end
  end
end
