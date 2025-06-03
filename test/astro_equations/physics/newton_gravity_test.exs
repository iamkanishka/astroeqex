
defmodule AstroEquations.Physics.NewtonGravityTest do
  use ExUnit.Case
  doctest  AstroEquations.Physics.NewtonGravity

  describe "force/3" do
    test "calculates gravitational force" do
      assert_in_delta  AstroEquations.Physics.NewtonGravity.force(5.972e24, 7.348e22, 3.844e8),
                      1.98e20, 1.0e18
    end

    test "raises error for zero distance" do
      assert_raise ArgumentError, fn ->
         AstroEquations.Physics.NewtonGravity.force(1, 1, 0)
      end
    end
  end

  describe "potential/2" do
    test "calculates gravitational potential" do
      assert_in_delta  AstroEquations.Physics.NewtonGravity.potential(5.972e24, 6.371e6),
                      -6.25e7, 1.0e6
    end
  end

  describe "field/2" do
    test "calculates gravitational field" do
      assert_in_delta  AstroEquations.Physics.NewtonGravity.field(5.972e24, 6.371e6),
                      9.81, 0.1
    end
  end

  describe "potential_energy/3" do
    test "calculates gravitational potential energy" do
      assert_in_delta  AstroEquations.Physics.NewtonGravity.potential_energy(1000, 5.972e24, 6.371e6),
                      -6.25e10, 1.0e9
    end
  end

  describe "approximate_potential_energy/3" do
    test "calculates mgh approximation" do
      assert  AstroEquations.Physics.NewtonGravity.approximate_potential_energy(1000, 9.81, 100) == 981000.0
    end
  end

  describe "orbital_period/3" do
    test "calculates orbital period (Earth-Moon)" do
      # ~27.3 days in seconds
      assert_in_delta  AstroEquations.Physics.NewtonGravity.orbital_period(3.844e8, 5.972e24, 7.348e22),
                      2.36e6, 1.0e5
    end
  end
end
