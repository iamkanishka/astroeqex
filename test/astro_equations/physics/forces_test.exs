defmodule AstroEquations.Physics.ForcesTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.Forces

  describe "buoyancy/2" do
    test "calculates buoyant force" do
      assert Physics.buoyancy(5) == 49.05
      assert Physics.buoyancy(10, 9.8) == 98.0
    end
  end

  describe "buoyancy_from_density/3" do
    test "calculates buoyant force from density and volume" do
      assert Physics.buoyancy_from_density(1000, 0.005) == 49.05
      assert Physics.buoyancy_from_density(800, 0.01, 10) == 80.0
    end
  end

  describe "friction calculations" do
    test "kinetic friction" do
      assert Physics.kinetic_friction(0.3, 10) == 3.0
    end

    test "static friction" do
      assert Physics.static_friction(0.4, 10) == 4.0
    end
  end
end
