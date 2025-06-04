defmodule AstroEquations.Physics.SpecialRelativityTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.SpecialRelativity

  @c 299_792_458

  describe "gamma_factor/2" do
    test "calculates gamma factor for v=0" do
      assert AstroEquations.Physics.SpecialRelativity.gamma_factor(0) == 1.0
    end

    test "calculates gamma factor for v=0.5c" do
      v = 0.5 * @c
      assert_in_delta AstroEquations.Physics.SpecialRelativity.gamma_factor(v), 1.1547, 0.0001
    end
  end

  describe "time_dilation/3" do
    test "proper time equals dilated time at v=0" do
      assert AstroEquations.Physics.SpecialRelativity.time_dilation(1, 0) == 1.0
    end

    test "time dilation at v=0.5c" do
      v = 0.5 * @c
      assert_in_delta AstroEquations.Physics.SpecialRelativity.time_dilation(1, v), 1.1547, 0.0001
    end
  end

  describe "length_contraction/3" do
    test "proper length equals contracted length at v=0" do
      assert AstroEquations.Physics.SpecialRelativity.length_contraction(1, 0) == 1.0
    end

    test "length contraction at v=0.5c" do
      v = 0.5 * @c
      assert_in_delta AstroEquations.Physics.SpecialRelativity.length_contraction(1, v), 0.8660, 0.0001
    end
  end

  describe "rest_energy/2" do
    test "rest energy calculation" do
      assert AstroEquations.Physics.SpecialRelativity.rest_energy(1) == 1 * @c ** 2
    end
  end

  describe "relative_velocity/3" do
    test "relative velocity when u=v" do
      assert AstroEquations.Physics.SpecialRelativity.relative_velocity(100_000_000, 100_000_000) == 0.0
    end

    test "relative velocity when u=0.8c and v=0.5c" do
      u = 0.8 * @c
      v = 0.5 * @c
      expected = (u - v) / (1 - (u * v) / @c ** 2)
      assert AstroEquations.Physics.SpecialRelativity.relative_velocity(u, v) == expected
    end
  end

  # Additional tests for other functions would follow the same pattern
end
