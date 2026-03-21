defmodule AstroEquations.AstrophysicsAndAstronomy.BlackHoleTest do
  use ExUnit.Case, async: true

  alias AstroEquations.AstrophysicsAndAstronomy.BlackHole

  @solar_mass 1.989e30
  @g 6.67430e-11
  @c 2.99792458e8

  describe "schwarzschild_radius/1" do
    test "1 solar mass → ≈ 2953 m" do
      assert_in_delta BlackHole.schwarzschild_radius(@solar_mass), 2953.25, 1.0
    end

    test "scales linearly with mass" do
      r1 = BlackHole.schwarzschild_radius(@solar_mass)
      r10 = BlackHole.schwarzschild_radius(10 * @solar_mass)
      assert_in_delta r10 / r1, 10.0, 1.0e-10
    end

    test "zero mass gives zero radius" do
      assert BlackHole.schwarzschild_radius(0) == 0.0
    end

    test "formula: 2GM/c²" do
      mass = @solar_mass
      expected = 2 * @g * mass / @c ** 2
      assert_in_delta BlackHole.schwarzschild_radius(mass), expected, 1.0e-6
    end
  end

  describe "schwarzschild_radius_solar/1" do
    test "1 solar mass matches schwarzschild_radius(@solar_mass)" do
      assert_in_delta BlackHole.schwarzschild_radius_solar(1),
                      BlackHole.schwarzschild_radius(@solar_mass),
                      1.0e-6
    end

    test "10 solar masses → ≈ 29_532 m" do
      assert_in_delta BlackHole.schwarzschild_radius_solar(10), 29_532.5, 1.0
    end
  end

  describe "hawking_temperature/1" do
    test "solar mass BH is extremely cold (~6.2e-8 K)" do
      t = BlackHole.hawking_temperature(@solar_mass)
      assert t > 0
      assert t < 1.0e-6
    end

    test "smaller mass → higher temperature (inverse relationship)" do
      t_solar = BlackHole.hawking_temperature(@solar_mass)
      t_small = BlackHole.hawking_temperature(1.0e20)
      assert t_small > t_solar
    end
  end

  describe "evaporation_time/1" do
    test "positive for any positive mass" do
      assert BlackHole.evaporation_time(1.0e12) > 0
    end

    test "scales as M³" do
      t1 = BlackHole.evaporation_time(1.0e10)
      t2 = BlackHole.evaporation_time(2.0e10)
      assert_in_delta t2 / t1, 8.0, 0.001
    end
  end

  describe "photon_sphere_radius/1" do
    test "equals 1.5 × Schwarzschild radius" do
      r_s = BlackHole.schwarzschild_radius(@solar_mass)
      r_ph = BlackHole.photon_sphere_radius(@solar_mass)
      assert_in_delta r_ph / r_s, 1.5, 1.0e-10
    end
  end

  describe "isco_radius/1" do
    test "equals 3 × Schwarzschild radius" do
      r_s = BlackHole.schwarzschild_radius(@solar_mass)
      r_isco = BlackHole.isco_radius(@solar_mass)
      assert_in_delta r_isco / r_s, 3.0, 1.0e-10
    end

    test "ISCO > photon sphere" do
      r_ph = BlackHole.photon_sphere_radius(@solar_mass)
      r_isco = BlackHole.isco_radius(@solar_mass)
      assert r_isco > r_ph
    end
  end

  describe "bekenstein_hawking_entropy/1" do
    test "positive for positive mass" do
      assert BlackHole.bekenstein_hawking_entropy(@solar_mass) > 0
    end

    test "scales as M² (area ∝ r_s² ∝ M²)" do
      s1 = BlackHole.bekenstein_hawking_entropy(@solar_mass)
      s2 = BlackHole.bekenstein_hawking_entropy(2 * @solar_mass)
      assert_in_delta s2 / s1, 4.0, 0.001
    end
  end

  describe "kerr_spin_parameter/2" do
    test "zero angular momentum → a* = 0" do
      assert BlackHole.kerr_spin_parameter(@solar_mass, 0.0) == 0.0
    end

    test "very large J is clamped to 1" do
      assert BlackHole.kerr_spin_parameter(@solar_mass, 1.0e100) == 1.0
    end

    test "result is in [0, 1]" do
      a = BlackHole.kerr_spin_parameter(@solar_mass, 1.0e47)
      assert a >= 0.0 and a <= 1.0
    end
  end

  describe "gravitational_time_dilation/2" do
    test "factor is 1 at infinity (very large r)" do
      assert_in_delta BlackHole.gravitational_time_dilation(@solar_mass, 1.0e20), 1.0, 1.0e-6
    end

    test "factor is 0 at event horizon" do
      r_s = BlackHole.schwarzschild_radius(@solar_mass)
      assert BlackHole.gravitational_time_dilation(@solar_mass, r_s) == 0.0
    end

    test "factor is between 0 and 1 outside horizon" do
      f = BlackHole.gravitational_time_dilation(@solar_mass, 6.957e8)
      assert f > 0.0 and f < 1.0
    end
  end
end
