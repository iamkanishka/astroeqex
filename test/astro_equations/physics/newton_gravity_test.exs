defmodule AstroEquations.Physics.NewtonGravityTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.NewtonGravity

  @g 6.67430e-11
  @m_earth 5.972e24
  @r_earth 6.371e6
  @m_sun 1.989e30
  @m_moon 7.348e22
  @r_moon 3.844e8

  describe "force/3" do
    test "Earth-Moon gravitational force" do
      f = NewtonGravity.force(@m_earth, @m_moon, @r_moon)
      assert_in_delta f, 1.982e20, 1.0e18
    end

    test "force is symmetric" do
      f1 = NewtonGravity.force(1.0e24, 7.0e22, 3.8e8)
      f2 = NewtonGravity.force(7.0e22, 1.0e24, 3.8e8)
      assert_in_delta f1, f2, 1.0
    end

    test "force increases as r decreases (inverse square)" do
      f1 = NewtonGravity.force(1.0e24, 1.0e20, 1.0e8)
      f2 = NewtonGravity.force(1.0e24, 1.0e20, 2.0e8)
      assert_in_delta f1 / f2, 4.0, 0.001
    end

    test "raises on non-positive distance" do
      assert_raise ArgumentError, fn -> NewtonGravity.force(1.0e20, 1.0e20, 0) end
    end
  end

  describe "field/2" do
    test "Earth's surface gravity ≈ 9.82 m/s²" do
      g = NewtonGravity.field(@m_earth, @r_earth)
      assert_in_delta g, 9.82, 0.05
    end

    test "inverse square: field at 2r is ¼ of field at r" do
      g1 = NewtonGravity.field(1.0e24, 1.0e7)
      g2 = NewtonGravity.field(1.0e24, 2.0e7)
      assert_in_delta g1 / g2, 4.0, 0.001
    end
  end

  describe "potential/2" do
    test "gravitational potential is always negative" do
      assert NewtonGravity.potential(@m_earth, @r_earth) < 0
    end

    test "potential increases (becomes less negative) with distance" do
      phi_r = NewtonGravity.potential(@m_earth, @r_earth)
      phi_2r = NewtonGravity.potential(@m_earth, 2 * @r_earth)
      assert phi_2r > phi_r
    end
  end

  describe "potential_energy/3" do
    test "Earth-Moon potential energy is negative" do
      u = NewtonGravity.potential_energy(@m_earth, @m_moon, @r_moon)
      assert u < 0
    end

    test "inversely proportional to separation" do
      u1 = NewtonGravity.potential_energy(1.0e24, 1.0e20, 1.0e8)
      u2 = NewtonGravity.potential_energy(1.0e24, 1.0e20, 2.0e8)
      assert_in_delta u1 / u2, 2.0, 0.001
    end
  end

  describe "orbital_period/3" do
    test "Moon's orbital period ≈ 27.3 days" do
      t = NewtonGravity.orbital_period(@r_moon, @m_earth, @m_moon)
      days = t / 86_400
      assert_in_delta days, 27.3, 0.5
    end

    test "scales as a^(3/2)" do
      t1 = NewtonGravity.orbital_period(1.0e8, 1.0e24, 0)
      t2 = NewtonGravity.orbital_period(2.0e8, 1.0e24, 0)
      assert_in_delta t2 / t1, 2.0 ** 1.5, 0.001
    end
  end

  describe "escape_velocity/2" do
    test "Earth escape velocity ≈ 11.2 km/s" do
      v_esc = NewtonGravity.escape_velocity(@m_earth, @r_earth)
      assert_in_delta v_esc / 1000, 11.2, 0.1
    end

    test "escape velocity decreases with distance" do
      v1 = NewtonGravity.escape_velocity(@m_earth, @r_earth)
      v2 = NewtonGravity.escape_velocity(@m_earth, 2 * @r_earth)
      assert v2 < v1
    end
  end

  describe "circular_orbit_speed/2" do
    test "ISS orbit (~400 km): ~7.7 km/s" do
      r_iss = @r_earth + 4.0e5
      v = NewtonGravity.circular_orbit_speed(@m_earth, r_iss)
      assert_in_delta v / 1000, 7.7, 0.1
    end

    test "v_circular < v_escape at same radius" do
      v_c = NewtonGravity.circular_orbit_speed(@m_earth, @r_earth)
      v_e = NewtonGravity.escape_velocity(@m_earth, @r_earth)
      assert v_e > v_c
      assert_in_delta v_e / v_c, :math.sqrt(2), 0.001
    end
  end

  describe "surface_gravity/2" do
    test "Earth's surface gravity ≈ 9.82 m/s²" do
      g = NewtonGravity.surface_gravity(@m_earth, @r_earth)
      assert_in_delta g, 9.82, 0.05
    end

    test "Moon's surface gravity ≈ 1.62 m/s²" do
      g = NewtonGravity.surface_gravity(@m_moon, 1.737e6)
      assert_in_delta g, 1.62, 0.05
    end
  end

  describe "tidal_acceleration/3" do
    test "positive value (differential acceleration)" do
      assert NewtonGravity.tidal_acceleration(@m_moon, @r_moon, @r_earth) > 0
    end

    test "scales as 1/r³" do
      a1 = NewtonGravity.tidal_acceleration(1.0e22, 1.0e9, 1.0e6)
      a2 = NewtonGravity.tidal_acceleration(1.0e22, 2.0e9, 1.0e6)
      assert_in_delta a1 / a2, 8.0, 0.001
    end
  end

  describe "roche_limit/3" do
    test "positive value" do
      assert NewtonGravity.roche_limit(1.737e6, @m_earth, @m_moon) > 0
    end

    test "Moon's Roche limit from Earth ≈ 9500 km" do
      d = NewtonGravity.roche_limit(1.737e6, @m_earth, @m_moon)
      assert_in_delta d / 1000, 9500, 500
    end
  end

  describe "hill_sphere/3" do
    test "Earth's Hill sphere ≈ 1.5 million km" do
      h = NewtonGravity.hill_sphere(1.496e11, @m_earth, @m_sun)
      assert_in_delta h / 1.0e9, 1500.0, 100.0
    end
  end

  describe "vis_viva/3" do
    test "circular orbit: a = r → v = √(GM/r)" do
      m = @m_earth
      r = @r_earth + 4.0e5
      v_vv = NewtonGravity.vis_viva(m, r, r)
      v_c = NewtonGravity.circular_orbit_speed(m, r)
      assert_in_delta v_vv, v_c, 1.0
    end
  end
end
