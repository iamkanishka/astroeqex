defmodule AstroEquations.Mathematics.GeometryTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Mathematics.Geometry

  @pi :math.pi()
  @tol 1.0e-8

  describe "deg_to_rad/1 and rad_to_deg/1" do
    test "0° = 0 rad" do
      assert_in_delta Geometry.deg_to_rad(0), 0.0, @tol
    end

    test "180° = π rad" do
      assert_in_delta Geometry.deg_to_rad(180), @pi, @tol
    end

    test "round trip" do
      assert_in_delta Geometry.rad_to_deg(Geometry.deg_to_rad(45)), 45.0, @tol
    end
  end

  describe "hours_to_deg/1 and deg_to_hours/1" do
    test "12h = 180°" do
      assert_in_delta Geometry.hours_to_deg(12), 180.0, @tol
    end

    test "round trip" do
      assert_in_delta Geometry.deg_to_hours(Geometry.hours_to_deg(6)), 6.0, @tol
    end
  end

  describe "arcsec_to_rad/1 and rad_to_arcsec/1" do
    test "206_265 arcsec = 1 rad" do
      assert_in_delta Geometry.arcsec_to_rad(206_265.0), 1.0, 1.0e-5
    end

    test "round trip" do
      r = 0.001
      assert_in_delta Geometry.arcsec_to_rad(Geometry.rad_to_arcsec(r)), r, @tol
    end
  end

  describe "angular_separation/4" do
    test "same point: separation = 0" do
      assert_in_delta Geometry.angular_separation(1.0, 0.5, 1.0, 0.5), 0.0, @tol
    end

    test "π/2 separation on equator" do
      assert_in_delta Geometry.angular_separation(0.0, 0.0, @pi / 2, 0.0), @pi / 2, @tol
    end

    test "antipodal points: separation = π" do
      assert_in_delta Geometry.angular_separation(0.0, 0.0, @pi, 0.0), @pi, @tol
    end
  end

  describe "angular_separation_deg/4" do
    test "same point: 0°" do
      assert_in_delta Geometry.angular_separation_deg(45.0, 30.0, 45.0, 30.0), 0.0, @tol
    end

    test "result in degrees" do
      sep_deg = Geometry.angular_separation_deg(0.0, 0.0, 90.0, 0.0)
      assert_in_delta sep_deg, 90.0, 0.001
    end
  end

  describe "position_angle/4" do
    test "due east: position angle = π/2" do
      pa = Geometry.position_angle(0.0, 0.0, 0.1, 0.0)
      assert_in_delta pa, @pi / 2, 0.01
    end

    test "result in [0, 2π)" do
      pa = Geometry.position_angle(1.0, 0.5, 1.1, 0.6)
      assert pa >= 0.0 and pa < 2 * @pi
    end
  end

  describe "solid_angle_cone/1" do
    test "full hemisphere: θ = π/2 → Ω = 2π" do
      assert_in_delta Geometry.solid_angle_cone(@pi / 2), 2 * @pi, @tol
    end

    test "full sphere: θ → π → Ω → 4π" do
      assert_in_delta Geometry.solid_angle_cone(@pi), 4 * @pi, @tol
    end

    test "zero angle: Ω = 0" do
      assert_in_delta Geometry.solid_angle_cone(0), 0.0, @tol
    end
  end

  describe "periapsis/2 and apoapsis/2" do
    test "periapsis = a(1-e): a=1, e=0.5 → 0.5" do
      assert_in_delta Geometry.periapsis(1.0, 0.5), 0.5, @tol
    end

    test "apoapsis = a(1+e): a=1, e=0.5 → 1.5" do
      assert_in_delta Geometry.apoapsis(1.0, 0.5), 1.5, @tol
    end

    test "circular orbit: periapsis = apoapsis = a" do
      assert_in_delta Geometry.periapsis(1.5, 0.0), 1.5, @tol
      assert_in_delta Geometry.apoapsis(1.5, 0.0), 1.5, @tol
    end
  end

  describe "eccentric_anomaly/3" do
    test "M=0 → E=0 (periapsis)" do
      assert_in_delta Geometry.eccentric_anomaly(0.0, 0.5), 0.0, 1.0e-10
    end

    test "M=π → E=π (apoapsis)" do
      assert_in_delta Geometry.eccentric_anomaly(@pi, 0.5), @pi, 1.0e-8
    end

    test "satisfies Kepler's equation: M = E - e sin(E)" do
      m = 1.5
      e = 0.4
      big_e = Geometry.eccentric_anomaly(m, e)
      m_reconstructed = big_e - e * :math.sin(big_e)
      assert_in_delta m_reconstructed, m, 1.0e-8
    end
  end

  describe "true_anomaly/2" do
    test "E=0 → ν=0 (periapsis)" do
      assert_in_delta Geometry.true_anomaly(0.0, 0.5), 0.0, @tol
    end

    test "E=π → ν=π (apoapsis)" do
      assert_in_delta Geometry.true_anomaly(@pi, 0.5), @pi, 1.0e-8
    end
  end

  describe "orbital_radius/3" do
    test "at periapsis (ν=0): r = a(1-e)" do
      r = Geometry.orbital_radius(1.0, 0.5, 0.0)
      assert_in_delta r, Geometry.periapsis(1.0, 0.5), @tol
    end

    test "at apoapsis (ν=π): r = a(1+e)" do
      r = Geometry.orbital_radius(1.0, 0.5, @pi)
      assert_in_delta r, Geometry.apoapsis(1.0, 0.5), 1.0e-8
    end
  end

  describe "microlensing_magnification/1" do
    test "positive for all u > 0" do
      assert Geometry.microlensing_magnification(0.5) > 0
    end

    test "magnification → ∞ as u → 0" do
      mu_small = Geometry.microlensing_magnification(0.01)
      mu_large = Geometry.microlensing_magnification(1.0)
      assert mu_small > mu_large
    end

    test "A(u=1) ≈ 1.342" do
      assert_in_delta Geometry.microlensing_magnification(1.0), 1.342, 0.001
    end
  end

  describe "equatorial_to_cartesian/2 and cartesian_to_equatorial/3" do
    test "round-trip preserves RA and Dec" do
      ra = 1.5
      dec = 0.3
      {x, y, z} = Geometry.equatorial_to_cartesian(ra, dec)
      {ra2, dec2} = Geometry.cartesian_to_equatorial(x, y, z)
      assert_in_delta ra2, ra, @tol
      assert_in_delta dec2, dec, @tol
    end

    test "unit vector: |r| = 1" do
      {x, y, z} = Geometry.equatorial_to_cartesian(1.0, 0.5)
      norm = :math.sqrt(x * x + y * y + z * z)
      assert_in_delta norm, 1.0, @tol
    end
  end
end

defmodule AstroEquations.Mathematics.NotationTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Mathematics.Notation

  @tol 1.0e-8

  describe "physical constants" do
    test "speed of light ≈ 2.998×10⁸ m/s" do
      assert_in_delta Notation.speed_of_light(), 2.99792458e8, 1.0
    end

    test "gravitational constant" do
      assert_in_delta Notation.gravitational_constant(), 6.67430e-11, 1.0e-15
    end

    test "Boltzmann constant" do
      assert_in_delta Notation.boltzmann_constant(), 1.380649e-23, 1.0e-30
    end

    test "solar mass" do
      assert_in_delta Notation.solar_mass(), 1.989e30, 1.0e26
    end
  end

  describe "unit conversions" do
    test "pc_to_m: 1 pc ≈ 3.086×10¹⁶ m" do
      assert_in_delta Notation.pc_to_m(1), 3.085677581e16, 1.0e8
    end

    test "au_to_m: 1 AU ≈ 1.496×10¹¹ m" do
      assert_in_delta Notation.au_to_m(1), 1.495978707e11, 1.0e4
    end

    test "ev_to_j: 1 eV ≈ 1.602×10⁻¹⁹ J" do
      assert_in_delta Notation.ev_to_j(1), 1.602176634e-19, 1.0e-26
    end

    test "j_to_ev round-trips ev_to_j" do
      j = 1.0e-18
      assert_in_delta Notation.ev_to_j(Notation.j_to_ev(j)), j, j * @tol
    end

    test "deg_to_rad: 180° = π" do
      assert_in_delta Notation.deg_to_rad(180), :math.pi(), @tol
    end

    test "rad_to_deg: π = 180°" do
      assert_in_delta Notation.rad_to_deg(:math.pi()), 180.0, @tol
    end

    test "k_to_celsius: 273.15 K = 0°C" do
      assert_in_delta Notation.k_to_celsius(273.15), 0.0, @tol
    end

    test "celsius_to_k: 0°C = 273.15 K" do
      assert_in_delta Notation.celsius_to_k(0.0), 273.15, @tol
    end
  end

  describe "sexagesimal conversions" do
    test "hms_to_deg: 12h 30m 0s = 187.5°" do
      assert_in_delta Notation.hms_to_deg(12, 30, 0), 187.5, @tol
    end

    test "dms_to_deg: 30°15'0\" = 30.25°" do
      assert_in_delta Notation.dms_to_deg(30, 15, 0), 30.25, @tol
    end

    test "dms_to_deg: negative declination" do
      assert Notation.dms_to_deg(-30, 0, 0) < 0
    end

    test "deg_to_hms round-trips hms_to_deg" do
      h = 5
      m = 34
      s = 32.0
      deg = Notation.hms_to_deg(h, m, s)
      {h2, m2, s2} = Notation.deg_to_hms(deg)
      assert h2 == h
      assert m2 == m
      assert_in_delta s2, s, 0.01
    end

    test "deg_to_dms round-trips dms_to_deg" do
      d = 22
      m = 0
      s = 52.0
      deg = Notation.dms_to_deg(d, m, s)
      {d2, m2, s2} = Notation.deg_to_dms(deg)
      assert d2 == d
      assert m2 == m
      assert_in_delta s2, s, 0.01
    end
  end

  describe "Julian Date" do
    test "J2000.0 = 2_451_545.0" do
      assert_in_delta Notation.j2000(), 2_451_545.0, @tol
    end

    test "calendar_to_jd: 2000-01-01 12:00 UTC = J2000.0" do
      jd = Notation.calendar_to_jd(2000, 1, 1, 12, 0, 0)
      assert_in_delta jd, 2_451_545.0, @tol
    end

    test "jd_to_mjd: MJD = JD - 2_400_000.5" do
      jd = 2_451_545.0
      assert_in_delta Notation.jd_to_mjd(jd), 51_544.5, @tol
    end

    test "mjd_to_jd round-trips jd_to_mjd" do
      jd = 2_451_545.0
      assert_in_delta Notation.mjd_to_jd(Notation.jd_to_mjd(jd)), jd, @tol
    end

    test "julian_centuries: J2000.0 → T = 0" do
      assert_in_delta Notation.julian_centuries(Notation.j2000()), 0.0, @tol
    end
  end
end
