defmodule AstroEquations.AstrophysicsAndAstronomy.AstrometryTest do
  use ExUnit.Case, async: true
  doctest AstroEquations.AstrophysicsAndAstronomy.Astrometry

  alias AstroEquations.AstrophysicsAndAstronomy.Astrometry

  # ---------------------------------------------------------------------------
  # Redshift
  # ---------------------------------------------------------------------------

  describe "redshift/2" do
    test "Hα emission line shifted from 656 nm to 700 nm" do
      assert_in_delta Astrometry.redshift(700, 656), 0.0671, 0.0001
    end

    test "no shift returns zero" do
      assert Astrometry.redshift(500, 500) == 0.0
    end

    test "blueshift gives negative z" do
      assert Astrometry.redshift(600, 656) < 0
    end
  end

  describe "redshift_ratio/2" do
    test "is numerically equivalent to redshift/2" do
      assert_in_delta Astrometry.redshift_ratio(700, 656),
                      Astrometry.redshift(700, 656),
                      1.0e-12
    end
  end

  describe "recession_velocity/2" do
    test "z = 0.1 gives ~30_000 km/s" do
      assert_in_delta Astrometry.recession_velocity(0.1), 29_979.2458, 0.001
    end

    test "zero redshift means zero velocity" do
      assert Astrometry.recession_velocity(0.0) == 0.0
    end
  end

  describe "relativistic_velocity/2" do
    test "low-z approaches non-relativistic limit" do
      v_rel = Astrometry.relativistic_velocity(0.01)
      v_nr = Astrometry.recession_velocity(0.01)
      assert_in_delta v_rel, v_nr, 50.0
    end

    test "z = 1 gives v ≈ 0.6c in km/s" do
      assert_in_delta Astrometry.relativistic_velocity(1.0), 179_875.47, 1.0
    end
  end

  # ---------------------------------------------------------------------------
  # Magnitude & Flux
  # ---------------------------------------------------------------------------

  describe "apparent_magnitude_diff/2" do
    test "equal fluxes give zero magnitude difference" do
      assert_in_delta Astrometry.apparent_magnitude_diff(1.0, 1.0), 0.0, 1.0e-10
    end

    test "flux ratio 100 gives 5 mag difference" do
      assert_in_delta Astrometry.apparent_magnitude_diff(1.0, 100.0), -5.0, 1.0e-10
    end
  end

  describe "absolute_magnitude/2" do
    test "object at 10 pc: M = m" do
      assert_in_delta Astrometry.absolute_magnitude(5.0, 10.0), 5.0, 1.0e-10
    end

    test "object at 100 pc: M = m - 5" do
      assert_in_delta Astrometry.absolute_magnitude(10.0, 100.0), 5.0, 1.0e-10
    end
  end

  describe "distance_modulus/1" do
    test "10 pc gives μ = 0" do
      assert_in_delta Astrometry.distance_modulus(10), 0.0, 1.0e-10
    end

    test "100 pc gives μ = 5" do
      assert_in_delta Astrometry.distance_modulus(100), 5.0, 1.0e-10
    end

    test "1000 pc gives μ = 10" do
      assert_in_delta Astrometry.distance_modulus(1000), 10.0, 1.0e-10
    end
  end

  describe "distance_from_modulus/1" do
    test "μ = 0 gives 10 pc" do
      assert_in_delta Astrometry.distance_from_modulus(0.0), 10.0, 1.0e-6
    end

    test "μ = 5 gives 100 pc" do
      assert_in_delta Astrometry.distance_from_modulus(5.0), 100.0, 1.0e-6
    end

    test "round-trips with distance_modulus" do
      d = 3500.0
      assert_in_delta Astrometry.distance_from_modulus(Astrometry.distance_modulus(d)), d, 1.0e-6
    end
  end

  describe "flux_ratio_from_magnitudes/2" do
    test "5 mag difference gives flux ratio 100" do
      assert_in_delta Astrometry.flux_ratio_from_magnitudes(10, 15), 100.0, 1.0e-8
    end

    test "0 mag difference gives flux ratio 1" do
      assert_in_delta Astrometry.flux_ratio_from_magnitudes(10, 10), 1.0, 1.0e-10
    end
  end

  describe "flux_from_magnitude/2" do
    test "magnitude 0 gives F0" do
      assert_in_delta Astrometry.flux_from_magnitude(0, 1.0), 1.0, 1.0e-10
    end

    test "magnitude 5 gives F0 × 10^(-2)" do
      expected = 1.0 * 10 ** (-0.4 * 5)
      assert_in_delta Astrometry.flux_from_magnitude(5, 1.0), expected, 1.0e-10
    end
  end

  describe "color_index/2" do
    test "equal fluxes give zero colour index" do
      assert_in_delta Astrometry.color_index(1.0, 1.0), 0.0, 1.0e-10
    end

    test "colour index is positive when F_f1 < F_f2" do
      assert Astrometry.color_index(0.5, 1.0) > 0
    end
  end

  # ---------------------------------------------------------------------------
  # Parallax & Proper Motion
  # ---------------------------------------------------------------------------

  describe "distance_from_parallax/1" do
    test "1 arcsec parallax → 1 pc" do
      assert_in_delta Astrometry.distance_from_parallax(1.0), 1.0, 1.0e-10
    end

    test "0.1 arcsec parallax → 10 pc" do
      assert_in_delta Astrometry.distance_from_parallax(0.1), 10.0, 1.0e-10
    end

    test "Proxima Centauri: π ≈ 0.7687 arcsec → ~1.30 pc" do
      assert_in_delta Astrometry.distance_from_parallax(0.7687), 1.3009, 0.001
    end
  end

  describe "parallax_from_distance/1" do
    test "1 pc → 1 arcsec parallax" do
      assert_in_delta Astrometry.parallax_from_distance(1.0), 1.0, 1.0e-10
    end

    test "round-trips with distance_from_parallax" do
      p = 0.314

      assert_in_delta Astrometry.parallax_from_distance(Astrometry.distance_from_parallax(p)),
                      p,
                      1.0e-10
    end
  end

  describe "proper_motion_total/2" do
    test "3, 4 arcsec/yr → 5 arcsec/yr (Pythagorean triple)" do
      assert_in_delta Astrometry.proper_motion_total(3.0, 4.0), 5.0, 1.0e-10
    end

    test "zero proper motion" do
      assert Astrometry.proper_motion_total(0.0, 0.0) == 0.0
    end
  end

  describe "transverse_velocity/2" do
    test "μ = 1 arcsec/yr at 1 pc → 4.74 km/s" do
      assert_in_delta Astrometry.transverse_velocity(1.0, 1.0), 4.74, 0.01
    end
  end

  describe "space_velocity/2" do
    test "3, 4 km/s components → 5 km/s total" do
      assert_in_delta Astrometry.space_velocity(3.0, 4.0), 5.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Angular sizes & distances
  # ---------------------------------------------------------------------------

  describe "angular_diameter/2" do
    test "1 m at 206_265 m → 1 arcsec" do
      assert_in_delta Astrometry.angular_diameter(1.0, 206_265.0), 1.0, 1.0e-6
    end
  end

  describe "physical_size/2" do
    test "inverse of angular_diameter" do
      assert_in_delta Astrometry.physical_size(1.0, 206_265.0), 1.0, 1.0e-6
    end
  end

  describe "hubble_distance/2" do
    test "z = 0 → zero distance" do
      assert Astrometry.hubble_distance(0.0) == 0.0
    end

    test "z = 0.1, H₀ = 70 → ~428 Mpc" do
      assert_in_delta Astrometry.hubble_distance(0.1, 70.0), 428.28, 0.5
    end
  end

  # ---------------------------------------------------------------------------
  # Metallicity
  # ---------------------------------------------------------------------------

  describe "metallicity/2" do
    test "solar ratio gives [Fe/H] = 0" do
      assert_in_delta Astrometry.metallicity(0.017, 0.017), 0.0, 1.0e-10
    end

    test "half solar Fe/H gives negative metallicity" do
      assert Astrometry.metallicity(0.0085, 0.017) < 0
    end
  end

  describe "feh_to_z/1" do
    test "[Fe/H] = 0 → Z ≈ 0.017 (solar)" do
      assert_in_delta Astrometry.feh_to_z(0.0), 0.017, 1.0e-6
    end

    test "[Fe/H] = 1 → Z ≈ 0.17" do
      assert_in_delta Astrometry.feh_to_z(1.0), 0.17, 0.001
    end
  end
end
