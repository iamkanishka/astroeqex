defmodule AstroEquations.AstrophysicsAndAstronomy.AstrometryTest do
  use ExUnit.Case
  alias AstroEquations.AstrophysicsAndAstronomy.Astrometry

  describe "redshift/2" do
    test "calculates redshift from wavelength difference" do
      assert Astrometry.redshift(700, 656) |> Float.round(6) == 0.067073
    end
  end

  describe "redshift_ratio/2" do
    test "calculates redshift from wavelength ratio" do
      assert Astrometry.redshift_ratio(700, 656) |> Float.round(6) == 0.067073
    end
  end

  describe "apparent_magnitude_diff/4" do
    test "calculates magnitude difference from flux ratio" do
      # Test with known values where F/F0 = 1/6.31 should give Δm = 2.5
      assert Astrometry.apparent_magnitude_diff(12, 10, 1.0, 6.31) |> Float.round(6) == -2.5
    end
  end

  describe "absolute_magnitude/2" do
    test "calculates absolute magnitude from apparent magnitude and distance" do
      assert Astrometry.absolute_magnitude(5, 100) == 0.0
      assert Astrometry.absolute_magnitude(10, 10) == 10.0
    end
  end

  describe "flux_ratio_from_magnitudes/2" do
    test "calculates flux ratio from magnitude difference" do
      assert Astrometry.flux_ratio_from_magnitudes(10, 15) == 100.0
      assert Astrometry.flux_ratio_from_magnitudes(15, 10) |> Float.round(6) == 0.01
    end
  end

  describe "flux_from_magnitude/2" do
    test "calculates flux from magnitude" do
      assert Astrometry.flux_from_magnitude(0, 1.0) == 1.0
      assert Astrometry.flux_from_magnitude(1, 1.0) |> Float.round(4) == 0.3981
    end
  end

  describe "color_index/2" do
    test "calculates color index from two fluxes" do
      assert Astrometry.color_index(1.0, 2.0) |> Float.round(4) == 0.7526
      assert Astrometry.color_index(2.0, 1.0) |> Float.round(4) == -0.7526
    end
  end

  describe "metallicity/2" do
    test "calculates metallicity relative to solar" do
      assert Astrometry.metallicity(0.316, 1.0) |> Float.round(6) == -0.500000
      assert Astrometry.metallicity(1.0, 1.0) == 0.0
    end
  end
end
