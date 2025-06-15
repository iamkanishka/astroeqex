defmodule AstroEquations.AstrophysicsAndAstronomy.GalaxiesTest do
  use ExUnit.Case
  import AstroEquations.AstrophysicsAndAstronomy.Galaxies

  describe "hubble_classify/2" do
    test "perfectly round galaxy is E0" do
      assert hubble_classify(10, 10) == "E0"
    end

    test "moderately flattened galaxy" do
      assert hubble_classify(10, 7) == "E3"
    end

    test "maximum flattening is E7" do
      assert hubble_classify(10, 3) == "E7"
      assert hubble_classify(10, 2) == "E7"  # Even flatter still E7
    end
  end

  describe "sersic_profile/4" do
    test "central surface brightness" do
      assert sersic_profile(100, 0, 1, 1) |> Float.round(4) == 100.0
    end

    test "exponential profile (n=1) at r_e" do
      # At r_e, brightness should be i0/e
      assert sersic_profile(100, 1, 1, 1) |> Float.round(4) == 100.0
    end

    test "de Vaucouleurs profile (n=4)" do
      assert sersic_profile(100, 1, 1, 4) |> Float.round(4) == 100.0
      assert sersic_profile(100, 2, 1, 4) |> Float.round(4) == 13.5335
    end
  end

  describe "disk_density/5" do
    test "central density" do
      assert disk_density(1, 0, 0, 1, 1) == 1.0
    end

    test "exponential falloff in both dimensions" do
      # At 1 scale length/height, density should be 1/e in each dimension
      assert disk_density(1, 1, 0, 1, 1) |> Float.round(4) == 0.3679
      assert disk_density(1, 0, 1, 1, 1) |> Float.round(4) == 0.3679
      assert disk_density(1, 1, 1, 1, 1) |> Float.round(4) == 0.1353
    end

    test "different scale lengths" do
      assert disk_density(1, 2, 0, 2, 1) |> Float.round(4) == 0.3679
      assert disk_density(1, 0, 2, 1, 2) |> Float.round(4) == 0.3679
    end
  end
end
