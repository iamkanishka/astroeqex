defmodule AstroEquations.AstrophysicsAndAstronomy.BlackHoleTest do
  use ExUnit.Case
  import AstroEquations.AstrophysicsAndAstronomy.BlackHole

  describe "schwarzschild_radius/1" do
    test "calculates for solar mass" do
      # Expected value for 1 solar mass: ~2953 meters
      assert schwarzschild_radius(1.989e30) |> Float.round(4) == 2953.2501
    end

    test "radius increases linearly with mass" do
      r1 = schwarzschild_radius(1.989e30)
      r2 = schwarzschild_radius(2 * 1.989e30)
      assert Float.round(r2/r1, 4) == 2.0
    end
  end

  describe "schwarzschild_radius_solar/1" do
    test "calculates for 1 solar mass" do
      assert schwarzschild_radius_solar(1) |> Float.round(4) == 2953.2501
    end

    test "calculates for 10 solar masses" do
      assert schwarzschild_radius_solar(10) |> Float.round(4) == 29532.5008
    end
  end
end
