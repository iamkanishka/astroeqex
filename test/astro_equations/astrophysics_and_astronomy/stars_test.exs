defmodule AstroEquations.AstrophysicsAndAstronomy.StarsTest do
  use ExUnit.Case, async: true

  alias AstroEquations.AstrophysicsAndAstronomy.Stars

  @solar_mass 1.989e30
  @solar_radius 6.957e8
  @solar_luminosity 3.828e26
  @g 6.67430e-11

  describe "hydrostatic_equilibrium/3" do
    test "returns negative value (inward pressure gradient)" do
      result = Stars.hydrostatic_equilibrium(@solar_mass, 1408.0, @solar_radius)
      assert result < 0
    end

    test "larger mass increases magnitude" do
      r1 = Stars.hydrostatic_equilibrium(@solar_mass, 1408.0, @solar_radius)
      r2 = Stars.hydrostatic_equilibrium(2 * @solar_mass, 1408.0, @solar_radius)
      assert abs(r2) > abs(r1)
    end
  end

  describe "mass_conservation/2" do
    test "always positive" do
      assert Stars.mass_conservation(1408.0, @solar_radius) > 0
    end

    test "scales as r²" do
      r1 = Stars.mass_conservation(1000.0, 1.0e8)
      r2 = Stars.mass_conservation(1000.0, 2.0e8)
      assert_in_delta r2 / r1, 4.0, 1.0e-10
    end
  end

  describe "energy_equation/3" do
    test "returns positive luminosity gradient" do
      result = Stars.energy_equation(1408.0, 1.934e-7, @solar_radius)
      assert result > 0
    end
  end

  describe "kelvin_helmholtz_timescale/3" do
    test "solar KH timescale is ~30 Myr" do
      t = Stars.kelvin_helmholtz_timescale(1.0, 1.0, 1.0)
      assert t > 1.0e7
      assert t < 1.0e8
    end

    test "more luminous star has shorter timescale" do
      t1 = Stars.kelvin_helmholtz_timescale(1.0, 1.0, 1.0)
      t10 = Stars.kelvin_helmholtz_timescale(1.0, 1.0, 10.0)
      assert t10 < t1
    end
  end

  describe "nuclear_timescale/1" do
    test "solar-mass star lives ~1 Gyr" do
      assert_in_delta Stars.nuclear_timescale(1.0), 1.0e9, 1.0e6
    end

    test "more massive star has shorter timescale (t ∝ M⁻³)" do
      t1 = Stars.nuclear_timescale(1.0)
      t2 = Stars.nuclear_timescale(2.0)
      assert_in_delta t2 / t1, 0.125, 0.001
    end
  end

  describe "dynamical_timescale/2" do
    test "solar dynamical timescale is positive" do
      assert Stars.dynamical_timescale(1.0, 1.0) > 0
    end

    test "timescale in seconds (shorter than KH)" do
      t_dyn = Stars.dynamical_timescale(1.0, 1.0)
      t_kh = Stars.kelvin_helmholtz_timescale(1.0, 1.0, 1.0) * 3.15576e7
      assert t_dyn < t_kh
    end
  end

  describe "eddington_luminosity/1" do
    test "solar Eddington luminosity is ~32_000 L☉" do
      assert_in_delta Stars.eddington_luminosity(1.0), 32_000.0, 1000.0
    end

    test "scales linearly with mass" do
      l1 = Stars.eddington_luminosity(1.0)
      l10 = Stars.eddington_luminosity(10.0)
      assert_in_delta l10 / l1, 10.0, 0.01
    end
  end

  describe "mass_luminosity/1" do
    test "1 M☉ → 1 L☉" do
      assert_in_delta Stars.mass_luminosity(1.0), 1.0, 0.01
    end

    test "10 M☉ → ~3162 L☉" do
      assert_in_delta Stars.mass_luminosity(10.0), 3162.3, 1.0
    end

    test "luminosity increases with mass throughout" do
      l1 = Stars.mass_luminosity(1.0)
      l5 = Stars.mass_luminosity(5.0)
      assert l5 > l1
    end

    test "0.3 M☉ (low-mass branch)" do
      assert Stars.mass_luminosity(0.3) < 1.0
    end
  end

  describe "chandrasekhar_mass/1" do
    test "default (C/O WD, μ_e=2) → ~1.46 M☉" do
      assert_in_delta Stars.chandrasekhar_mass(), 1.4575, 0.01
    end

    test "higher mean molecular weight → lower limit" do
      assert Stars.chandrasekhar_mass(2.0) > Stars.chandrasekhar_mass(4.0)
    end
  end

  describe "wien_peak_wavelength/1" do
    test "solar surface (5778 K) peaks near 500 nm" do
      lambda = Stars.wien_peak_wavelength(5778)
      assert_in_delta lambda, 5.015e-7, 1.0e-9
    end

    test "hotter star peaks at shorter wavelength" do
      assert Stars.wien_peak_wavelength(10_000) < Stars.wien_peak_wavelength(5000)
    end
  end

  describe "stefan_boltzmann_luminosity/2" do
    test "solar values give ~1 L☉" do
      l = Stars.stefan_boltzmann_luminosity(@solar_radius, 5778)
      assert_in_delta l / @solar_luminosity, 1.0, 0.05
    end

    test "doubles with sqrt(2) increase in radius (at fixed T)" do
      l1 = Stars.stefan_boltzmann_luminosity(1.0e9, 5000)
      l2 = Stars.stefan_boltzmann_luminosity(:math.sqrt(2) * 1.0e9, 5000)
      assert_in_delta l2 / l1, 2.0, 0.001
    end
  end

  describe "jeans_mass/3" do
    test "typical molecular cloud gives positive Jeans mass" do
      assert Stars.jeans_mass(10, 1.0e-17) > 0
    end

    test "higher temperature → larger Jeans mass" do
      m10 = Stars.jeans_mass(10, 1.0e-17)
      m30 = Stars.jeans_mass(30, 1.0e-17)
      assert m30 > m10
    end

    test "higher density → smaller Jeans mass" do
      m_lo = Stars.jeans_mass(10, 1.0e-17)
      m_hi = Stars.jeans_mass(10, 1.0e-16)
      assert m_hi < m_lo
    end
  end
end
