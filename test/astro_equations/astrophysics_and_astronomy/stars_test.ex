defmodule AstroEquations.AstrophysicsAndAstronomy.StarsTest do
  use ExUnit.Case
  import AstroEquations.AstrophysicsAndAstronomy.Stars

  # Test values approximating conditions at the surface of the Sun
  # kg
  @solar_mass 1.989e30
  # m
  @solar_radius 6.957e8
  # kg/m^3 (estimated)
  @surface_density 1.408e3
  # K
  @surface_temp 5.778e3
  # W
  @solar_luminosity 3.828e26
  # W/kg (estimated)
  @energy_gen_rate 1.934e-7
  # m^2/kg (estimated)
  @opacity 0.04

  describe "hydrostatic_equilibrium/3" do
    test "calculates pressure gradient at solar surface" do
      result = hydrostatic_equilibrium(@solar_mass, @surface_density, @solar_radius)
      assert_in_delta result, -3.445e4, 0.001e4
    end

    test "gradient becomes more negative with higher mass" do
      normal = hydrostatic_equilibrium(@solar_mass, @surface_density, @solar_radius)
      higher_mass = hydrostatic_equilibrium(2 * @solar_mass, @surface_density, @solar_radius)
      assert higher_mass < normal
    end
  end

  describe "mass_conservation/2" do
    test "calculates mass gradient at solar surface" do
      result = mass_conservation(@surface_density, @solar_radius)
      assert_in_delta result, 8.571e19, 0.001e19
    end

    test "gradient increases with higher density" do
      normal = mass_conservation(@surface_density, @solar_radius)
      higher_density = mass_conservation(2 * @surface_density, @solar_radius)
      assert higher_density > normal
    end
  end

  describe "energy_equation/3" do
    test "calculates luminosity gradient at solar surface" do
      result = energy_equation(@surface_density, @energy_gen_rate, @solar_radius)
      assert_in_delta result, 8.231e7, 0.001e7
    end

    test "gradient increases with higher energy generation" do
      normal = energy_equation(@surface_density, @energy_gen_rate, @solar_radius)
      higher_energy = energy_equation(@surface_density, 2 * @energy_gen_rate, @solar_radius)
      assert higher_energy > normal
    end
  end

  describe "radiative_transport/5" do
    test "calculates temperature gradient at solar surface" do
      result =
        radiative_transport(
          @opacity,
          @surface_density,
          @surface_temp,
          @solar_luminosity,
          @solar_radius
        )

      assert_in_delta result, -0.01034, 0.00001
    end

    test "gradient becomes steeper with higher opacity" do
      normal =
        radiative_transport(
          @opacity,
          @surface_density,
          @surface_temp,
          @solar_luminosity,
          @solar_radius
        )

      higher_opacity =
        radiative_transport(
          2 * @opacity,
          @surface_density,
          @surface_temp,
          @solar_luminosity,
          @solar_radius
        )

      assert higher_opacity < normal
    end
  end


   describe "kelvin_helmholtz_timescale/3" do
    test "calculates KH timescale for Sun" do
      result = kelvin_helmholtz_timescale(1, 1, 1) |> round()
      assert result == 31_484_441 # ~31.5 million years
    end

    test "timescale decreases with higher luminosity" do
      normal = kelvin_helmholtz_timescale(1, 1, 1)
      higher_lum = kelvin_helmholtz_timescale(1, 1, 2)
      assert higher_lum < normal
    end
  end

  describe "nuclear_timescale/1" do
    test "calculates nuclear timescale for Sun" do
      assert nuclear_timescale(1) == 1.0e9
    end

    test "timescale decreases with mass^-3" do
      assert nuclear_timescale(2) == 1.0e9 / 8
    end
  end

  describe "gravitational_potential_energy/2" do
    test "calculates potential energy for Sun" do
      result = gravitational_potential_energy(1, 1) |> round()
      assert result == -3.7909166e41
    end

    test "energy becomes more negative with higher mass" do
      normal = gravitational_potential_energy(1, 1)
      higher_mass = gravitational_potential_energy(2, 1)
      assert higher_mass < normal
    end
  end

  describe "eddington_luminosity/1" do
    test "calculates Eddington luminosity for Sun" do
      assert eddington_luminosity(1) |> round() == 32000
    end

    test "luminosity scales linearly with mass" do
      assert eddington_luminosity(2) |> round() == 64000
    end
  end

  describe "eddington_mass/1" do
    test "calculates Eddington mass for Sun's Eddington luminosity" do
      assert eddington_mass(32000) |> Float.round(6) == 1.0
    end

    test "mass scales linearly with luminosity" do
      assert eddington_mass(64000) |> Float.round(6) == 2.0
    end
  end

  describe "eddington_mass_loss_rate/1" do
    test "calculates Eddington mass loss rate for Sun" do
      assert eddington_mass_loss_rate(1) |> Float.round(8) == 2.4e-8
    end

    test "rate scales linearly with mass" do
      assert eddington_mass_loss_rate(2) |> Float.round(8) == 4.8e-8
    end
  end

  describe "mass_luminosity/1" do
    test "low mass stars (M < 0.43)" do
      assert mass_luminosity(0.1) |> Float.round(6) == 0.000719
    end

    test "Sun-like stars (0.43 < M < 2)" do
      assert mass_luminosity(1) |> Float.round(4) == 1.0
    end

    test "intermediate mass stars (2 < M < 20)" do
      assert mass_luminosity(10) |> Float.round(4) == 3162.2777
    end

    test "very massive stars (M > 55)" do
      assert mass_luminosity(100) |> Float.round(4) == 32000.0
    end
  end

end
