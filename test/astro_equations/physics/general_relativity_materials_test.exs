defmodule AstroEquations.Physics.GeneralRelativityTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.GeneralRelativity, as: GR

  @g 6.67430e-11
  @c 2.99792458e8
  @pi :math.pi()
  @tol 1.0e-8

  # ---------------------------------------------------------------------------
  # Metrics
  # ---------------------------------------------------------------------------

  describe "minkowski_metric/0" do
    test "returns 4×4 matrix" do
      m = GR.minkowski_metric()
      assert length(m) == 4
      assert Enum.all?(m, fn row -> length(row) == 4 end)
    end

    test "signature: η₀₀ = -1, ηᵢᵢ = +1" do
      m = GR.minkowski_metric()
      assert Enum.at(Enum.at(m, 0), 0) == -1
      assert Enum.at(Enum.at(m, 1), 1) == 1
      assert Enum.at(Enum.at(m, 2), 2) == 1
      assert Enum.at(Enum.at(m, 3), 3) == 1
    end

    test "off-diagonal elements are zero" do
      m = GR.minkowski_metric()

      Enum.each(0..3, fn i ->
        Enum.each(0..3, fn j ->
          if i != j, do: assert(Enum.at(Enum.at(m, i), j) == 0)
        end)
      end)
    end
  end

  describe "minkowski_interval/5" do
    test "timelike event (dt=1, dx=0): ds² = -c²" do
      ds2 = GR.minkowski_interval(1, 0, 0, 0, 1)
      assert_in_delta ds2, -1.0, @tol
    end

    test "null (light-like): ds² = 0 when dx = c dt" do
      c = 1.0
      ds2 = GR.minkowski_interval(1.0, 1.0, 0.0, 0.0, c)
      assert_in_delta ds2, 0.0, @tol
    end

    test "spacelike event (dt=0): ds² > 0" do
      ds2 = GR.minkowski_interval(0, 1, 0, 0)
      assert ds2 > 0
    end
  end

  describe "schwarzschild_metric/5" do
    test "returns 4×4 matrix" do
      m = GR.schwarzschild_metric(1, 3, @pi / 2)
      assert length(m) == 4
    end

    test "at equator (θ=π/2): angular component g₃₃ = r²" do
      r = 10.0
      m = GR.schwarzschild_metric(1, r, @pi / 2)
      assert_in_delta Enum.at(Enum.at(m, 3), 3), r * r, @tol
    end

    test "far from mass (large r): approaches flat Minkowski" do
      # At r → ∞, g₀₀ → -1
      m = GR.schwarzschild_metric(1, 1.0e10, @pi / 2)
      g00 = Enum.at(Enum.at(m, 0), 0)
      assert_in_delta g00, -1.0, 1.0e-9
    end
  end

  describe "flrw_interval/4" do
    test "static universe (a=1): same as Minkowski for flat space" do
      ds2 = GR.flrw_interval(1.0, 1.0, 1.0, 1.0)
      ds2_mink = GR.minkowski_interval(1.0, 1.0, 0.0, 0.0, 1.0)
      assert_in_delta ds2, ds2_mink, @tol
    end

    test "expanding universe (a>1): comoving distances stretched" do
      ds2_a1 = GR.flrw_interval(0, 1.0, 1.0, 1.0)
      ds2_a2 = GR.flrw_interval(0, 2.0, 1.0, 1.0)
      assert ds2_a2 > ds2_a1
    end
  end

  # ---------------------------------------------------------------------------
  # Tensor Algebra
  # ---------------------------------------------------------------------------

  describe "inverse_metric/1" do
    test "Minkowski metric is its own inverse" do
      m = GR.minkowski_metric()
      inv = GR.inverse_metric(m)
      assert inv == m
    end

    test "diagonal element inversion" do
      diag = [[2, 0], [0, 3]]
      inv = GR.inverse_metric(diag)
      assert_in_delta Enum.at(Enum.at(inv, 0), 0), 0.5, @tol
      assert_in_delta Enum.at(Enum.at(inv, 1), 1), 1.0 / 3.0, @tol
    end
  end

  describe "four_vector_product/3" do
    test "timelike unit vector: a·a = -1 (Minkowski)" do
      m = GR.minkowski_metric()
      a = [1, 0, 0, 0]
      assert_in_delta GR.four_vector_product(a, a, m), -1.0, @tol
    end

    test "spacelike unit vector: a·a = +1 (Minkowski)" do
      m = GR.minkowski_metric()
      a = [0, 1, 0, 0]
      assert_in_delta GR.four_vector_product(a, a, m), 1.0, @tol
    end

    test "orthogonal four-vectors: a·b = 0" do
      m = GR.minkowski_metric()
      a = [1, 0, 0, 0]
      b = [0, 1, 0, 0]
      assert_in_delta GR.four_vector_product(a, b, m), 0.0, @tol
    end
  end

  # ---------------------------------------------------------------------------
  # GR Physical Effects
  # ---------------------------------------------------------------------------

  describe "gravitational_time_dilation/2" do
    test "at large r: factor approaches 1" do
      f = GR.gravitational_time_dilation(1.989e30, 1.0e15)
      assert_in_delta f, 1.0, 1.0e-6
    end

    test "at event horizon: factor = 0" do
      r_s = 2 * @g * 1.989e30 / @c ** 2
      assert GR.gravitational_time_dilation(1.989e30, r_s) == 0.0
    end

    test "factor is between 0 and 1 outside horizon" do
      f = GR.gravitational_time_dilation(1.989e30, 6.957e8)
      assert f > 0.0 and f < 1.0
    end
  end

  describe "orbital_precession/3" do
    test "Mercury perihelion precession: ~43 arcsec/century" do
      m_sun = 1.989e30
      a = 5.79e10
      e = 0.206
      # per orbit, then scale to century
      delta_phi = GR.orbital_precession(m_sun, a, e)
      period_s = 88.0 * 24 * 3600
      century_s = 100 * 365.25 * 24 * 3600
      orbits_per_century = century_s / period_s
      arcsec = delta_phi * 180 / @pi * 3600 * orbits_per_century
      assert_in_delta arcsec, 43.0, 2.0
    end

    test "positive precession for any mass and orbit" do
      assert GR.orbital_precession(1.989e30, 1.0e11, 0.3) > 0
    end
  end

  describe "light_deflection/2" do
    test "solar limb deflection: ~1.75 arcsec" do
      m_sun = 1.989e30
      r_sun = 6.957e8
      alpha = GR.light_deflection(m_sun, r_sun)
      arcsec = alpha * 180 / @pi * 3600
      assert_in_delta arcsec, 1.75, 0.05
    end

    test "positive deflection" do
      assert GR.light_deflection(1.989e30, 6.957e8) > 0
    end
  end

  # ---------------------------------------------------------------------------
  # Cosmology
  # ---------------------------------------------------------------------------

  describe "critical_density/1" do
    test "positive for any H > 0" do
      h0 = 2.27e-18
      assert GR.critical_density(h0) > 0
    end

    test "ρ_c scales as H²" do
      h1 = 2.27e-18
      h2 = 2 * h1
      rho1 = GR.critical_density(h1)
      rho2 = GR.critical_density(h2)
      assert_in_delta rho2 / rho1, 4.0, @tol
    end

    test "H₀=67 km/s/Mpc → ρ_c ≈ 9.47e-27 kg/m³" do
      h0_si = 67.0 * 1000 / 3.0856e22
      rho_c = GR.critical_density(h0_si)
      assert_in_delta rho_c, 9.47e-27, 1.0e-28
    end
  end

  describe "density_parameter/2" do
    test "Ω = 1 for ρ = ρ_c" do
      h0 = 2.27e-18
      rho_c = GR.critical_density(h0)
      assert_in_delta GR.density_parameter(rho_c, rho_c), 1.0, @tol
    end

    test "Ω > 1 for overcritical density" do
      assert GR.density_parameter(2.0e-26, 1.0e-26) == 2.0
    end
  end

  describe "hubble_from_density/1" do
    test "round-trips with critical_density" do
      h0 = 2.27e-18
      rho_c = GR.critical_density(h0)
      h_back = GR.hubble_from_density(rho_c)
      assert_in_delta h_back, h0, 1.0e-25
    end
  end

  describe "lookback_time/2" do
    test "z=0 → lookback time = 0" do
      h0 = 2.27e-18
      assert_in_delta GR.lookback_time(0.0, h0), 0.0, @tol
    end

    test "higher z → longer lookback time" do
      h0 = 2.27e-18
      t1 = GR.lookback_time(0.5, h0)
      t2 = GR.lookback_time(1.0, h0)
      assert t2 > t1
    end
  end
end

defmodule AstroEquations.Physics.MaterialsTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Materials

  @tol 1.0e-8

  describe "density/2" do
    test "ρ = m/V: 10 kg / 2 m³ = 5 kg/m³" do
      assert_in_delta Materials.density(10, 2), 5.0, @tol
    end

    test "water: 1 kg / 0.001 m³ = 1000 kg/m³" do
      assert_in_delta Materials.density(1.0, 0.001), 1000.0, @tol
    end

    test "raises on zero volume" do
      assert_raise ArgumentError, fn -> Materials.density(10, 0) end
    end
  end

  describe "continuous_density/3" do
    test "uniform density function: dm/dV = constant" do
      rho_0 = 2.5
      mass_fn = fn v -> rho_0 * v end
      assert_in_delta Materials.continuous_density(mass_fn, 2.0), rho_0, 1.0e-6
    end

    test "linear increase in density" do
      mass_fn = fn v -> v * v end
      rho_at_2 = Materials.continuous_density(mass_fn, 2.0)
      assert_in_delta rho_at_2, 4.0, 1.0e-5
    end
  end

  describe "number_density/2" do
    test "n = N/V: 1e23 particles in 0.001 m³ → 1e26 m⁻³" do
      assert_in_delta Materials.number_density(1.0e23, 0.001), 1.0e26, 1.0e18
    end
  end

  describe "normal_stress/2 and normal_strain/2" do
    test "stress = F/A: 1000 N / 0.01 m² = 100 kPa" do
      assert_in_delta Materials.normal_stress(1000, 0.01), 100_000.0, @tol
    end

    test "strain = ΔL/L₀: extension 1 mm / 1 m = 0.001" do
      assert_in_delta Materials.normal_strain(0.001, 1.0), 0.001, @tol
    end
  end

  describe "youngs_modulus/2" do
    test "E = σ/ε: 200 MPa / 0.001 = 200 GPa" do
      assert_in_delta Materials.youngs_modulus(200.0e6, 0.001), 200.0e9, 1.0
    end
  end

  describe "extension/4" do
    test "ΔL = FL₀/(EA): 1 kN, 1 m, 200 GPa, 1e-4 m² → 0.05 mm" do
      delta_l = Materials.extension(1000, 1.0, 200.0e9, 1.0e-4)
      assert_in_delta delta_l, 5.0e-5, 1.0e-8
    end

    test "stiffer material → less extension" do
      dl_steel = Materials.extension(1000, 1.0, 200.0e9, 1.0e-4)
      dl_alu = Materials.extension(1000, 1.0, 70.0e9, 1.0e-4)
      assert dl_alu > dl_steel
    end
  end

  describe "bulk_modulus/3" do
    test "K > 0 for incompressible materials" do
      k = Materials.bulk_modulus(1.0e6, -1.0e-6, 1.0e-3)
      assert k > 0
    end
  end

  describe "elastic_energy_density/2" do
    test "u = ½Eε²: always positive" do
      assert Materials.elastic_energy_density(200.0e9, 0.001) > 0
    end

    test "scales as ε²" do
      u1 = Materials.elastic_energy_density(200.0e9, 0.001)
      u2 = Materials.elastic_energy_density(200.0e9, 0.002)
      assert_in_delta u2 / u1, 4.0, @tol
    end
  end

  describe "linear_expansion/3" do
    test "ΔL = αL₀ΔT: steel 12e-6/K, 1 m, 100 K → 1.2 mm" do
      assert_in_delta Materials.linear_expansion(12.0e-6, 1.0, 100), 1.2e-3, 1.0e-10
    end

    test "proportional to original length" do
      dl1 = Materials.linear_expansion(12.0e-6, 1.0, 100)
      dl2 = Materials.linear_expansion(12.0e-6, 2.0, 100)
      assert_in_delta dl2 / dl1, 2.0, @tol
    end
  end

  describe "heat_conduction/4" do
    test "Q/t = kAΔT/d: copper k=401, A=0.01, ΔT=100, d=0.01 → 40100 W" do
      assert_in_delta Materials.heat_conduction(401, 0.01, 100, 0.01), 40_100.0, 1.0
    end

    test "thicker material → less heat flow" do
      q1 = Materials.heat_conduction(401, 0.01, 100, 0.01)
      q2 = Materials.heat_conduction(401, 0.01, 100, 0.02)
      assert q2 < q1
    end
  end

  describe "reynolds_number/4" do
    test "Re = ρvL/η: air at 30 m/s, 10 cm chord → ~200_000" do
      re = Materials.reynolds_number(1.225, 30, 0.1, 1.81e-5)
      assert_in_delta re, 203_000.0, 5000.0
    end

    test "Re > 4000 → turbulent flow (approximate)" do
      re = Materials.reynolds_number(1000, 1.0, 0.01, 1.0e-3)
      assert re > 4000
    end
  end

  describe "speed_of_sound_solid/2" do
    test "steel: v ≈ 5000 m/s" do
      v = Materials.speed_of_sound_solid(200.0e9, 7900.0)
      assert_in_delta v, 5032.0, 10.0
    end

    test "denser material → slower speed (at same E)" do
      v1 = Materials.speed_of_sound_solid(200.0e9, 7900.0)
      v2 = Materials.speed_of_sound_solid(200.0e9, 15_800.0)
      assert v2 < v1
    end
  end

  describe "speed_of_sound_gas/4" do
    test "air at 20°C: v ≈ 343 m/s" do
      v = Materials.speed_of_sound_gas(1.4, 293.0, 0.029)
      assert_in_delta v, 343.0, 5.0
    end

    test "hotter air → faster sound" do
      v_cold = Materials.speed_of_sound_gas(1.4, 273.0, 0.029)
      v_hot = Materials.speed_of_sound_gas(1.4, 373.0, 0.029)
      assert v_hot > v_cold
    end
  end

  describe "polytropic_pressure/3" do
    test "P = Kρ^Γ: positive for positive inputs" do
      assert Materials.polytropic_pressure(1.0e5, 1000.0, 1.4) > 0
    end

    test "higher density → higher pressure (Γ > 0)" do
      p1 = Materials.polytropic_pressure(1.0e5, 1000.0, 1.4)
      p2 = Materials.polytropic_pressure(1.0e5, 2000.0, 1.4)
      assert p2 > p1
    end
  end
end
