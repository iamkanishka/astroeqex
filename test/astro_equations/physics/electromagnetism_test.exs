defmodule AstroEquations.Physics.ElectromagnetismTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Electromagnetism, as: EM

  @eps0 8.8541878128e-12
  @mu0 1.25663706212e-6
  @pi :math.pi()
  @tol 1.0e-8

  # ---------------------------------------------------------------------------
  # Maxwell's Equations
  # ---------------------------------------------------------------------------

  describe "gauss_law/2" do
    test "Φ = Q/ε₀: 1 C → ~1.13e11 N·m²/C" do
      assert_in_delta EM.gauss_law(1.0), 1.0 / @eps0, 1.0
    end

    test "zero charge → zero flux" do
      assert EM.gauss_law(0.0) == 0.0
    end

    test "scales linearly with charge" do
      phi1 = EM.gauss_law(1.0e-9)
      phi2 = EM.gauss_law(2.0e-9)
      assert_in_delta phi2 / phi1, 2.0, @tol
    end
  end

  describe "gauss_law_magnetism/0" do
    test "always returns 0 (no magnetic monopoles)" do
      assert EM.gauss_law_magnetism() == 0
    end
  end

  describe "faraday_law/1" do
    test "curl E = -dB/dt: returns negation" do
      assert_in_delta EM.faraday_law(3.0), -3.0, @tol
    end

    test "zero rate → zero curl" do
      assert EM.faraday_law(0.0) == 0.0
    end
  end

  describe "faraday_emf/1" do
    test "ε = -dΦ/dt: 0.05 Wb/s → -0.05 V" do
      assert_in_delta EM.faraday_emf(0.05), -0.05, @tol
    end
  end

  describe "ampere_law/4" do
    test "pure current: ∇×B = μ₀J" do
      j = 1.0e6
      result = EM.ampere_law(j, 0.0)
      assert_in_delta result, @mu0 * j, 1.0e-15
    end

    test "pure displacement current: ∇×B = μ₀ε₀ ∂E/∂t" do
      dE_dt = 1.0e12
      result = EM.ampere_law(0.0, dE_dt)
      assert_in_delta result, @mu0 * @eps0 * dE_dt, 1.0e-20
    end
  end

  # ---------------------------------------------------------------------------
  # Lorentz Force
  # ---------------------------------------------------------------------------

  describe "lorentz_force_point/4" do
    test "electric force only (B=0): F = qE" do
      q = 1.602e-19
      e = {1.0e4, 0.0, 0.0}
      v = {0.0, 0.0, 0.0}
      b = {0.0, 0.0, 0.0}
      {fx, fy, fz} = EM.lorentz_force_point(q, e, v, b)
      assert_in_delta fx, q * 1.0e4, 1.0e-30
      assert_in_delta fy, 0.0, @tol
      assert_in_delta fz, 0.0, @tol
    end

    test "magnetic force only (E=0, v⊥B): F = qv×B" do
      q = 1.0
      e = {0.0, 0.0, 0.0}
      v = {0.0, 1.0, 0.0}
      b = {0.0, 0.0, 1.0}
      {fx, _, _} = EM.lorentz_force_point(q, e, v, b)
      assert_in_delta fx, -1.0, @tol
    end

    test "charge sign flips force direction" do
      e = {1.0e3, 0.0, 0.0}
      {fx_pos, _, _} = EM.lorentz_force_point(1.0, e, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0})
      {fx_neg, _, _} = EM.lorentz_force_point(-1.0, e, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0})
      assert_in_delta fx_pos, -fx_neg, @tol
    end
  end

  describe "cyclotron_radius/4" do
    test "r = mv⊥/(|q|B): positive value" do
      r = EM.cyclotron_radius(9.109e-31, 1.0e6, 1.602e-19, 0.01)
      assert r > 0
    end

    test "larger B → smaller radius" do
      r1 = EM.cyclotron_radius(1.0e-27, 1.0e5, 1.6e-19, 0.1)
      r2 = EM.cyclotron_radius(1.0e-27, 1.0e5, 1.6e-19, 1.0)
      assert r2 < r1
    end
  end

  # ---------------------------------------------------------------------------
  # Electric Fields and Potentials
  # ---------------------------------------------------------------------------

  describe "electric_field_point/3" do
    test "E = q/(4πε₀r²): 1 nC at 1 m → 8.99 V/m" do
      e = EM.electric_field_point(1.0e-9, 1.0)
      assert_in_delta e, 8.988e0, 0.01
    end

    test "inverse square law: r doubles → E quarters" do
      e1 = EM.electric_field_point(1.0e-9, 1.0)
      e2 = EM.electric_field_point(1.0e-9, 2.0)
      assert_in_delta e1 / e2, 4.0, @tol
    end

    test "negative charge → negative field (pointing inward)" do
      assert EM.electric_field_point(-1.0e-9, 1.0) < 0
    end
  end

  describe "electric_potential_point/3" do
    test "V = q/(4πε₀r): 1 nC at 0.1 m → ~89.9 V" do
      v = EM.electric_potential_point(1.0e-9, 0.1)
      assert_in_delta v, 89.88, 0.1
    end

    test "scales as 1/r" do
      v1 = EM.electric_potential_point(1.0e-9, 1.0)
      v2 = EM.electric_potential_point(1.0e-9, 2.0)
      assert_in_delta v1 / v2, 2.0, @tol
    end
  end

  describe "potential_energy/4" do
    test "like charges: positive PE (repulsive)" do
      assert EM.potential_energy(1.0e-9, 1.0e-9, 0.1) > 0
    end

    test "unlike charges: negative PE (attractive)" do
      assert EM.potential_energy(1.0e-9, -1.0e-9, 0.1) < 0
    end

    test "U = kq₁q₂/r: symmetric in charges" do
      u1 = EM.potential_energy(2.0e-9, 3.0e-9, 0.5)
      u2 = EM.potential_energy(3.0e-9, 2.0e-9, 0.5)
      assert_in_delta u1, u2, @tol
    end
  end

  describe "field_energy_density/2" do
    test "u = ½ε₀E²: positive" do
      assert EM.field_energy_density(1000.0) > 0
    end

    test "scales as E²" do
      u1 = EM.field_energy_density(100.0)
      u2 = EM.field_energy_density(200.0)
      assert_in_delta u2 / u1, 4.0, @tol
    end
  end

  # ---------------------------------------------------------------------------
  # Circuit Theory
  # ---------------------------------------------------------------------------

  describe "ohms_law/2" do
    test "V = IR: 2 A × 5 Ω = 10 V" do
      assert_in_delta EM.ohms_law(2, 5), 10.0, @tol
    end
  end

  describe "electrical_power/2" do
    test "P = IV: 2 A × 10 V = 20 W" do
      assert_in_delta EM.electrical_power(2, 10), 20.0, @tol
    end
  end

  describe "electrical_power_from_resistance/2" do
    test "P = I²R: 2 A, 5 Ω → 20 W" do
      assert_in_delta EM.electrical_power_from_resistance(2, 5), 20.0, @tol
    end
  end

  describe "electrical_power_from_voltage/2" do
    test "P = V²/R: 10 V, 5 Ω → 20 W" do
      assert_in_delta EM.electrical_power_from_voltage(10, 5), 20.0, @tol
    end

    test "consistent with P = IV when I = V/R" do
      v = 12.0
      r = 4.0
      i = v / r

      assert_in_delta EM.electrical_power_from_voltage(v, r),
                      EM.electrical_power(i, v),
                      @tol
    end
  end

  describe "series_resistance/1" do
    test "R = Σ Rᵢ: [10, 20, 30] → 60 Ω" do
      assert_in_delta EM.series_resistance([10, 20, 30]), 60.0, @tol
    end

    test "single resistor unchanged" do
      assert_in_delta EM.series_resistance([47.0]), 47.0, @tol
    end
  end

  describe "parallel_resistance/1" do
    test "two equal R → R/2: [10, 10] → 5 Ω" do
      assert_in_delta EM.parallel_resistance([10, 10]), 5.0, @tol
    end

    test "parallel R always less than smallest" do
      r_min = 10.0
      r_parallel = EM.parallel_resistance([r_min, 100.0, 1000.0])
      assert r_parallel < r_min
    end
  end

  describe "lc_resonant_frequency/2" do
    test "omega₀ = 1/√(LC): L=1mH, C=10μF → 10_000 rad/s" do
      assert_in_delta EM.lc_resonant_frequency(1.0e-3, 10.0e-6), 10_000.0, 0.1
    end

    test "higher L or C → lower resonant frequency" do
      omega1 = EM.lc_resonant_frequency(1.0e-3, 1.0e-6)
      omega2 = EM.lc_resonant_frequency(2.0e-3, 1.0e-6)
      assert omega2 < omega1
    end
  end

  describe "rc_time_constant/2" do
    test "τ = RC: R=10kΩ, C=100μF → 1.0 s" do
      assert_in_delta EM.rc_time_constant(10_000, 100.0e-6), 1.0, @tol
    end
  end

  describe "rl_time_constant/2" do
    test "τ = L/R: L=1H, R=100Ω → 0.01 s" do
      assert_in_delta EM.rl_time_constant(1.0, 100), 0.01, @tol
    end
  end

  # ---------------------------------------------------------------------------
  # Capacitors
  # ---------------------------------------------------------------------------

  describe "capacitance/2" do
    test "C = Q/V: 1e-6 C, 10 V → 0.1 μF" do
      assert_in_delta EM.capacitance(1.0e-6, 10), 1.0e-7, @tol
    end
  end

  describe "capacitor_energy/2" do
    test "U = ½CV²: C=1μF, V=10 V → 50 μJ" do
      assert_in_delta EM.capacitor_energy(1.0e-6, 10), 5.0e-5, @tol
    end

    test "scales as V²" do
      u1 = EM.capacitor_energy(1.0e-6, 10.0)
      u2 = EM.capacitor_energy(1.0e-6, 20.0)
      assert_in_delta u2 / u1, 4.0, @tol
    end
  end

  # ---------------------------------------------------------------------------
  # Magnetic Fields
  # ---------------------------------------------------------------------------

  describe "wire_magnetic_field/3" do
    test "B = μ₀I/(2πr): 1 A, 0.1 m → 2 μT" do
      b = EM.wire_magnetic_field(1.0, 0.1)
      assert_in_delta b, 2.0e-6, 1.0e-10
    end

    test "inverse distance law: double r → half B" do
      b1 = EM.wire_magnetic_field(1.0, 0.1)
      b2 = EM.wire_magnetic_field(1.0, 0.2)
      assert_in_delta b1 / b2, 2.0, @tol
    end
  end

  describe "solenoid_field/3" do
    test "B = μ₀nI: n=1000/m, I=2A → 2.51 mT" do
      b = EM.solenoid_field(1000, 2.0)
      assert_in_delta b, @mu0 * 1000 * 2.0, 1.0e-12
    end
  end

  describe "inductor_energy/2" do
    test "U = ½LI²: L=0.1H, I=2A → 0.2 J" do
      assert_in_delta EM.inductor_energy(0.1, 2.0), 0.2, @tol
    end
  end

  describe "magnetic_energy_density/2" do
    test "u = B²/(2μ₀): positive" do
      assert EM.magnetic_energy_density(1.0) > 0
    end

    test "scales as B²" do
      u1 = EM.magnetic_energy_density(1.0)
      u2 = EM.magnetic_energy_density(2.0)
      assert_in_delta u2 / u1, 4.0, @tol
    end
  end

  # ---------------------------------------------------------------------------
  # Materials
  # ---------------------------------------------------------------------------

  describe "electric_susceptibility/1" do
    test "χₑ = εᵣ - 1: εᵣ=2 → χₑ=1" do
      assert_in_delta EM.electric_susceptibility(2.0), 1.0, @tol
    end

    test "vacuum (εᵣ=1): χₑ = 0" do
      assert_in_delta EM.electric_susceptibility(1.0), 0.0, @tol
    end
  end

  describe "relative_permittivity/2" do
    test "εᵣ = ε/ε₀: ε = 2ε₀ → εᵣ = 2" do
      assert_in_delta EM.relative_permittivity(2 * @eps0), 2.0, @tol
    end
  end

  # ---------------------------------------------------------------------------
  # EM Waves
  # ---------------------------------------------------------------------------

  describe "free_space_impedance/2" do
    test "Z₀ = √(μ₀/ε₀) ≈ 376.7 Ω" do
      z0 = EM.free_space_impedance()
      assert_in_delta z0, 376.73, 0.1
    end
  end

  describe "wave_speed/2" do
    test "vacuum: v = 1/√(μ₀ε₀) = c" do
      v = EM.wave_speed(@mu0, @eps0)
      assert_in_delta v, 2.998e8, 1.0e5
    end

    test "medium with εᵣ=4: v = c/2" do
      v = EM.wave_speed(@mu0, 4 * @eps0)
      assert_in_delta v, 2.998e8 / 2, 1.0e5
    end
  end

  describe "skin_depth/3" do
    test "copper at 1 MHz: δ ≈ 66 μm" do
      omega = 2 * @pi * 1.0e6
      delta = EM.skin_depth(@mu0, 5.8e7, omega)
      assert_in_delta delta * 1.0e6, 66.0, 2.0
    end

    test "higher frequency → smaller skin depth" do
      omega1 = 2 * @pi * 1.0e6
      omega2 = 2 * @pi * 4.0e6
      d1 = EM.skin_depth(@mu0, 5.8e7, omega1)
      d2 = EM.skin_depth(@mu0, 5.8e7, omega2)
      assert d2 < d1
    end
  end

  describe "brewsters_angle/2" do
    test "glass (n=1.5): θ_B ≈ 56.3°" do
      theta_b = EM.brewsters_angle(1.0, 1.5)
      assert_in_delta theta_b * 180 / @pi, 56.3, 0.1
    end

    test "symmetric material: 45°" do
      assert_in_delta EM.brewsters_angle(1.0, 1.0) * 180 / @pi, 45.0, @tol
    end
  end
end
