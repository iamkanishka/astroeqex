defmodule AstroEquations.Physics.ThermodynamicsTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Thermodynamics

  @kb 1.380649e-23
  @h 6.62607015e-34
  @c 299_792_458
  @sigma 5.670374419e-8
  @wien_b 2.897771955e-3

  # ---------------------------------------------------------------------------
  # Ideal Gas Law
  # ---------------------------------------------------------------------------

  describe "ideal_gas_law/1" do
    test "solve for pressure" do
      result = Thermodynamics.ideal_gas_law(p: nil, v: 0.001, n: 6.022e23, k_b: @kb, t: 300)
      assert_in_delta result.p, @kb * 6.022e23 * 300 / 0.001, 1.0
    end

    test "solve for volume" do
      result = Thermodynamics.ideal_gas_law(p: 101_325, v: nil, n: 6.022e23, k_b: @kb, t: 273.15)
      assert result.v > 0
    end

    test "solve for temperature" do
      result = Thermodynamics.ideal_gas_law(p: 101_325, v: 0.0224, n: 6.022e23, k_b: @kb, t: nil)
      assert_in_delta result.t, 273.15, 1.0
    end

    test "pV = NkT consistency" do
      n = 1.0e24
      t = 300.0
      v = 0.01
      p_calc = n * @kb * t / v
      result = Thermodynamics.ideal_gas_law(p: nil, v: v, n: n, k_b: @kb, t: t)
      assert_in_delta result.p, p_calc, 1.0
    end
  end

  describe "rms_speed/3" do
    test "nitrogen at 300 K (~515 m/s)" do
      m_n2 = 4.65e-26
      assert_in_delta Thermodynamics.rms_speed(300, m_n2), 515.0, 5.0
    end

    test "higher temperature → higher RMS speed" do
      m = 4.65e-26
      v_lo = Thermodynamics.rms_speed(300, m)
      v_hi = Thermodynamics.rms_speed(1200, m)
      assert v_hi > v_lo
    end

    test "v_rms > v_p > 0" do
      m = 4.65e-26
      t = 300
      assert Thermodynamics.rms_speed(t, m) > Thermodynamics.most_probable_speed(t, m)
    end
  end

  describe "most_probable_speed/3" do
    test "positive for any positive T and m" do
      assert Thermodynamics.most_probable_speed(300, 1.0e-26) > 0
    end
  end

  describe "mean_kinetic_energy/2" do
    test "⟨KE⟩ = 3kT/2 at 300 K" do
      expected = 1.5 * @kb * 300
      assert_in_delta Thermodynamics.mean_kinetic_energy(300), expected, 1.0e-25
    end
  end

  # ---------------------------------------------------------------------------
  # First Law & Work
  # ---------------------------------------------------------------------------

  describe "first_law/2" do
    test "ΔU = Q - W: 500 J in, 200 J work → 300 J stored" do
      assert_in_delta Thermodynamics.first_law(500, 200), 300.0, 1.0e-10
    end

    test "adiabatic: Q = 0, ΔU = -W" do
      assert_in_delta Thermodynamics.first_law(0, 100), -100.0, 1.0e-10
    end
  end

  describe "isobaric_work/2" do
    test "W = PΔV: 101_325 Pa, ΔV = 0.001 m³" do
      assert_in_delta Thermodynamics.isobaric_work(101_325, 0.001), 101.325, 0.001
    end
  end

  describe "isothermal_work/5" do
    test "expansion (V₂ > V₁) does positive work" do
      w = Thermodynamics.isothermal_work(6.022e23, 300, 0.001, 0.002)
      assert w > 0
    end

    test "compression (V₂ < V₁) does negative work" do
      w = Thermodynamics.isothermal_work(6.022e23, 300, 0.002, 0.001)
      assert w < 0
    end
  end

  describe "adiabatic_pressure/4" do
    test "compression increases pressure" do
      p2 = Thermodynamics.adiabatic_pressure(101_325, 1.0, 0.5, 1.4)
      assert p2 > 101_325
    end

    test "expansion decreases pressure" do
      p2 = Thermodynamics.adiabatic_pressure(101_325, 1.0, 2.0, 1.4)
      assert p2 < 101_325
    end
  end

  # ---------------------------------------------------------------------------
  # Entropy & Second Law
  # ---------------------------------------------------------------------------

  describe "entropy/2" do
    test "S = kB ln Ω > 0 for Ω > 1" do
      assert Thermodynamics.entropy(10) > 0
    end

    test "S = 0 when Ω = 1 (only one microstate)" do
      assert_in_delta Thermodynamics.entropy(1), 0.0, 1.0e-30
    end

    test "more microstates → higher entropy" do
      assert Thermodynamics.entropy(100) > Thermodynamics.entropy(10)
    end
  end

  describe "carnot_efficiency/2" do
    test "η = 1 - T_cold/T_hot: 300 K / 500 K → 40%" do
      assert_in_delta Thermodynamics.carnot_efficiency(500, 300), 0.4, 1.0e-10
    end

    test "η = 0 when T_hot = T_cold" do
      assert_in_delta Thermodynamics.carnot_efficiency(300, 300), 0.0, 1.0e-10
    end

    test "efficiency always < 1 for finite temperatures" do
      assert Thermodynamics.carnot_efficiency(1000, 300) < 1.0
    end
  end

  describe "cop_refrigerator/2" do
    test "COP = T_cold/(T_hot - T_cold)" do
      assert_in_delta Thermodynamics.cop_refrigerator(300, 250), 5.0, 1.0e-10
    end

    test "COP > 0" do
      assert Thermodynamics.cop_refrigerator(400, 300) > 0
    end
  end

  # ---------------------------------------------------------------------------
  # Black-body Radiation
  # ---------------------------------------------------------------------------

  describe "photon_energy/2" do
    test "E = hf: f = 1e15 Hz → ~6.63e-19 J" do
      assert_in_delta Thermodynamics.photon_energy(1.0e15), @h * 1.0e15, 1.0e-30
    end
  end

  describe "wiens_displacement/2" do
    test "solar surface (5778 K) peaks near 501 nm" do
      lambda = Thermodynamics.wiens_displacement(5778)
      assert_in_delta lambda, @wien_b / 5778, 1.0e-12
    end

    test "hotter → shorter wavelength" do
      l1 = Thermodynamics.wiens_displacement(5000)
      l2 = Thermodynamics.wiens_displacement(10_000)
      assert l2 < l1
    end
  end

  describe "stefan_boltzmann/2" do
    test "I = σT⁴: T = 5000 K" do
      assert_in_delta Thermodynamics.stefan_boltzmann(5000), @sigma * 5000 ** 4, 1.0
    end

    test "hotter → much more intensity (T⁴ scaling)" do
      i1 = Thermodynamics.stefan_boltzmann(1000)
      i2 = Thermodynamics.stefan_boltzmann(2000)
      assert_in_delta i2 / i1, 16.0, 0.001
    end
  end

  describe "planck_wavelength/5" do
    test "positive for all physical inputs" do
      assert Thermodynamics.planck_wavelength(500.0e-9, 5778) > 0
    end

    test "peaks near Wien displacement" do
      t = 5778
      lambda_wien = @wien_b / t
      b_peak = Thermodynamics.planck_wavelength(lambda_wien, t)
      b_off = Thermodynamics.planck_wavelength(lambda_wien * 2, t)
      assert b_peak > b_off
    end
  end

  describe "newton_cooling/4" do
    test "at t=0, T = T₀" do
      assert_in_delta Thermodynamics.newton_cooling(100, 20, 0.1, 0), 100.0, 1.0e-10
    end

    test "approaches T_env asymptotically" do
      t_long = Thermodynamics.newton_cooling(100, 20, 0.1, 1000)
      assert_in_delta t_long, 20.0, 0.001
    end

    test "monotonically decreasing when T₀ > T_env" do
      t1 = Thermodynamics.newton_cooling(100, 20, 0.1, 5)
      t2 = Thermodynamics.newton_cooling(100, 20, 0.1, 10)
      assert t2 < t1
    end
  end

  # ---------------------------------------------------------------------------
  # Maxwell-Boltzmann
  # ---------------------------------------------------------------------------

  describe "maxwell_boltzmann/4" do
    test "distribution is non-negative" do
      for v <- [100, 300, 500, 1000] do
        assert Thermodynamics.maxwell_boltzmann(v, 4.65e-26, 300) >= 0
      end
    end

    test "peaks near most probable speed" do
      m = 4.65e-26
      t = 300
      v_p = Thermodynamics.most_probable_speed(t, m)
      f_peak = Thermodynamics.maxwell_boltzmann(v_p, m, t)
      f_off = Thermodynamics.maxwell_boltzmann(v_p * 2, m, t)
      assert f_peak > f_off
    end
  end
end
