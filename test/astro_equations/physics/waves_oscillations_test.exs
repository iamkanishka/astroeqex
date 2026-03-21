defmodule AstroEquations.Physics.WavesTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Waves

  @pi :math.pi()
  @tol 1.0e-8

  describe "wave_number/1" do
    test "k = 2π/λ: λ = 2 → k = π" do
      assert_in_delta Waves.wave_number(2.0), @pi, @tol
    end

    test "k = 2π/λ: λ = 1 m → k = 2π rad/m" do
      assert_in_delta Waves.wave_number(1.0), 2 * @pi, @tol
    end
  end

  describe "wave_velocity/2" do
    test "v = fλ: 440 Hz × 0.78 m ≈ 343.2 m/s" do
      assert_in_delta Waves.wave_velocity(440, 0.780), 343.2, 0.1
    end

    test "v = fλ: proportional to both frequency and wavelength" do
      v1 = Waves.wave_velocity(100, 1.0)
      v2 = Waves.wave_velocity(200, 1.0)
      assert_in_delta v2 / v1, 2.0, @tol
    end
  end

  describe "frequency/2 and wavelength/2" do
    test "f = v/λ: 340 m/s, 1 m → 340 Hz" do
      assert_in_delta Waves.frequency(340, 1.0), 340.0, @tol
    end

    test "λ = v/f: round trip" do
      v = 340.0
      f = 440.0
      lambda = Waves.wavelength(v, f)
      assert_in_delta Waves.frequency(v, lambda), f, @tol
    end
  end

  describe "period/1" do
    test "T = 1/f: 440 Hz → T ≈ 2.27 ms" do
      assert_in_delta Waves.period(440), 1.0 / 440, @tol
    end
  end

  describe "angular_frequency/1" do
    test "from period: omega = 2π/T" do
      assert_in_delta Waves.angular_frequency(period: 2.0), @pi, @tol
    end

    test "from frequency: omega = 2πf" do
      assert_in_delta Waves.angular_frequency(frequency: 50), 2 * @pi * 50, @tol
    end

    test "raises on missing keyword" do
      assert_raise ArgumentError, fn -> Waves.angular_frequency(amplitude: 1.0) end
    end
  end

  describe "wave_function/6" do
    test "at t=0, x=0, φ=0: y = A sin(0) = 0" do
      assert_in_delta Waves.wave_function(1.0, 1.0, 1.0, 0, 0), 0.0, @tol
    end

    test "bounded by amplitude" do
      a = 2.0

      for t <- [0.1, 0.5, 1.0, 3.0] do
        y = Waves.wave_function(a, 5.0, 3.0, t, 0.5)
        assert abs(y) <= a + @tol
      end
    end
  end

  describe "standing_wave/5" do
    test "node at x = 0 (sin(0) = 0)" do
      assert_in_delta Waves.standing_wave(1.0, @pi, 2 * @pi, 0.0, 0.5), 0.0, @tol
    end
  end

  describe "beat_frequency/2" do
    test "|440 - 444| = 4 Hz beats" do
      assert_in_delta Waves.beat_frequency(440, 444), 4.0, @tol
    end

    test "symmetric: |f₁ - f₂| = |f₂ - f₁|" do
      assert_in_delta Waves.beat_frequency(444, 440), Waves.beat_frequency(440, 444), @tol
    end
  end

  describe "single_slit_minimum/3" do
    test "angle is positive for positive order" do
      assert Waves.single_slit_minimum(1, 500.0e-9, 0.0001) > 0
    end

    test "larger slit → smaller angle" do
      θ_narrow = Waves.single_slit_minimum(1, 500.0e-9, 0.0001)
      θ_wide = Waves.single_slit_minimum(1, 500.0e-9, 0.001)
      assert θ_wide < θ_narrow
    end
  end

  describe "grating_resolving_power/2" do
    test "R = mN: m=1, N=500 → R=500" do
      assert Waves.grating_resolving_power(1, 500) == 500
    end
  end

  describe "intensity_decibels/2" do
    test "threshold intensity (1e-12 W/m²) → 0 dB" do
      assert_in_delta Waves.intensity_decibels(1.0e-12), 0.0, @tol
    end

    test "10× intensity → +10 dB" do
      db1 = Waves.intensity_decibels(1.0e-6)
      db2 = Waves.intensity_decibels(1.0e-5)
      assert_in_delta db2 - db1, 10.0, 1.0e-6
    end
  end

  describe "intensity_from_decibels/2" do
    test "round-trip with intensity_decibels" do
      i = 1.0e-7
      assert_in_delta Waves.intensity_from_decibels(Waves.intensity_decibels(i)), i, i * 1.0e-10
    end
  end

  describe "string_harmonic/3" do
    test "fundamental: n=1, v=340, L=0.68 → 250 Hz" do
      assert_in_delta Waves.string_harmonic(1, 340, 0.68), 250.0, 0.1
    end

    test "second harmonic is double the fundamental" do
      f1 = Waves.string_harmonic(1, 340, 0.68)
      f2 = Waves.string_harmonic(2, 340, 0.68)
      assert_in_delta f2, 2 * f1, @tol
    end
  end

  describe "closed_pipe_harmonic/3" do
    test "only odd harmonics" do
      f1 = Waves.closed_pipe_harmonic(1, 340, 0.25)
      f2 = Waves.closed_pipe_harmonic(2, 340, 0.25)
      assert_in_delta f2 / f1, 3.0, @tol
    end
  end

  describe "doppler/4" do
    test "approaching observer: f > f_source" do
      f_src = 440.0
      f_obs = Waves.doppler(f_src, 20.0, 0.0)
      assert f_obs > f_src
    end

    test "receding observer: f < f_source" do
      f_src = 440.0
      f_obs = Waves.doppler(f_src, -20.0, 0.0)
      assert f_obs < f_src
    end

    test "stationary observer and source: no shift" do
      assert_in_delta Waves.doppler(440.0, 0.0, 0.0), 440.0, @tol
    end
  end

  describe "snells_law/3" do
    test "normal incidence (θ₁=0): no refraction" do
      assert_in_delta Waves.snells_law(1.0, 0.0, 1.5), 0.0, @tol
    end

    test "glass (n=1.5) bends ray toward normal (θ₂ < θ₁)" do
      θ1 = @pi / 4
      θ2 = Waves.snells_law(1.0, θ1, 1.5)
      assert θ2 < θ1
    end
  end

  describe "critical_angle/2" do
    test "glass→air: sin(θ_c) = n₂/n₁ = 1/1.5 → θ_c ≈ 41.8°" do
      θ_c = Waves.critical_angle(1.5, 1.0)
      assert_in_delta θ_c * 180 / @pi, 41.8, 0.1
    end
  end

  describe "malus_law/2" do
    test "at θ=0: I = I₀" do
      assert_in_delta Waves.malus_law(100.0, 0.0), 100.0, @tol
    end

    test "at θ=π/2: I = 0" do
      assert_in_delta Waves.malus_law(100.0, @pi / 2), 0.0, 1.0e-12
    end

    test "at θ=π/4: I = I₀/2" do
      assert_in_delta Waves.malus_law(100.0, @pi / 4), 50.0, @tol
    end
  end

  describe "spherical_wave_intensity/2" do
    test "I ∝ 1/r²: double distance → quarter intensity" do
      i1 = Waves.spherical_wave_intensity(100, 1.0)
      i2 = Waves.spherical_wave_intensity(100, 2.0)
      assert_in_delta i1 / i2, 4.0, @tol
    end
  end
end

defmodule AstroEquations.Physics.OscillationsTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Oscillations

  @pi :math.pi()
  @tol 1.0e-10

  describe "force/2" do
    test "F = -kx: k=10, x=0.5 → -5 N" do
      assert_in_delta Oscillations.force(10, 0.5), -5.0, @tol
    end

    test "restoring: negative for positive displacement" do
      assert Oscillations.force(10, 1.0) < 0
    end

    test "zero displacement gives zero force" do
      assert Oscillations.force(100, 0) == 0.0
    end
  end

  describe "potential_energy/2" do
    test "U = ½kx²: k=10, x=0.5 → 1.25 J" do
      assert_in_delta Oscillations.potential_energy(10, 0.5), 1.25, @tol
    end

    test "always non-negative" do
      assert Oscillations.potential_energy(10, -2.0) >= 0
    end
  end

  describe "angular_frequency/2" do
    test "omega = √(k/m): k=10, m=2.5 → 2 rad/s" do
      assert_in_delta Oscillations.angular_frequency(10, 2.5), 2.0, @tol
    end

    test "stiffer spring → higher frequency" do
      omega1 = Oscillations.angular_frequency(10, 1.0)
      omega2 = Oscillations.angular_frequency(40, 1.0)
      assert omega2 > omega1
    end
  end

  describe "spring_period/2" do
    test "T = 2π √(m/k): m=0.5, k=10 → T ≈ 1.405 s" do
      t = Oscillations.spring_period(0.5, 10)
      assert_in_delta t, 2 * @pi * :math.sqrt(0.5 / 10), @tol
    end
  end

  describe "pendulum_period/2" do
    test "T = 2π √(L/g): L=1 m → T ≈ 2.007 s" do
      assert_in_delta Oscillations.pendulum_period(1.0), 2.0065, 0.001
    end

    test "longer pendulum has longer period" do
      assert Oscillations.pendulum_period(2.0) > Oscillations.pendulum_period(1.0)
    end
  end

  describe "shm_displacement/4" do
    test "x(0) = A at t=0, φ=0" do
      a = 0.05
      assert_in_delta Oscillations.shm_displacement(a, 10, 0, 0), a, @tol
    end

    test "periodic: x(T) = x(0)" do
      a = 0.1
      omega = 3.0
      t_period = 2 * @pi / omega
      x0 = Oscillations.shm_displacement(a, omega, 0)
      xT = Oscillations.shm_displacement(a, omega, t_period)
      assert_in_delta x0, xT, 1.0e-12
    end
  end

  describe "shm_velocity_from_position/3" do
    test "maximum speed at x=0" do
      a = 0.1
      omega = 5.0
      v_max = Oscillations.shm_velocity_from_position(a, omega, 0)
      assert_in_delta v_max, a * omega, @tol
    end

    test "zero speed at x=A" do
      a = 0.1
      assert_in_delta Oscillations.shm_velocity_from_position(a, 5.0, a), 0.0, @tol
    end
  end

  describe "shm_total_energy/2" do
    test "E = ½kA²: k=10, A=0.1 → 0.05 J" do
      assert_in_delta Oscillations.shm_total_energy(10, 0.1), 0.05, @tol
    end
  end

  describe "damping_ratio/3" do
    test "critical damping: ζ = 1 when b = 2√(km)" do
      k = 10.0
      m = 1.0
      b = 2 * :math.sqrt(k * m)
      assert_in_delta Oscillations.damping_ratio(b, k, m), 1.0, @tol
    end

    test "underdamped: ζ < 1" do
      assert Oscillations.damping_ratio(0.5, 10, 1.0) < 1.0
    end
  end

  describe "damped_frequency/2" do
    test "omega_d < omega₀ (damping reduces frequency)" do
      omega0 = 10.0
      ζ = 0.3
      wd = Oscillations.damped_frequency(omega0, ζ)
      assert omegad < omega0
    end

    test "omega_d = omega₀ √(1 - ζ²)" do
      omega0 = 10.0
      ζ = 0.3
      omegad = Oscillations.damped_frequency(omega0, ζ)
      assert_in_delta omegad, omega0 * :math.sqrt(1 - ζ ** 2), @tol
    end
  end

  describe "lc_angular_frequency/2" do
    test "omega = 1/√(LC): L=1 mH, C=10 μF → omega = 10_000 rad/s" do
      assert_in_delta Oscillations.lc_angular_frequency(1.0e-3, 10.0e-6), 10_000.0, 0.1
    end
  end

  describe "beat_frequency/2" do
    test "beat frequency is |f₁ - f₂|" do
      assert_in_delta Oscillations.beat_frequency(440, 444), 4.0, @tol
    end

    test "no beats when frequencies are equal" do
      assert Oscillations.beat_frequency(440, 440) == 0.0
    end
  end
end
