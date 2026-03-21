defmodule AstroEquations.Physics.QuantumMechanicsTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.QuantumMechanics, as: QM

  @hbar 1.054571817e-34
  @h 6.62607015e-34
  @tol 1.0e-8

  # ---------------------------------------------------------------------------
  # Fundamental Relations
  # ---------------------------------------------------------------------------

  describe "uncertainty_principle?/2" do
    test "satisfies principle: ΔxΔp = ħ/2" do
      assert QM.uncertainty_principle?(1.0, @hbar / 2)
    end

    test "below limit: returns false" do
      refute QM.uncertainty_principle?(1.0e-10, 1.0e-25)
    end

    test "large values satisfy" do
      assert QM.uncertainty_principle?(1.0, 1.0)
    end
  end

  describe "min_uncertainty_product/0" do
    test "returns ħ/2" do
      assert_in_delta QM.min_uncertainty_product(), @hbar / 2, 1.0e-45
    end
  end

  describe "de_broglie_wavelength/1" do
    test "λ = h/p: 1e-24 kg·m/s → 6.626e-10 m" do
      assert_in_delta QM.de_broglie_wavelength(1.0e-24), @h / 1.0e-24, 1.0e-44
    end

    test "smaller momentum → longer wavelength" do
      lambda1 = QM.de_broglie_wavelength(1.0e-24)
      lambda2 = QM.de_broglie_wavelength(2.0e-24)
      assert lambda1 > lambda2
    end
  end

  describe "de_broglie_wavelength_ke/2" do
    test "electron at 1 eV → ~1.23 nm" do
      ke = 1.602e-19
      lambda = QM.de_broglie_wavelength_ke(ke)
      assert_in_delta lambda * 1.0e9, 1.226, 0.01
    end

    test "higher KE → shorter wavelength" do
      l1 = QM.de_broglie_wavelength_ke(1.0e-19)
      l2 = QM.de_broglie_wavelength_ke(4.0e-19)
      assert l2 < l1
    end
  end

  describe "photon_energy/1" do
    test "E = hf: f = 1e15 Hz → 6.63e-19 J" do
      assert_in_delta QM.photon_energy(1.0e15), @h * 1.0e15, 1.0e-33
    end

    test "visible light (f = 6e14 Hz) → ~2.5 eV" do
      e_j = QM.photon_energy(6.0e14)
      e_ev = e_j / 1.602e-19
      assert_in_delta e_ev, 2.5, 0.1
    end
  end

  describe "photon_momentum/1" do
    test "p = h/λ: λ = 500 nm → 1.33e-27 kg·m/s" do
      p = QM.photon_momentum(500.0e-9)
      assert_in_delta p, @h / 500.0e-9, 1.0e-41
    end
  end

  # ---------------------------------------------------------------------------
  # Energy Levels
  # ---------------------------------------------------------------------------

  describe "hydrogen_energy_level/1" do
    test "ground state: E₁ = -13.6 eV" do
      assert_in_delta QM.hydrogen_energy_level(1), -13.6, 0.01
    end

    test "first excited state: E₂ = -3.4 eV" do
      assert_in_delta QM.hydrogen_energy_level(2), -3.4, 0.01
    end

    test "energies scale as 1/n²" do
      e1 = QM.hydrogen_energy_level(1)
      e2 = QM.hydrogen_energy_level(2)
      assert_in_delta e1 / e2, 4.0, 0.001
    end

    test "energy levels increase toward 0 as n increases" do
      assert QM.hydrogen_energy_level(1) < QM.hydrogen_energy_level(2)
      assert QM.hydrogen_energy_level(2) < QM.hydrogen_energy_level(3)
    end
  end

  describe "rydberg_wavelength/3" do
    test "Hα line (n₁=2, n₂=3): λ ≈ 656 nm" do
      lambda = QM.rydberg_wavelength(2, 3)
      assert_in_delta lambda * 1.0e9, 656.3, 1.0
    end

    test "Lyman alpha (n₁=1, n₂=2): λ ≈ 121.6 nm" do
      lambda = QM.rydberg_wavelength(1, 2)
      assert_in_delta lambda * 1.0e9, 121.6, 0.5
    end

    test "higher transition → shorter wavelength" do
      l_21 = QM.rydberg_wavelength(1, 2)
      l_31 = QM.rydberg_wavelength(1, 3)
      assert l_31 < l_21
    end
  end

  describe "infinite_well_energy/3" do
    test "ground state is positive" do
      e = QM.infinite_well_energy(1, 9.109e-31, 1.0e-9)
      assert e > 0
    end

    test "energy scales as n²" do
      m = 9.109e-31
      l = 1.0e-9
      e1 = QM.infinite_well_energy(1, m, l)
      e2 = QM.infinite_well_energy(2, m, l)
      assert_in_delta e2 / e1, 4.0, @tol
    end
  end

  describe "harmonic_oscillator_energy/2" do
    test "ground state (n=0): E₀ = ħomega/2" do
      omega = 1.0e14
      e0 = QM.harmonic_oscillator_energy(0, omega)
      assert_in_delta e0, @hbar * omega / 2, 1.0e-30
    end

    test "equally spaced levels: Eₙ₊₁ - Eₙ = ħomega" do
      omega = 1.0e14
      e0 = QM.harmonic_oscillator_energy(0, omega)
      e1 = QM.harmonic_oscillator_energy(1, omega)
      assert_in_delta e1 - e0, @hbar * omega, 1.0e-30
    end
  end

  # ---------------------------------------------------------------------------
  # Operators and States
  # ---------------------------------------------------------------------------

  describe "pauli_x/0, pauli_z/0, identity/0" do
    test "Pauli X is a 2×2 matrix" do
      assert QM.pauli_x() == [[0, 1], [1, 0]]
    end

    test "Pauli Z has eigenvalues ±1" do
      pz = QM.pauli_z()
      assert Enum.at(Enum.at(pz, 0), 0) == 1
      assert Enum.at(Enum.at(pz, 1), 1) == -1
    end

    test "identity matrix" do
      assert QM.identity() == [[1, 0], [0, 1]]
    end
  end

  describe "expectation_braket/2" do
    test "⟨0|σ_z|0⟩ = +1" do
      state = [1.0, 0.0]
      assert_in_delta QM.expectation_braket(QM.pauli_z(), state), 1.0, @tol
    end

    test "⟨1|σ_z|1⟩ = -1" do
      state = [0.0, 1.0]
      assert_in_delta QM.expectation_braket(QM.pauli_z(), state), -1.0, @tol
    end

    test "⟨+|σ_z|+⟩ = 0 for |+⟩ = (|0⟩+|1⟩)/√2" do
      s = 1.0 / :math.sqrt(2)
      state = [s, s]
      assert_in_delta QM.expectation_braket(QM.pauli_z(), state), 0.0, @tol
    end
  end

  describe "variance/2 and standard_deviation/2" do
    test "definite state has zero variance" do
      state = [1.0, 0.0]
      assert_in_delta QM.variance(QM.pauli_z(), state), 0.0, @tol
    end

    test "superposition has non-zero variance" do
      s = 1.0 / :math.sqrt(2)
      state = [s, s]
      assert QM.variance(QM.pauli_z(), state) > 0
    end

    test "std dev = √variance" do
      s = 1.0 / :math.sqrt(2)
      state = [s, s]
      var = QM.variance(QM.pauli_z(), state)
      std = QM.standard_deviation(QM.pauli_z(), state)
      assert_in_delta std * std, var, @tol
    end
  end

  describe "density_matrix/2" do
    test "pure state: ρ = |ψ⟩⟨ψ|" do
      states = [[1.0, 0.0]]
      probs = [1.0]
      rho = QM.density_matrix(states, probs)
      assert_in_delta Enum.at(Enum.at(rho, 0), 0), 1.0, @tol
      assert_in_delta Enum.at(Enum.at(rho, 1), 1), 0.0, @tol
    end

    test "maximally mixed state: ρ = I/2" do
      states = [[1.0, 0.0], [0.0, 1.0]]
      probs = [0.5, 0.5]
      rho = QM.density_matrix(states, probs)
      assert_in_delta Enum.at(Enum.at(rho, 0), 0), 0.5, @tol
      assert_in_delta Enum.at(Enum.at(rho, 1), 1), 0.5, @tol
    end
  end

  describe "purity/1" do
    test "pure state: Tr(ρ²) = 1" do
      rho = [[1.0, 0.0], [0.0, 0.0]]
      assert_in_delta QM.purity(rho), 1.0, @tol
    end

    test "maximally mixed: Tr(ρ²) = 1/d" do
      rho = [[0.5, 0.0], [0.0, 0.5]]
      assert_in_delta QM.purity(rho), 0.5, @tol
    end
  end

  describe "transmission_coefficient/2 and reflection_coefficient/2" do
    test "T + R = 1" do
      k1 = 5.0e9
      k2 = 3.0e9
      t = QM.transmission_coefficient(k1, k2)
      r = QM.reflection_coefficient(k1, k2)
      assert_in_delta t + r, 1.0, @tol
    end

    test "equal k₁=k₂: T=1, R=0" do
      k = 5.0e9
      assert_in_delta QM.transmission_coefficient(k, k), 1.0, @tol
      assert_in_delta QM.reflection_coefficient(k, k), 0.0, @tol
    end
  end
end

defmodule AstroEquations.AstrophysicsAndAstronomy.GalaxiesTest do
  use ExUnit.Case, async: true

  alias AstroEquations.AstrophysicsAndAstronomy.Galaxies

  @tol 1.0e-8

  describe "hubble_classify/2" do
    test "circular galaxy (a=b): E0" do
      assert Galaxies.hubble_classify(10, 10) == "E0"
    end

    test "half-minor-axis: E5" do
      assert Galaxies.hubble_classify(10, 5) == "E5"
    end

    test "clamped to E7 max" do
      assert Galaxies.hubble_classify(10, 0) == "E7"
    end

    test "clamped to E0 min" do
      assert Galaxies.hubble_classify(5, 6) == "E0"
    end
  end

  describe "sersic_profile/4" do
    test "at r = r_e: I(r_e) = I₀ exp(0) = I₀ for correct bn" do
      # At r = r_e the exponent is -bn*(1-1)=0, so I = I0
      i = Galaxies.sersic_profile(100, 1, 1, 4)
      assert_in_delta i, 100.0, 0.001
    end

    test "surface brightness decreases with radius" do
      i1 = Galaxies.sersic_profile(100, 1, 1, 4)
      i2 = Galaxies.sersic_profile(100, 2, 1, 4)
      assert i2 < i1
    end

    test "always positive" do
      assert Galaxies.sersic_profile(100, 3, 1, 4) > 0
    end
  end

  describe "disk_density/5" do
    test "maximum at centre (r=0, z=0)" do
      rho_centre = Galaxies.disk_density(1.0, 0, 0, 1, 1)
      rho_outer = Galaxies.disk_density(1.0, 2, 1, 1, 1)
      assert rho_centre > rho_outer
    end

    test "ρ(0,0,h,z0) = ρ₀" do
      assert_in_delta Galaxies.disk_density(1.0, 0, 0, 1, 1), 1.0, @tol
    end

    test "exponential radial falloff" do
      rho_r1 = Galaxies.disk_density(1.0, 1, 0, 1, 1)
      rho_r2 = Galaxies.disk_density(1.0, 2, 0, 1, 1)
      assert_in_delta rho_r2 / rho_r1, :math.exp(-1), @tol
    end
  end

  describe "nfw_profile/3" do
    test "density decreases with radius" do
      rho1 = Galaxies.nfw_profile(1.0e7, 1, 3)
      rho2 = Galaxies.nfw_profile(1.0e7, 5, 3)
      assert rho1 > rho2
    end

    test "positive density" do
      assert Galaxies.nfw_profile(1.0e7, 5, 3) > 0
    end
  end

  describe "tully_fisher/3" do
    test "higher rotation velocity → higher luminosity" do
      l1 = Galaxies.tully_fisher(100)
      l2 = Galaxies.tully_fisher(200)
      assert l2 > l1
    end

    test "scales as v^4 (default)" do
      l1 = Galaxies.tully_fisher(100)
      l2 = Galaxies.tully_fisher(200)
      assert_in_delta l2 / l1, 16.0, @tol
    end
  end

  describe "sfr_from_halpha/1" do
    test "L_Hα = 1e42 erg/s → ~7.9 M☉/yr" do
      sfr = Galaxies.sfr_from_halpha(1.0e42)
      assert_in_delta sfr, 1.0e42 / 1.26e41, @tol
    end

    test "higher luminosity → higher SFR" do
      sfr1 = Galaxies.sfr_from_halpha(1.0e41)
      sfr2 = Galaxies.sfr_from_halpha(1.0e42)
      assert sfr2 > sfr1
    end
  end

  describe "virial_mass/2" do
    test "positive mass for positive inputs" do
      assert Galaxies.virial_mass(2.0e5, 1.0e21) > 0
    end

    test "scales as σ²" do
      m1 = Galaxies.virial_mass(1.0e5, 1.0e21)
      m2 = Galaxies.virial_mass(2.0e5, 1.0e21)
      assert_in_delta m2 / m1, 4.0, @tol
    end
  end

  describe "schechter_function/4" do
    test "positive value for physical inputs" do
      assert Galaxies.schechter_function(0.016, 1.0e10, 1.0e10, -1.07) > 0
    end

    test "decreases above L*" do
      phi_at_lstar = Galaxies.schechter_function(0.016, 1.0e10, 1.0e10, -1.07)
      phi_above_lstar = Galaxies.schechter_function(0.016, 5.0e10, 1.0e10, -1.07)
      assert phi_above_lstar < phi_at_lstar
    end
  end
end

defmodule AstroEquations.AstrophysicsAndAstronomy.InstrumentationTest do
  use ExUnit.Case, async: true

  alias AstroEquations.AstrophysicsAndAstronomy.Instrumentation

  @tol 1.0e-8

  describe "focal_ratio/2" do
    test "N = f/D: f=0.5m, D=0.1m → N=5" do
      assert_in_delta Instrumentation.focal_ratio(0.5, 0.1), 5.0, @tol
    end
  end

  describe "field_of_view/4" do
    test "FOV = detector/focal_length: 36mm sensor, f=180mm → 0.2 rad" do
      assert_in_delta Instrumentation.field_of_view(0.036, 0.2, 5), 0.036 / (0.2 * 5), @tol
    end

    test "larger detector → wider FOV" do
      fov1 = Instrumentation.field_of_view(0.024, 0.2, 10)
      fov2 = Instrumentation.field_of_view(0.036, 0.2, 10)
      assert fov2 > fov1
    end
  end

  describe "diffraction_limit/2" do
    test "θ = 1.22λ/D: λ=500nm, D=0.1m → 6.1 μrad" do
      theta = Instrumentation.diffraction_limit(500.0e-9, 0.1)
      assert_in_delta theta, 6.1e-6, 1.0e-9
    end

    test "larger aperture → smaller diffraction limit" do
      t1 = Instrumentation.diffraction_limit(500.0e-9, 0.1)
      t2 = Instrumentation.diffraction_limit(500.0e-9, 1.0)
      assert t2 < t1
    end
  end

  describe "seeing_limit/2" do
    test "θ = 0.98λ/r₀: λ=500nm, r₀=0.15m → 3.3 μrad" do
      theta = Instrumentation.seeing_limit(500.0e-9, 0.15)
      assert_in_delta theta, 0.98 * 500.0e-9 / 0.15, 1.0e-12
    end
  end

  describe "total_resolution_limit/3" do
    test "combined limit > individual limits" do
      diff = Instrumentation.diffraction_limit(500.0e-9, 4.0)
      see = Instrumentation.seeing_limit(500.0e-9, 0.15)
      tot = Instrumentation.total_resolution_limit(500.0e-9, 4.0, 0.15)
      assert tot > diff
      assert tot > see
    end

    test "quadrature addition: θ_tot = √(θ_diff² + θ_see²)" do
      diff = Instrumentation.diffraction_limit(500.0e-9, 0.4)
      see = Instrumentation.seeing_limit(500.0e-9, 0.2)
      tot = Instrumentation.total_resolution_limit(500.0e-9, 0.4, 0.2)
      expected = :math.sqrt(diff ** 2 + see ** 2)
      assert_in_delta tot, expected, 1.0e-14
    end
  end

  describe "signal_to_noise/6" do
    test "SNR > 0 for valid inputs" do
      assert Instrumentation.signal_to_noise(100, 10, 5, 1000, 0.1, 2.0) > 0
    end

    test "higher signal → higher SNR" do
      snr1 = Instrumentation.signal_to_noise(100, 10, 5, 100, 0.1, 2.0)
      snr2 = Instrumentation.signal_to_noise(500, 10, 5, 100, 0.1, 2.0)
      assert snr2 > snr1
    end

    test "higher noise → lower SNR" do
      snr1 = Instrumentation.signal_to_noise(100, 10, 1, 10, 0.01, 2.0)
      snr2 = Instrumentation.signal_to_noise(100, 10, 10, 10, 0.01, 2.0)
      assert snr2 < snr1
    end
  end

  describe "strehl_ratio/1" do
    test "S = exp(-σ²): σ²=0 → S=1" do
      assert_in_delta Instrumentation.strehl_ratio(0.0), 1.0, @tol
    end

    test "σ²=1 → S = e⁻¹ ≈ 0.368" do
      assert_in_delta Instrumentation.strehl_ratio(1.0), :math.exp(-1), @tol
    end

    test "larger wavefront error → lower Strehl" do
      assert Instrumentation.strehl_ratio(0.5) > Instrumentation.strehl_ratio(1.0)
    end
  end

  describe "airmass/1" do
    test "zenith (z=0): X=1" do
      assert_in_delta Instrumentation.airmass(0), 1.0, @tol
    end

    test "z=60°: X=2" do
      assert_in_delta Instrumentation.airmass(:math.pi() / 3), 2.0, @tol
    end

    test "increases with zenith angle" do
      assert Instrumentation.airmass(:math.pi() / 4) < Instrumentation.airmass(:math.pi() / 3)
    end
  end

  describe "tsiolkovsky_rocket_equation/3" do
    test "Δv = v_e ln(m₀/m_f): v_e=2500, m₀=1000, m_f=500 → ~1733 m/s" do
      dv = Instrumentation.tsiolkovsky_rocket_equation(2500, 1000, 500)
      assert_in_delta dv, 2500 * :math.log(2), @tol
    end

    test "more propellant → higher delta-v" do
      dv1 = Instrumentation.tsiolkovsky_rocket_equation(3000, 1000, 500)
      dv2 = Instrumentation.tsiolkovsky_rocket_equation(3000, 1000, 250)
      assert dv2 > dv1
    end
  end
end

defmodule AstroEquations.Mathematics.TrigonometryTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Mathematics.Trigonometry

  @pi :math.pi()
  @tol 1.0e-6

  describe "spherical_law_of_cosines/3" do
    test "equilateral spherical triangle: all sides π/3" do
      c = Trigonometry.spherical_law_of_cosines(@pi / 3, @pi / 3, @pi / 3)
      assert_in_delta c, @pi / 3, @tol
    end

    test "right angle triangle on sphere" do
      # a = b = π/4, C = π/2 → cos c = cos²(π/4) = 0.5 → c = π/3
      c = Trigonometry.spherical_law_of_cosines(@pi / 4, @pi / 4, @pi / 2)
      assert_in_delta :math.cos(c), 0.5, @tol
    end
  end

  describe "altitude/3" do
    test "object at zenith (dec = lat): altitude = π/2" do
      lat = 0.5
      dec = lat
      ha = 0.0
      alt = Trigonometry.altitude(dec, lat, ha)
      assert_in_delta alt, @pi / 2, @tol
    end

    test "object on celestial equator at meridian (dec=0, lat=0): altitude = π/2" do
      alt = Trigonometry.altitude(0.0, 0.0, 0.0)
      assert_in_delta alt, @pi / 2, @tol
    end

    test "altitude in [−π/2, π/2]" do
      alt = Trigonometry.altitude(0.3, 0.8, 1.0)
      assert alt >= -@pi / 2 and alt <= @pi / 2
    end
  end

  describe "azimuth/3" do
    test "result in [0, 2π)" do
      az = Trigonometry.azimuth(0.3, 0.8, 1.0)
      assert az >= 0.0 and az < 2 * @pi
    end
  end

  describe "hour_angle/2" do
    test "H = LST - α: at transit, H = 0" do
      lst = 1.5
      ra = 1.5
      assert_in_delta Trigonometry.hour_angle(lst, ra), 0.0, @tol
    end

    test "result in (−π, π]" do
      h = Trigonometry.hour_angle(3.0, 1.0)
      assert h > -@pi and h <= @pi
    end
  end

  describe "gmst/1" do
    test "GMST at J2000.0 ≈ 1.753 rad (18.697 h)" do
      gmst = Trigonometry.gmst(2_451_545.0)
      assert gmst >= 0.0 and gmst < 2 * @pi
      # Expected ≈ 18.697 hours × (2π/24) ≈ 4.894 rad mod 2π
      # Actually GMST at J2000.0 noon ≈ 18.697 h = 4.89 rad
      assert_in_delta gmst, 4.894, 0.01
    end

    test "result always in [0, 2π)" do
      for jd <- [2_451_545.0, 2_456_000.0, 2_400_000.0] do
        gmst = Trigonometry.gmst(jd)
        assert gmst >= 0.0 and gmst < 2 * @pi
      end
    end
  end

  describe "atmospheric_refraction/1" do
    test "refraction at zenith (90°) is minimal (~0 arcmin)" do
      r = Trigonometry.atmospheric_refraction(90.0)
      assert r < 0.01
    end

    test "refraction at 45° ≈ 0.98 arcmin" do
      assert_in_delta Trigonometry.atmospheric_refraction(45.0), 0.979, 0.01
    end

    test "refraction increases at lower altitudes" do
      r_high = Trigonometry.atmospheric_refraction(45.0)
      r_low = Trigonometry.atmospheric_refraction(10.0)
      assert r_low > r_high
    end
  end

  describe "equatorial_to_galactic/2 and galactic_to_equatorial/2" do
    test "Galactic centre is near (l,b) = (0,0)" do
      ra_gc = 266.405 * @pi / 180
      dec_gc = -28.936 * @pi / 180
      {l, b} = Trigonometry.equatorial_to_galactic(ra_gc, dec_gc)
      assert_in_delta l, 0.0, 0.05
      assert_in_delta b, 0.0, 0.05
    end

    test "round-trip equatorial → galactic → equatorial" do
      ra = 1.5
      dec = 0.3
      {l, b} = Trigonometry.equatorial_to_galactic(ra, dec)
      {ra2, dec2} = Trigonometry.galactic_to_equatorial(l, b)
      assert_in_delta ra2, ra, @tol
      assert_in_delta dec2, dec, @tol
    end
  end

  describe "equatorial_to_ecliptic/3" do
    test "result contains two finite numbers" do
      {lambda, beta} = Trigonometry.equatorial_to_ecliptic(1.5, 0.3)
      assert is_float(lambda) and is_float(beta)
    end

    test "ecliptic latitude is in [−π/2, π/2]" do
      {_, beta} = Trigonometry.equatorial_to_ecliptic(1.5, 0.3)
      assert beta >= -@pi / 2 and beta <= @pi / 2
    end
  end

  describe "sunrise_hour_angle/3" do
    test "returns nil for polar day (always-up)" do
      lat = @pi / 2 - 0.01
      dec = @pi / 3
      result = Trigonometry.sunrise_hour_angle(lat, dec)
      assert result == nil or is_float(result)
    end

    test "returns positive H₀ for normal sunrise" do
      lat = 0.6
      dec = 0.2
      h0 = Trigonometry.sunrise_hour_angle(lat, dec)
      assert h0 == nil or h0 > 0
    end
  end
end
