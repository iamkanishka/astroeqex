defmodule AstroEquations.Physics.SpecialRelativityTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.SpecialRelativity, as: SR

  @c 299_792_458
  @tolerance 1.0e-6

  describe "gamma_factor/2" do
    test "at v = 0, γ = 1" do
      assert_in_delta SR.gamma_factor(0), 1.0, @tolerance
    end

    test "at v = 0.5c, γ ≈ 1.1547" do
      assert_in_delta SR.gamma_factor(0.5 * @c), 1.0 / :math.sqrt(0.75), @tolerance
    end

    test "at v = 0.9c, γ ≈ 2.294" do
      assert_in_delta SR.gamma_factor(0.9 * @c), 2.2942, 0.0001
    end

    test "always ≥ 1" do
      for beta <- [0.0, 0.1, 0.5, 0.9, 0.99] do
        assert SR.gamma_factor(beta * @c) >= 1.0
      end
    end
  end

  describe "beta/2" do
    test "at v = 0, β = 0" do
      assert SR.beta(0) == 0.0
    end

    test "at v = c, β = 1" do
      assert_in_delta SR.beta(@c), 1.0, @tolerance
    end

    test "β = v / c" do
      v = 0.6 * @c
      assert_in_delta SR.beta(v), 0.6, @tolerance
    end
  end

  describe "rapidity/2" do
    test "zero velocity → zero rapidity" do
      assert_in_delta SR.rapidity(0), 0.0, @tolerance
    end

    test "round-trips with speed_from_rapidity" do
      v = 0.7 * @c
      assert_in_delta SR.speed_from_rapidity(SR.rapidity(v)), v, 1.0
    end

    test "rapidity is additive for collinear boosts" do
      v1 = 0.6 * @c
      v2 = 0.5 * @c
      phi1 = SR.rapidity(v1)
      phi2 = SR.rapidity(v2)
      v_combined = SR.relative_velocity(-v2, -v1)
      assert_in_delta SR.rapidity(abs(v_combined)), phi1 + phi2, 0.001
    end
  end

  describe "time_dilation/3" do
    test "at v = 0, t = t₀" do
      assert_in_delta SR.time_dilation(1.0, 0), 1.0, @tolerance
    end

    test "moving clock runs slow: t > t₀" do
      t0 = 1.0
      t = SR.time_dilation(t0, 0.9 * @c)
      assert t > t0
    end

    test "t = γ t₀" do
      t0 = 5.0
      v = 0.6 * @c
      assert_in_delta SR.time_dilation(t0, v), SR.gamma_factor(v) * t0, @tolerance
    end
  end

  describe "length_contraction/3" do
    test "at v = 0, L = L₀" do
      assert_in_delta SR.length_contraction(1.0, 0), 1.0, @tolerance
    end

    test "moving objects are shorter: L < L₀" do
      l0 = 1.0
      l = SR.length_contraction(l0, 0.9 * @c)
      assert l < l0
    end

    test "time_dilation × length_contraction = L₀ × t₀ (invariance check)" do
      gamma = SR.gamma_factor(0.6 * @c)
      t = SR.time_dilation(1.0, 0.6 * @c)
      l = SR.length_contraction(1.0, 0.6 * @c)
      assert_in_delta t, gamma, @tolerance
      assert_in_delta l, 1.0 / gamma, @tolerance
    end
  end

  describe "relative_velocity/3" do
    test "velocities add sub-luminally" do
      u_prime = SR.relative_velocity(0.9 * @c, 0.9 * @c)
      assert u_prime < @c
    end

    test "0 + 0 = 0" do
      assert_in_delta SR.relative_velocity(0, 0), 0.0, @tolerance
    end

    test "classical limit: small velocities" do
      u = 100.0
      v = 200.0
      assert_in_delta SR.relative_velocity(u, v), u - v, 0.001
    end
  end

  describe "relativistic_momentum/3" do
    test "p = γmv at v = 0 gives 0" do
      assert SR.relativistic_momentum(1.0, 0) == 0.0
    end

    test "greater than classical momentum for v > 0" do
      m = 1.0
      v = 0.8 * @c
      p_rel = SR.relativistic_momentum(m, v)
      p_cl = m * v
      assert p_rel > p_cl
    end
  end

  describe "rest_energy/2" do
    test "E₀ = mc²: 1 kg → 8.988e16 J" do
      assert_in_delta SR.rest_energy(1.0), 8.9875e16, 1.0e12
    end
  end

  describe "total_energy/3" do
    test "at v = 0, total energy equals rest energy" do
      assert_in_delta SR.total_energy(1.0, 0), SR.rest_energy(1.0), 1.0
    end

    test "total_energy > rest_energy for v > 0" do
      assert SR.total_energy(1.0, 0.5 * @c) > SR.rest_energy(1.0)
    end
  end

  describe "kinetic_energy/3" do
    test "at v = 0, kinetic energy is 0" do
      assert_in_delta SR.kinetic_energy(1.0, 0), 0.0, 1.0
    end

    test "equals total_energy - rest_energy" do
      m = 1.0
      v = 0.6 * @c
      ke = SR.kinetic_energy(m, v)
      delta = SR.total_energy(m, v) - SR.rest_energy(m)
      assert_in_delta ke, delta, 1.0
    end
  end

  describe "energy_momentum_relation/3" do
    test "massless particle: E = pc" do
      p = 1.0e-22
      assert_in_delta SR.energy_momentum_relation(p, 0.0), p * @c, 1.0
    end

    test "at rest (p = 0): E = mc²" do
      m = 1.0e-27
      assert_in_delta SR.energy_momentum_relation(0.0, m), m * @c ** 2, 1.0
    end
  end

  describe "relativistic_doppler/3" do
    test "approaching source (v > 0): frequency blueshifted" do
      f0 = 1.0e14
      f = SR.relativistic_doppler(f0, 0.5 * @c)
      assert f > f0
    end

    test "receding source (v < 0): frequency redshifted" do
      f0 = 1.0e14
      f = SR.relativistic_doppler(f0, -0.5 * @c)
      assert f < f0
    end

    test "v = 0: no Doppler shift" do
      f0 = 5.0e14
      assert_in_delta SR.relativistic_doppler(f0, 0), f0, @tolerance
    end
  end

  describe "spacetime_interval/5" do
    test "purely temporal event is timelike (negative)" do
      assert SR.spacetime_interval(1.0, 0, 0, 0) < 0
    end

    test "purely spatial event is spacelike (positive)" do
      assert SR.spacetime_interval(0, 1.0, 0, 0) > 0
    end

    test "light-like event gives 0" do
      assert_in_delta SR.spacetime_interval(1.0, @c, 0, 0), 0.0, 1.0
    end
  end

  describe "lorentz_boost/6" do
    test "v = 0 leaves coordinates unchanged" do
      {ct_p, x_p, y_p, z_p} = SR.lorentz_boost(1.0, 2.0, 3.0, 4.0, 0)
      assert_in_delta ct_p, 1.0, @tolerance
      assert_in_delta x_p, 2.0, @tolerance
      assert_in_delta y_p, 3.0, @tolerance
      assert_in_delta z_p, 4.0, @tolerance
    end

    test "y and z components are unchanged" do
      {_, _, y_p, z_p} = SR.lorentz_boost(1.0, 0.0, 5.0, 6.0, 0.5 * @c)
      assert_in_delta y_p, 5.0, @tolerance
      assert_in_delta z_p, 6.0, @tolerance
    end
  end

  describe "proper_time/5" do
    test "stationary clock: proper time = coordinate time" do
      tau = SR.proper_time(0, 1.0, fn _ -> 0 end)
      assert_in_delta tau, 1.0, 0.001
    end

    test "moving clock records less proper time" do
      v = 0.9 * @c
      tau = SR.proper_time(0, 1.0, fn _ -> v end)
      assert tau < 1.0
    end

    test "proper time ≈ t₀ / γ for constant velocity" do
      v = 0.6 * @c
      tau = SR.proper_time(0, 1.0, fn _ -> v end)
      expected = 1.0 / SR.gamma_factor(v)
      assert_in_delta tau, expected, 0.001
    end
  end
end
