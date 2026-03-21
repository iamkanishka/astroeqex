defmodule AstroEquations.Physics.ForcesTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Forces

  @tol 1.0e-8
  @g 9.81

  describe "newtons_second_law/2" do
    test "F = ma: 70 kg × 9.81 = 686.7 N" do
      assert_in_delta Forces.newtons_second_law(70, @g), 686.7, 0.01
    end

    test "zero mass gives zero force" do
      assert Forces.newtons_second_law(0, 100) == 0.0
    end
  end

  describe "weight/2" do
    test "70 kg at g = 9.81 → 686.7 N" do
      assert_in_delta Forces.weight(70), 686.7, 0.01
    end

    test "custom g" do
      assert_in_delta Forces.weight(70, 1.62), 113.4, 0.1
    end
  end

  describe "gravitational_force/3" do
    test "Earth-Moon force ≈ 2e20 N" do
      f = Forces.gravitational_force(5.972e24, 7.348e22, 3.844e8)
      assert_in_delta f, 1.982e20, 1.0e18
    end

    test "symmetric" do
      f1 = Forces.gravitational_force(1.0e20, 1.0e18, 1.0e9)
      f2 = Forces.gravitational_force(1.0e18, 1.0e20, 1.0e9)
      assert_in_delta f1, f2, 1.0
    end
  end

  describe "buoyancy/2" do
    test "5 kg displaced → 49.05 N" do
      assert_in_delta Forces.buoyancy(5), 49.05, 0.001
    end
  end

  describe "buoyancy_from_density/3" do
    test "1 L of water (1 kg) → 9.81 N" do
      assert_in_delta Forces.buoyancy_from_density(1000, 0.001), 9.81, 0.001
    end
  end

  describe "kinetic_friction/2" do
    test "μ=0.3, N=100 → 30 N" do
      assert_in_delta Forces.kinetic_friction(0.3, 100), 30.0, @tol
    end
  end

  describe "static_friction/2" do
    test "μ=0.5, N=100 → 50 N" do
      assert_in_delta Forces.static_friction(0.5, 100), 50.0, @tol
    end

    test "static > kinetic for same normal force" do
      n = 100.0
      assert Forces.static_friction(0.5, n) > Forces.kinetic_friction(0.3, n)
    end
  end

  describe "spring_force/2" do
    test "F = -kx: k=10, x=2 → -20 N" do
      assert_in_delta Forces.spring_force(10, 2), -20.0, @tol
    end

    test "restoring: negative for positive x" do
      assert Forces.spring_force(5, 1) < 0
    end
  end

  describe "centripetal_force/3" do
    test "F = mv²/r: m=1000, v=30, r=50 → 18000 N" do
      assert_in_delta Forces.centripetal_force(1000, 30, 50), 18_000.0, 0.1
    end
  end

  describe "quadratic_drag/4" do
    test "positive drag force" do
      assert Forces.quadratic_drag(0.47, 1.225, 0.04, 30) > 0
    end

    test "scales as v²" do
      d1 = Forces.quadratic_drag(0.47, 1.225, 0.04, 10)
      d2 = Forces.quadratic_drag(0.47, 1.225, 0.04, 20)
      assert_in_delta d2 / d1, 4.0, @tol
    end
  end

  describe "terminal_velocity/5" do
    test "positive for a falling object" do
      assert Forces.terminal_velocity(75, 0.47, 1.225, 0.6) > 0
    end

    test "heavier object falls faster" do
      v1 = Forces.terminal_velocity(50, 0.47, 1.225, 0.6)
      v2 = Forces.terminal_velocity(100, 0.47, 1.225, 0.6)
      assert v2 > v1
    end
  end

  describe "hydrostatic_pressure/3" do
    test "P = ρgh: 10 m depth in water → 98_100 Pa" do
      assert_in_delta Forces.hydrostatic_pressure(1000, 10), 98_100.0, 1.0
    end
  end

  describe "normal_force_incline/3" do
    test "flat surface (θ=0): N = mg" do
      assert_in_delta Forces.normal_force_incline(10, 0), 10 * @g, 0.001
    end

    test "vertical surface (θ=π/2): N = 0" do
      assert_in_delta Forces.normal_force_incline(10, :math.pi() / 2), 0.0, 1.0e-6
    end
  end

  describe "incline_force_parallel/3" do
    test "flat surface: F_// = 0" do
      assert_in_delta Forces.incline_force_parallel(10, 0), 0.0, @tol
    end

    test "vertical: F_// = mg" do
      assert_in_delta Forces.incline_force_parallel(10, :math.pi() / 2), 10 * @g, 1.0e-6
    end
  end

  describe "impulse/2" do
    test "J = Ft: 50 N × 3 s = 150 N·s" do
      assert_in_delta Forces.impulse(50, 3), 150.0, @tol
    end
  end

  describe "momentum/2" do
    test "p = mv: 10 kg × 5 m/s = 50 kg·m/s" do
      assert_in_delta Forces.momentum(10, 5), 50.0, @tol
    end
  end
end

defmodule AstroEquations.Physics.EnergyTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Energy

  @tol 1.0e-8
  @c 2.99792458e8

  describe "work/2" do
    test "parallel force and displacement: W = Fd" do
      assert_in_delta Energy.work([10, 0, 0], [5, 0, 0]), 50.0, @tol
    end

    test "perpendicular: W = 0" do
      assert_in_delta Energy.work([0, 10, 0], [5, 0, 0]), 0.0, @tol
    end

    test "3-D dot product" do
      assert_in_delta Energy.work([1, 2, 3], [4, 5, 6]), 32.0, @tol
    end
  end

  describe "work_angle/3" do
    test "θ=0: W = Fd" do
      assert_in_delta Energy.work_angle(10, 5, 0), 50.0, @tol
    end

    test "θ=π/2: W = 0" do
      assert_in_delta Energy.work_angle(10, 5, :math.pi() / 2), 0.0, 1.0e-14
    end
  end

  describe "kinetic_energy/2" do
    test "KE = ½mv²: m=4, v=5 → 50 J" do
      assert_in_delta Energy.kinetic_energy(4, 5), 50.0, @tol
    end

    test "zero velocity → zero KE" do
      assert Energy.kinetic_energy(100, 0) == 0.0
    end
  end

  describe "kinetic_energy_from_momentum/2" do
    test "KE = p²/(2m): p=10, m=2 → 25 J" do
      assert_in_delta Energy.kinetic_energy_from_momentum(10, 2), 25.0, @tol
    end

    test "consistent with KE = ½mv²" do
      m = 3.0
      v = 4.0
      p = m * v

      assert_in_delta Energy.kinetic_energy_from_momentum(p, m),
                      Energy.kinetic_energy(m, v),
                      @tol
    end
  end

  describe "rotational_kinetic_energy/2" do
    test "KE_rot = ½Iomega²: I=2, omega=3 → 9 J" do
      assert_in_delta Energy.rotational_kinetic_energy(2, 3), 9.0, @tol
    end
  end

  describe "gravitational_pe/3" do
    test "PE = mgh: m=10, g=9.81, h=5 → 490.5 J" do
      assert_in_delta Energy.gravitational_pe(10, 5), 490.5, 0.001
    end
  end

  describe "gravitational_pe_general/3" do
    test "U = -Gm₁m₂/r: always negative" do
      u = Energy.gravitational_pe_general(5.972e24, 7.348e22, 3.844e8)
      assert u < 0
    end
  end

  describe "spring_pe/2" do
    test "U = ½kx²: k=100, x=0.1 → 0.5 J" do
      assert_in_delta Energy.spring_pe(100, 0.1), 0.5, @tol
    end

    test "always non-negative" do
      assert Energy.spring_pe(50, -0.2) >= 0
    end
  end

  describe "power_from_work/2" do
    test "P = W/t: 100 J / 5 s = 20 W" do
      assert_in_delta Energy.power_from_work(100, 5), 20.0, @tol
    end
  end

  describe "mechanical_power/3" do
    test "P = Fv: F=50, v=10 → 500 W" do
      assert_in_delta Energy.mechanical_power(50, 10), 500.0, @tol
    end

    test "P = Fv cos θ: perpendicular gives 0" do
      assert_in_delta Energy.mechanical_power(50, 10, :math.pi() / 2), 0.0, 1.0e-12
    end
  end

  describe "work_energy_theorem/2" do
    test "W_net = ΔKE: KE increases with positive work" do
      ke_initial = 100.0
      w_net = 50.0
      assert_in_delta Energy.work_energy_theorem(ke_initial, w_net), 150.0, @tol
    end
  end

  describe "escape_velocity/2" do
    test "Earth escape velocity ≈ 11.2 km/s" do
      v = Energy.escape_velocity(5.972e24, 6.371e6)
      assert_in_delta v / 1000, 11.2, 0.1
    end

    test "v_esc = √(2GM/r)" do
      g = 6.67430e-11
      m = 5.972e24
      r = 6.371e6
      expected = :math.sqrt(2 * g * m / r)
      assert_in_delta Energy.escape_velocity(m, r), expected, 0.01
    end
  end

  describe "rest_energy/2" do
    test "E₀ = mc²: 1 kg → 8.988×10¹⁶ J" do
      assert_in_delta Energy.rest_energy(1.0), @c ** 2, 1.0e8
    end
  end

  describe "binding_energy/2" do
    test "ΔE = Δm c²: positive for positive mass defect" do
      assert Energy.binding_energy(1.0e-28) > 0
    end
  end

  describe "shm_total_energy/2" do
    test "E = ½kA²: k=50, A=0.2 → 1 J" do
      assert_in_delta Energy.shm_total_energy(50, 0.2), 1.0, @tol
    end
  end

  describe "shm_kinetic_energy/4" do
    test "maximum at x=0: KE = ½momega²A²" do
      m = 1.0
      omega = 10.0
      a = 0.1
      ke_max = Energy.shm_kinetic_energy(m, omega, a, 0)
      assert_in_delta ke_max, 0.5 * m * omega ** 2 * a ** 2, @tol
    end

    test "zero at x=A (all potential)" do
      assert_in_delta Energy.shm_kinetic_energy(1.0, 10.0, 0.1, 0.1), 0.0, @tol
    end
  end
end
