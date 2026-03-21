defmodule AstroEquations.Physics.MotionTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Physics.Motion

  # ---------------------------------------------------------------------------
  # Kinematics
  # ---------------------------------------------------------------------------

  describe "velocity/2" do
    test "10 m in 2 s → 5 m/s" do
      assert_in_delta Motion.velocity(10, 2), 5.0, 1.0e-10
    end

    test "negative displacement gives negative velocity" do
      assert Motion.velocity(-10, 2) < 0
    end
  end

  describe "acceleration/2" do
    test "20 m/s change over 5 s → 4 m/s²" do
      assert_in_delta Motion.acceleration(20, 5), 4.0, 1.0e-10
    end
  end

  describe "final_velocity/3" do
    test "v = u + at → 0 + 9.81×2 = 19.62 m/s" do
      assert_in_delta Motion.final_velocity(0, 9.81, 2), 19.62, 1.0e-10
    end

    test "deceleration to rest" do
      assert_in_delta Motion.final_velocity(20, -4, 5), 0.0, 1.0e-10
    end
  end

  describe "displacement_suvat/3" do
    test "s = ut + ½at²: free fall 3 s" do
      assert_in_delta Motion.displacement_suvat(0, 3, 9.81), 44.145, 0.001
    end

    test "uniform motion (a = 0)" do
      assert_in_delta Motion.displacement_suvat(10, 5, 0), 50.0, 1.0e-10
    end
  end

  describe "velocity_squared/3" do
    test "v² = u² + 2as — freefall 5 m from rest" do
      v = Motion.velocity_squared(0, 9.81, 5)
      assert_in_delta v, :math.sqrt(2 * 9.81 * 5), 1.0e-6
    end
  end

  describe "displacement_average/3" do
    test "s = (u+v)t/2" do
      assert_in_delta Motion.displacement_average(0, 20, 4), 40.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Newton's Laws & Momentum
  # ---------------------------------------------------------------------------

  describe "newtons_second_law/2" do
    test "F = ma: 5 kg × 2 m/s² = 10 N" do
      assert_in_delta Motion.newtons_second_law(5, 2), 10.0, 1.0e-10
    end

    test "zero mass gives zero force" do
      assert Motion.newtons_second_law(0, 9.81) == 0.0
    end
  end

  describe "momentum/2" do
    test "p = mv: 10 kg × 5 m/s = 50 kg·m/s" do
      assert_in_delta Motion.momentum(10, 5), 50.0, 1.0e-10
    end
  end

  describe "impulse/2" do
    test "J = Ft: 10 N × 5 s = 50 N·s" do
      assert_in_delta Motion.impulse(10, 5), 50.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Kinetic Energy
  # ---------------------------------------------------------------------------

  describe "kinetic_energy/2" do
    test "KE = ½mv²: 4 kg at 5 m/s = 50 J" do
      assert_in_delta Motion.kinetic_energy(4, 5), 50.0, 1.0e-10
    end

    test "zero velocity gives zero KE" do
      assert Motion.kinetic_energy(100, 0) == 0.0
    end
  end

  describe "rotational_kinetic_energy/2" do
    test "KE_rot = ½Iomega²: I=2, omega=3 → 9 J" do
      assert_in_delta Motion.rotational_kinetic_energy(2, 3), 9.0, 1.0e-10
    end
  end

  describe "total_kinetic_energy/4" do
    test "combines translational and rotational: 2 kg@3 m/s + I=4@5 rad/s = 9+50 = 59 J" do
      assert_in_delta Motion.total_kinetic_energy(2, 3, 4, 5), 59.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Centripetal Motion
  # ---------------------------------------------------------------------------

  describe "centripetal_force/3" do
    test "2 kg at 5 m/s radius 10 m → 5 N" do
      assert_in_delta Motion.centripetal_force(2, 5, 10), 5.0, 1.0e-10
    end
  end

  describe "centripetal_acceleration/2" do
    test "a_c = v²/r: 5 m/s, 10 m → 2.5 m/s²" do
      assert_in_delta Motion.centripetal_acceleration(5, 10), 2.5, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Moments of Inertia
  # ---------------------------------------------------------------------------

  describe "thin_disk_moment_of_inertia/2" do
    test "I = ½mR²: m=4, R=2 → 8 kg·m²" do
      assert_in_delta Motion.thin_disk_moment_of_inertia(4, 2), 8.0, 1.0e-10
    end
  end

  describe "thin_rod_center_moment_of_inertia/2" do
    test "I = mL²/12: m=6, L=2 → 2 kg·m²" do
      assert_in_delta Motion.thin_rod_center_moment_of_inertia(6, 2), 2.0, 1.0e-10
    end
  end

  describe "solid_sphere_moment_of_inertia/2" do
    test "I = 2mR²/5: m=5, R=1 → 2 kg·m²" do
      assert_in_delta Motion.solid_sphere_moment_of_inertia(5, 1), 2.0, 1.0e-10
    end
  end

  describe "parallel_axis/3" do
    test "I_cm + md²" do
      assert_in_delta Motion.parallel_axis(2.0, 3.0, 2.0), 14.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Projectile Motion
  # ---------------------------------------------------------------------------

  describe "projectile_range/3" do
    test "optimal launch angle 45° maximises range" do
      r45 = Motion.projectile_range(20, :math.pi() / 4)
      r30 = Motion.projectile_range(20, :math.pi() / 6)
      r60 = Motion.projectile_range(20, :math.pi() / 3)
      assert r45 >= r30
      assert r45 >= r60
    end

    test "symmetry: range at θ equals range at (90° - θ)" do
      r30 = Motion.projectile_range(20, :math.pi() / 6)
      r60 = Motion.projectile_range(20, :math.pi() / 3)
      assert_in_delta r30, r60, 1.0e-8
    end
  end

  describe "time_of_flight/3" do
    test "positive launch angle gives positive flight time" do
      assert Motion.time_of_flight(20, :math.pi() / 4) > 0
    end
  end

  describe "projectile_max_height/3" do
    test "90° launch gives maximum height" do
      h90 = Motion.projectile_max_height(20, :math.pi() / 2)
      h45 = Motion.projectile_max_height(20, :math.pi() / 4)
      assert h90 > h45
    end
  end

  # ---------------------------------------------------------------------------
  # SHM
  # ---------------------------------------------------------------------------

  describe "shm_displacement/4" do
    test "at t=0, φ=0: x = A" do
      assert_in_delta Motion.shm_displacement(0.1, 10, 0, 0), 0.1, 1.0e-10
    end

    test "displacement bounded by amplitude" do
      a = 0.05

      for t <- [0.0, 0.1, 0.5, 1.0, 2.0] do
        x = Motion.shm_displacement(a, 5, t)
        assert abs(x) <= a + 1.0e-10
      end
    end
  end

  describe "period_from_omega/1" do
    test "omega = 2π → T = 1 s" do
      assert_in_delta Motion.period_from_omega(2 * :math.pi()), 1.0, 1.0e-10
    end
  end

  describe "spring_angular_frequency/2" do
    test "k=10, m=2.5 → omega = 2 rad/s" do
      assert_in_delta Motion.spring_angular_frequency(10, 2.5), 2.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Analytical Mechanics
  # ---------------------------------------------------------------------------

  describe "lagrangian/2" do
    test "L = T - V: 50 - 20 = 30 J" do
      assert_in_delta Motion.lagrangian(50, 20), 30.0, 1.0e-10
    end
  end

  describe "hamiltons_equations/3" do
    test "SHM: {dp/dt, dq/dt} = {-omega²q, p}" do
      {dpdt, dqdt} = Motion.hamiltons_equations(2, 3, 1.5)
      assert_in_delta dpdt, -1.5 ** 2 * 3, 1.0e-10
      assert_in_delta dqdt, 2.0, 1.0e-10
    end
  end

  # ---------------------------------------------------------------------------
  # Kepler / Orbital
  # ---------------------------------------------------------------------------

  describe "keplers_third_law/2" do
    test "Earth's orbit (a = 1 AU, M = M☉) → T ≈ 1 year" do
      au = 1.496e11
      m_sun = 1.989e30
      t = Motion.keplers_third_law(au, m_sun)
      year = 3.15576e7
      assert_in_delta t / year, 1.0, 0.01
    end
  end

  describe "vis_viva/3" do
    test "circular orbit: v = √(GM/r) when a = r" do
      g = 6.67430e-11
      m = 5.972e24
      r = 6.771e6
      v_visviva = Motion.vis_viva(m, r, r)
      v_circular = Motion.orbital_speed(m, r)
      assert_in_delta v_visviva, v_circular, 1.0
    end
  end
end
