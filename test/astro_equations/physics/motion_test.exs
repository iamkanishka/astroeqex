defmodule AstroEquations.Physics.MotionTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.Motion

  test "velocity calculation" do
    assert Motion.velocity(10, 2) == 5.0
    assert Motion.velocity(15, 3) == 5.0
  end

  test "acceleration calculation" do
    assert Motion.acceleration(20, 5) == 4.0
    assert Motion.acceleration(30, 6) == 5.0
  end

  test "Newton's Second Law" do
    assert Motion.newtons_second_law(5, 2) == 10.0
    assert Motion.newtons_second_law(3, 4) == 12.0
  end

  test "momentum calculation" do
    assert Motion.momentum(10, 5) == 50.0
    assert Motion.momentum(7, 3) == 21.0
  end

  test "impulse calculation" do
    assert Motion.impulse(10, 5) == 50.0
    assert Motion.impulse(8, 4) == 32.0
  end

  test "centripetal force calculation" do
    assert Motion.centripetal_force(2, 5, 10) == 5.0
    assert Motion.centripetal_force(3, 6, 9) == 12.0
  end

  test "kinetic energy calculation" do
    assert Motion.kinetic_energy(4, 5) == 50.0
    assert Motion.kinetic_energy(2, 10) == 100.0
  end

  test "angular velocity calculation" do
    assert Motion.angular_velocity(:math.pi(), 2) == :math.pi() / 2
  end

  test "angular acceleration calculation" do
    assert Motion.angular_acceleration(4, 2) == 2.0
  end

  test "point mass moment of inertia" do
    assert Motion.point_mass_moment_of_inertia(2, 3) == 18.0
  end

  test "multiple point masses moment of inertia" do
    assert Motion.multiple_point_masses_moment_of_inertia([1, 2], [3, 4]) == 35.0
  end

  test "thin disk moment of inertia" do
    assert Motion.thin_disk_moment_of_inertia(4, 2) == 8.0
  end

  test "thin loop moment of inertia" do
    assert Motion.thin_loop_moment_of_inertia(3, 2) == 12.0
  end

  test "thin rod center moment of inertia" do
    assert Motion.thin_rod_center_moment_of_inertia(6, 2) == 2.0
  end

  test "thin rod end moment of inertia" do
    assert Motion.thin_rod_end_moment_of_inertia(6, 2) == 8.0
  end

  test "Motional kinetic energy" do
    assert Motion.kinetic_energy(2, 3) == 9.0
  end

  test "total kinetic energy" do
    assert Motion.total_kinetic_energy(2, 3, 4, 5) == 59.0
  end

  test "angular momentum" do
    assert Motion.angular_momentum(3, 4) == 12.0
  end

  test "torque calculation" do
    assert Motion.torque(10, 2) == 20.0
    assert Motion.torque(10, 2, :math.pi() / 4) == 10 * 2 * :math.sin(:math.pi() / 4)
  end


    test "vertical velocity calculation" do
    assert_in_delta Physics.AdvancedMotion.vertical_velocity(10, -9.8, 5), 8.28, 0.01
  end

  test "horizontal displacement calculation" do
    assert Physics.AdvancedMotion.horizontal_displacement(15, 3) == 45.0
  end

  test "vertical displacement calculation" do
    assert Physics.AdvancedMotion.vertical_displacement(20, 2, -9.8) == 20.4
  end

  test "lagrangian calculation" do
    assert Physics.AdvancedMotion.lagrangian(50, 20) == 30.0
  end

  test "generalized momentum approximation" do
    # Testing with a simple quadratic Lagrangian L = 0.5 * q̇²
    # ∂L/∂q̇ = q̇, so at q̇=2 should be 2
    result = Physics.AdvancedMotion.generalized_momentum(fn x -> 0.5 * x ** 2 end, 2)
    assert_in_delta result, 2.0, 0.001
  end

  test "hamiltonian calculation" do
    assert Physics.AdvancedMotion.hamiltonian([2, 3], [1.5, 2], 10) == 2.0
  end

  test "hamilton's equations for harmonic oscillator" do
    {dpdt, dqdt} = Physics.AdvancedMotion.hamiltons_equations(2, 3, 1.5)
    assert_in_delta dpdt, -6.75, 0.001
    assert dqdt == 2
  end

end
