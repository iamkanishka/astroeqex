defmodule AstroEquations.Physics.MotionTest do
  use ExUnit.Case
  doctest Physics.Motion

  test "velocity calculation" do
    assert Physics.Motion.velocity(10, 2) == 5.0
    assert Physics.Motion.velocity(15, 3) == 5.0
  end

  test "acceleration calculation" do
    assert Physics.Motion.acceleration(20, 5) == 4.0
    assert Physics.Motion.acceleration(30, 6) == 5.0
  end

  test "Newton's Second Law" do
    assert Physics.Motion.newtons_second_law(5, 2) == 10.0
    assert Physics.Motion.newtons_second_law(3, 4) == 12.0
  end

  test "momentum calculation" do
    assert Physics.Motion.momentum(10, 5) == 50.0
    assert Physics.Motion.momentum(7, 3) == 21.0
  end

  test "impulse calculation" do
    assert Physics.Motion.impulse(10, 5) == 50.0
    assert Physics.Motion.impulse(8, 4) == 32.0
  end

  test "centripetal force calculation" do
    assert Physics.Motion.centripetal_force(2, 5, 10) == 5.0
    assert Physics.Motion.centripetal_force(3, 6, 9) == 12.0
  end

  test "kinetic energy calculation" do
    assert Physics.Motion.kinetic_energy(4, 5) == 50.0
    assert Physics.Motion.kinetic_energy(2, 10) == 100.0
  end
end
