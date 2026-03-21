defmodule AstroEquations.Mathematics.CalculusTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Mathematics.Calculus

  @tol 1.0e-6

  # ---------------------------------------------------------------------------
  # Numerical Differentiation
  # ---------------------------------------------------------------------------

  describe "derivative/3" do
    test "d/dx sin(x) at x=0 ≈ cos(0) = 1" do
      assert_in_delta Calculus.derivative(&:math.sin/1, 0.0), 1.0, @tol
    end

    test "d/dx sin(x) at x=π/2 ≈ 0" do
      assert_in_delta Calculus.derivative(&:math.sin/1, :math.pi() / 2), 0.0, @tol
    end

    test "d/dx exp(x) at x=0 ≈ 1" do
      assert_in_delta Calculus.derivative(&:math.exp/1, 0.0), 1.0, @tol
    end

    test "d/dx exp(x) at x=1 ≈ e" do
      assert_in_delta Calculus.derivative(&:math.exp/1, 1.0), :math.exp(1), @tol
    end

    test "d/dx x² = 2x at x=3 → 6" do
      assert_in_delta Calculus.derivative(fn x -> x * x end, 3.0), 6.0, @tol
    end
  end

  describe "forward_difference/3" do
    test "sin'(0) ≈ 1" do
      assert_in_delta Calculus.forward_difference(&:math.sin/1, 0.0), 1.0, 1.0e-5
    end

    test "less accurate than central difference" do
      err_c = abs(Calculus.derivative(&:math.sin/1, 1.0) - :math.cos(1.0))
      err_f = abs(Calculus.forward_difference(&:math.sin/1, 1.0) - :math.cos(1.0))
      assert err_f >= err_c
    end
  end

  describe "backward_difference/3" do
    test "sin'(π/2) ≈ 0" do
      assert_in_delta Calculus.backward_difference(&:math.sin/1, :math.pi() / 2), 0.0, 1.0e-5
    end
  end

  describe "second_derivative/3" do
    test "d²/dx² sin(x) at x=0 ≈ 0" do
      assert_in_delta Calculus.second_derivative(&:math.sin/1, 0.0), 0.0, 1.0e-4
    end

    test "d²/dx² sin(x) at x=π/2 ≈ -1" do
      assert_in_delta Calculus.second_derivative(&:math.sin/1, :math.pi() / 2), -1.0, 1.0e-4
    end

    test "d²/dx² x³ at x=2 = 6x = 12" do
      f = fn x -> x * x * x end
      assert_in_delta Calculus.second_derivative(f, 2.0), 12.0, 1.0e-4
    end
  end

  # ---------------------------------------------------------------------------
  # Numerical Integration
  # ---------------------------------------------------------------------------

  describe "trapezoid/4" do
    test "∫₀^π sin(x) dx = 2" do
      result = Calculus.trapezoid(&:math.sin/1, 0, :math.pi(), 10_000)
      assert_in_delta result, 2.0, 1.0e-6
    end

    test "∫₀¹ x dx = 0.5" do
      result = Calculus.trapezoid(fn x -> x end, 0, 1, 10_000)
      assert_in_delta result, 0.5, 1.0e-6
    end

    test "∫₀¹ x² dx = 1/3" do
      result = Calculus.trapezoid(fn x -> x * x end, 0, 1, 10_000)
      assert_in_delta result, 1.0 / 3.0, 1.0e-6
    end
  end

  describe "simpsons/4" do
    test "∫₀^π sin(x) dx = 2 (very accurate)" do
      result = Calculus.simpsons(&:math.sin/1, 0, :math.pi(), 1000)
      assert_in_delta result, 2.0, 1.0e-10
    end

    test "∫₀¹ x³ dx = 0.25" do
      result = Calculus.simpsons(fn x -> x * x * x end, 0, 1, 1000)
      assert_in_delta result, 0.25, 1.0e-8
    end

    test "more accurate than trapezoid for same n" do
      f = fn x -> :math.sin(x) end
      trap = Calculus.trapezoid(f, 0, :math.pi(), 100)
      simp = Calculus.simpsons(f, 0, :math.pi(), 100)
      err_trap = abs(trap - 2.0)
      err_simp = abs(simp - 2.0)
      assert err_simp < err_trap
    end
  end

  describe "rectangle/4" do
    test "∫₀¹ 1 dx = 1" do
      result = Calculus.rectangle(fn _ -> 1.0 end, 0, 1, 1000)
      assert_in_delta result, 1.0, 1.0e-6
    end

    test "∫₀¹ x dx ≈ 0.5" do
      result = Calculus.rectangle(fn x -> x end, 0, 1, 10_000)
      assert_in_delta result, 0.5, 1.0e-4
    end
  end

  # ---------------------------------------------------------------------------
  # Root Finding
  # ---------------------------------------------------------------------------

  describe "bisection/5" do
    test "finds root of sin(x) in [3, 4] → π" do
      result = Calculus.bisection(&:math.sin/1, 3.0, 4.0)
      assert_in_delta result, :math.pi(), 1.0e-8
    end

    test "finds root of x² - 2 → √2" do
      result = Calculus.bisection(fn x -> x * x - 2 end, 1.0, 2.0)
      assert_in_delta result, :math.sqrt(2), 1.0e-8
    end

    test "returns error when f(a) and f(b) have same sign" do
      assert {:error, _} = Calculus.bisection(&:math.sin/1, 0.1, 0.5)
    end
  end

  describe "newton_raphson/4" do
    test "finds root of x² - 2 → √2" do
      result = Calculus.newton_raphson(fn x -> x * x - 2 end, 1.5)
      assert_in_delta result, :math.sqrt(2), 1.0e-8
    end

    test "finds root of sin(x) near π" do
      result = Calculus.newton_raphson(&:math.sin/1, 3.0)
      assert_in_delta result, :math.pi(), 1.0e-8
    end

    test "finds root of e^x - 2 → ln(2)" do
      result = Calculus.newton_raphson(fn x -> :math.exp(x) - 2 end, 1.0)
      assert_in_delta result, :math.log(2), 1.0e-8
    end
  end

  # ---------------------------------------------------------------------------
  # ODE Solvers
  # ---------------------------------------------------------------------------

  describe "euler_ode/5" do
    test "dy/dx = -y, y(0) = 1 → e^(-x) at x=1 ≈ 0.368" do
      results = Calculus.euler_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 10_000)
      {_, y_final} = List.last(results)
      assert_in_delta y_final, :math.exp(-1), 0.001
    end

    test "dy/dx = 0, y(0) = 5 → y stays 5" do
      results = Calculus.euler_ode(fn _x, _y -> 0 end, 0, 5.0, 1.0, 100)
      {_, y_final} = List.last(results)
      assert_in_delta y_final, 5.0, 1.0e-8
    end
  end

  describe "rk4_ode/5" do
    test "dy/dx = -y → e^(-x), much more accurate than Euler" do
      results = Calculus.rk4_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 100)
      {_, y_rk4} = List.last(results)
      assert_in_delta y_rk4, :math.exp(-1), 1.0e-8
    end

    test "dy/dx = x → x²/2 + C at x=2 starting from (0, 0)" do
      results = Calculus.rk4_ode(fn x, _y -> x end, 0, 0.0, 2.0, 1000)
      {_, y_final} = List.last(results)
      assert_in_delta y_final, 2.0, 0.001
    end

    test "returns list of {x, y} pairs" do
      results = Calculus.rk4_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 10)
      assert is_list(results)
      assert Enum.all?(results, fn {x, y} -> is_number(x) and is_number(y) end)
    end
  end
end
