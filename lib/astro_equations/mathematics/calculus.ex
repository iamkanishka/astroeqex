defmodule AstroEquations.Mathematics.Calculus do
  @moduledoc """
  Numerical calculus utilities commonly used in astronomical computations.

  Provides:
  - Numerical differentiation (first and second derivatives, forward/central/backward)
  - Numerical integration (rectangle, trapezoid, Simpson's rule, Gaussian quadrature)
  - Root finding (bisection, Newton-Raphson, secant method)
  - ODE integration (Euler, RK4)
  - Gradient and Laplacian (finite differences)
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "A real-valued function of one real variable."
  @type unary_fn :: (number() -> number())

  @typedoc "A real-valued function of two real variables (e.g. ODE right-hand side)."
  @type binary_fn :: (number(), number() -> number())

  @typedoc "Step size for finite-difference approximations. Must be positive."
  @type step :: float()

  # ---------------------------------------------------------------------------
  # Numerical Differentiation
  # ---------------------------------------------------------------------------

  @doc """
  Central difference first derivative: f'(x) ≈ (f(x+h) - f(x-h)) / (2h).

  Most accurate single-point finite difference for smooth functions.

  ## Parameters
    - f: Function f(x)
    - x: Point at which to evaluate
    - h: Step size (default: 1.0e-7)

  ## Returns
    Approximate derivative f'(x)

  ## Examples
      iex> Calculus.derivative(&:math.sin/1, 0.0) |> Float.round(6)
      1.0
  """

  @spec derivative((number -> number), number, number) :: float
  def derivative(f, x, h \\ 1.0e-7) do
    (f.(x + h) - f.(x - h)) / (2 * h)
  end

  @doc """
  Forward difference first derivative: f'(x) ≈ (f(x+h) - f(x)) / h.

  Lower accuracy but useful at boundary points.

  ## Parameters
    - f: Function f(x)
    - x: Evaluation point
    - h: Step size (default: 1.0e-7)

  ## Examples
      iex> Calculus.forward_difference(&:math.sin/1, 0.0) |> Float.round(4)
      1.0
  """

  @spec forward_difference((number -> number), number, number) :: float
  def forward_difference(f, x, h \\ 1.0e-7), do: (f.(x + h) - f.(x)) / h

  @doc """
  Backward difference first derivative: f'(x) ≈ (f(x) - f(x-h)) / h.

  ## Parameters
    - f: Function f(x)
    - x: Evaluation point
    - h: Step size (default: 1.0e-7)

  ## Examples
      iex> Calculus.backward_difference(&:math.sin/1, :math.pi()/2) |> Float.round(4)
      0.0
  """

  @spec backward_difference((number -> number), number, number) :: float
  def backward_difference(f, x, h \\ 1.0e-7), do: (f.(x) - f.(x - h)) / h

  @doc """
  Central difference second derivative: f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h².

  ## Parameters
    - f: Function f(x)
    - x: Evaluation point
    - h: Step size (default: 1.0e-5)

  ## Examples
      iex> Calculus.second_derivative(&:math.sin/1, 0.0) |> Float.round(6)
      0.0
  """

  @spec second_derivative((number -> number), number, number) :: float
  def second_derivative(f, x, h \\ 1.0e-5) do
    (f.(x + h) - 2 * f.(x) + f.(x - h)) / :math.pow(h, 2)
  end

  # ---------------------------------------------------------------------------
  # Numerical Integration
  # ---------------------------------------------------------------------------

  @doc """
  Trapezoidal rule integration over [a, b] with n equal subintervals.

  ∫_a^b f(x) dx ≈ h/2 * [f(a) + 2 f(x_1) + ... + 2 f(x_{n-1}) + f(b)]

  ## Parameters
    - f: Integrand f(x)
    - a: Lower bound
    - b: Upper bound
    - n: Number of subintervals (default: 1000)

  ## Returns
    Approximate integral

  ## Examples
      iex> Calculus.trapezoid(&:math.sin/1, 0, :math.pi) |> Float.round(6)
      2.0
  """

  @spec trapezoid((number -> number), number, number, pos_integer) :: float
  def trapezoid(f, a, b, n \\ 1_000) do
    h = (b - a) / n

    interior_sum =
      Enum.reduce(1..(n - 1), 0.0, fn i, acc ->
        acc + f.(a + i * h)
      end)

    h / 2 * (f.(a) + 2 * interior_sum + f.(b))
  end

  @doc """
  Simpson's 1/3 rule integration over [a, b] (n must be even).

  Higher accuracy than the trapezoidal rule for smooth functions.

  ## Parameters
    - f: Integrand f(x)
    - a: Lower bound
    - b: Upper bound
    - n: Number of subintervals — must be even (default: 1000)

  ## Returns
    Approximate integral

  ## Examples
      iex> Calculus.simpsons(&:math.sin/1, 0, :math.pi) |> Float.round(8)
      2.0
  """

  @spec simpsons((number -> number), number, number, pos_integer) :: float
  def simpsons(f, a, b, n \\ 1_000) do
    n_even = if rem(n, 2) == 0, do: n, else: n + 1
    h = (b - a) / n_even
    {odd_sum, even_sum} = simpsons_sums(f, a, h, n_even)
    h / 3 * (f.(a) + 4 * odd_sum + 2 * even_sum + f.(b))
  end

  defp simpsons_sums(f, a, h, n) do
    Enum.reduce(1..(n - 1), {0.0, 0.0}, fn i, {odd, even} ->
      x = a + i * h
      if rem(i, 2) == 1, do: {odd + f.(x), even}, else: {odd, even + f.(x)}
    end)
  end

  @doc """
  Rectangle (midpoint) rule integration over [a, b].

  ## Parameters
    - f: Integrand f(x)
    - a: Lower bound
    - b: Upper bound
    - n: Number of subintervals (default: 1000)

  ## Examples
      iex> Calculus.rectangle(&:math.sin/1, 0, :math.pi(), 1000) |> Float.round(4)
      2.0
  """

  @spec rectangle((number -> number), number, number, pos_integer) :: float
  def rectangle(f, a, b, n \\ 1_000) do
    h = (b - a) / n

    Enum.reduce(0..(n - 1), 0.0, fn i, acc ->
      x_mid = a + (i + 0.5) * h
      acc + f.(x_mid) * h
    end)
  end

  # ---------------------------------------------------------------------------
  # Root Finding
  # ---------------------------------------------------------------------------

  @doc """
  Bisection method for root finding on [a, b].

  Requires f(a) and f(b) to have opposite signs.

  ## Parameters
    - f:        Function f(x)
    - a:        Left bracket
    - b:        Right bracket
    - tol:      Tolerance (default: 1.0e-10)
    - max_iter: Maximum iterations (default: 1000)

  ## Returns
    Approximate root, or `{:error, reason}`

  ## Examples
      iex> Calculus.bisection(&:math.sin/1, 3.0, 4.0) |> Float.round(10)
      3.1415926536
  """

  @spec bisection((number -> number), number, number, number, pos_integer) ::
          float | {:error, String.t()}
  def bisection(f, a, b, tol \\ 1.0e-10, max_iter \\ 1_000) do
    if f.(a) * f.(b) > 0 do
      {:error, "f(a) and f(b) must have opposite signs"}
    else
      do_bisection(f, a, b, tol, max_iter, 0)
    end
  end

  defp do_bisection(_f, a, b, tol, _max_iter, _iter) when abs(b - a) < tol, do: (a + b) / 2

  defp do_bisection(_f, _a, _b, _tol, max_iter, iter) when iter >= max_iter,
    do: {:error, "Maximum iterations reached"}

  defp do_bisection(f, a, b, tol, max_iter, iter) do
    mid = (a + b) / 2

    if f.(a) * f.(mid) < 0,
      do: do_bisection(f, a, mid, tol, max_iter, iter + 1),
      else: do_bisection(f, mid, b, tol, max_iter, iter + 1)
  end

  @doc """
  Newton-Raphson root finding from an initial guess x₀.

  Iterates x_{n+1} = x_n - f(x_n) / f'(x_n) until convergence.

  ## Parameters
    - f:        Function f(x)
    - x0:       Initial guess
    - tol:      Tolerance (default: 1.0e-10)
    - max_iter: Maximum iterations (default: 100)

  ## Returns
    Approximate root, or `{:error, reason}`

  ## Examples
      iex> Calculus.newton_raphson(fn x -> x*x - 2 end, 1.5) |> Float.round(6)
      1.414214
  """

  @spec newton_raphson((number -> number), number, number, pos_integer) ::
          number() | {:error, String.t()}
  def newton_raphson(f, x0, tol \\ 1.0e-10, max_iter \\ 100) do
    do_newton(f, x0, tol, max_iter, 0)
  end

  defp do_newton(_f, x, tol, _max_iter, _iter) when abs(x) < tol, do: x

  defp do_newton(_f, _x, _tol, max_iter, iter) when iter >= max_iter,
    do: {:error, "Maximum iterations reached"}

  defp do_newton(f, x, tol, max_iter, iter) do
    fx = f.(x)
    fpx = derivative(f, x)

    if abs(fpx) < 1.0e-15,
      do: {:error, "Derivative near zero"},
      else: do_newton(f, x - fx / fpx, tol, max_iter, iter + 1)
  end

  # ---------------------------------------------------------------------------
  # ODE Integration
  # ---------------------------------------------------------------------------

  @doc """
  Euler method for first-order ODE: dy/dx = f(x, y).

  ## Parameters
    - f:     RHS function f(x, y)
    - x0:    Initial x value
    - y0:    Initial y value
    - x_end: End x value
    - n:     Number of steps (default: 1000)

  ## Returns
    List of `{x, y}` pairs at each step

  ## Examples
      iex> Calculus.euler_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 1000) |> Enum.count() > 0
      true
  """

  @spec euler_ode((number, number -> number), number, number, number, pos_integer) :: [
          {float, float}
        ]
  def euler_ode(f, x0, y0, x_end, n \\ 1_000) do
    h = (x_end - x0) / n

    Enum.scan(0..(n - 1), {x0, y0}, fn _, {x, y} ->
      {x + h, y + h * f.(x, y)}
    end)
  end

  @doc """
  Fourth-order Runge-Kutta ODE solver for dy/dx = f(x, y).

  ## Parameters
    - f:     RHS function f(x, y)
    - x0:    Initial x value
    - y0:    Initial y value
    - x_end: End x value
    - n:     Number of steps (default: 1000)

  ## Returns
    List of `{x, y}` pairs at each step

  ## Examples
      iex> Calculus.rk4_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 10) |> length()
      11
  """

  @spec rk4_ode((number, number -> number), number, number, number, pos_integer) :: [
          {float, float}
        ]
  def rk4_ode(f, x0, y0, x_end, n \\ 1_000) do
    h = (x_end - x0) / n
    Enum.scan(0..(n - 1), {x0, y0}, fn _, state -> rk4_step(f, state, h) end)
  end

  defp rk4_step(f, {x, y}, h) do
    k1 = f.(x, y)
    k2 = f.(x + h / 2, y + h * k1 / 2)
    k3 = f.(x + h / 2, y + h * k2 / 2)
    k4 = f.(x + h, y + h * k3)
    {x + h, y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6}
  end

  # ---------------------------------------------------------------------------
  # Multi-Dimensional Calculus
  # ---------------------------------------------------------------------------

  @doc """
  Numerical gradient vector at a point using central differences.

  Returns [∂f/∂x₁, ∂f/∂x₂, …] for an n-dimensional function.

  ## Parameters
    - f:  Function accepting a list of coordinates and returning a number
    - xs: Coordinate list [x₁, x₂, …]
    - h:  Step size (default: 1.0e-5)

  ## Examples
      iex> Calculus.gradient(fn [x, y] -> x*x + y*y end, [1.0, 2.0]) |> Enum.map(&Float.round(&1, 4))
      [2.0, 4.0]
  """

  @spec gradient((list(number) -> number), [number], number) :: [float]
  def gradient(f, xs, h \\ 1.0e-5) do
    Enum.with_index(xs)
    |> Enum.map(fn {_, i} ->
      xs_p = List.replace_at(xs, i, Enum.at(xs, i) + h)
      xs_m = List.replace_at(xs, i, Enum.at(xs, i) - h)
      (f.(xs_p) - f.(xs_m)) / (2 * h)
    end)
  end

  @doc """
  Numerical Laplacian (sum of second partial derivatives): Δf = Σ ∂²f/∂xᵢ²

  ## Parameters
    - f:  Function accepting a list of coordinates
    - xs: Coordinate list
    - h:  Step size (default: 1.0e-5)

  ## Examples
      iex> Calculus.laplacian(fn [x, y] -> x*x + y*y end, [1.0, 1.0]) |> Float.round(4)
      4.0
  """

  @spec laplacian((list(number) -> number), [number], number) :: float
  def laplacian(f, xs, h \\ 1.0e-5) do
    f0 = f.(xs)

    Enum.with_index(xs)
    |> Enum.reduce(0.0, fn {_, i}, acc ->
      xi = Enum.at(xs, i)
      fp = f.(List.replace_at(xs, i, xi + h))
      fm = f.(List.replace_at(xs, i, xi - h))
      acc + (fp - 2 * f0 + fm) / (h * h)
    end)
  end

  @doc """
  5-point Gaussian quadrature over [a, b] (exact for polynomials up to degree 9).

  ## Parameters
    - f: Integrand f(x)
    - a: Lower bound
    - b: Upper bound

  ## Examples
      iex> Calculus.gaussian_quadrature(&:math.sin/1, 0, :math.pi()) |> Float.round(6)
      2.0
  """

  @spec gaussian_quadrature((number -> number), number, number) :: float
  def gaussian_quadrature(f, a, b) do
    # 5-point Gauss-Legendre nodes and weights on [-1, 1]
    nodes = [0.0, 0.5384693101, -0.5384693101, 0.9061798459, -0.9061798459]
    weights = [0.5688888889, 0.4786286705, 0.4786286705, 0.2369268851, 0.2369268851]

    mid = (a + b) / 2
    half = (b - a) / 2

    Enum.zip(nodes, weights)
    |> Enum.reduce(0.0, fn {t, w}, acc ->
      acc + w * f.(mid + half * t)
    end)
    |> Kernel.*(half)
  end

  @doc """
  Secant method root finding — derivative-free alternative to Newton-Raphson.

  ## Parameters
    - f:        Function f(x)
    - x0, x1:   Two initial guesses
    - tol:      Tolerance (default: 1.0e-10)
    - max_iter: Maximum iterations (default: 100)

  ## Returns
    Approximate root, or `{:error, reason}`

  ## Examples
      iex> Calculus.secant(fn x -> x*x - 2 end, 1.0, 1.5) |> Float.round(6)
      1.414214
  """

  @spec secant((number -> number), number, number, number, pos_integer) ::
          number() | {:error, String.t()}
  def secant(f, x0, x1, tol \\ 1.0e-10, max_iter \\ 100) do
    do_secant(f, x0, x1, tol, max_iter, 0)
  end

  defp do_secant(_f, _x0, x1, tol, _max_iter, _iter) when abs(x1) < tol, do: x1

  defp do_secant(_f, _x0, _x1, _tol, max_iter, iter) when iter >= max_iter,
    do: {:error, "Maximum iterations reached"}

  defp do_secant(f, x0, x1, tol, max_iter, iter) do
    fx0 = f.(x0)
    fx1 = f.(x1)
    denom = fx1 - fx0

    if abs(denom) < 1.0e-15 do
      {:error, "Denominator near zero"}
    else
      x2 = x1 - fx1 * (x1 - x0) / denom
      do_secant(f, x1, x2, tol, max_iter, iter + 1)
    end
  end

  @doc """
  RK4 solver for a system of first-order ODEs: dy/dt = f(t, y) where y is a vector.

  ## Parameters
    - f:     RHS function f(t, [y₁, y₂, …]) → [dy₁/dt, dy₂/dt, …]
    - t0:    Initial time
    - y0:    Initial state vector
    - t_end: End time
    - n:     Number of steps (default: 1000)

  ## Returns
    List of {t, y_vector} pairs

  ## Examples
      iex> Calculus.rk4_system(fn _t, [y] -> [-y] end, 0, [1.0], 1.0, 10) |> length()
      11
  """

  @spec rk4_system(
          (number, [number] -> [number]),
          number,
          [number],
          number,
          pos_integer
        ) :: [{float, [float]}]
  def rk4_system(f, t0, y0, t_end, n \\ 1_000) do
    h = (t_end - t0) / n
    vec_add = fn a, b -> Enum.zip_with(a, b, fn x, y -> x + y end) end
    scale = fn s, v -> Enum.map(v, fn x -> s * x end) end

    Enum.scan(0..(n - 1), {t0, y0}, fn _, {t, y} ->
      k1 = f.(t, y)
      k2 = f.(t + h / 2, vec_add.(y, scale.(h / 2, k1)))
      k3 = f.(t + h / 2, vec_add.(y, scale.(h / 2, k2)))
      k4 = f.(t + h, vec_add.(y, scale.(h, k3)))

      step =
        Enum.zip_with([k1, k2, k3, k4], fn [a, b, c, d] -> h * (a + 2 * b + 2 * c + d) / 6 end)

      {t + h, vec_add.(y, step)}
    end)
  end
end
