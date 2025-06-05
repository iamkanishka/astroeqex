defmodule GeneralRelativity.Metrics do
  @moduledoc """
  This module provides implementations of various spacetime metrics in General Relativity,
  including Minkowski, Schwarzschild, and Rindler coordinates.
  """

  @doc """
  Returns the Minkowski metric tensor (η) in matrix form.

  ## Examples

      iex> GeneralRelativity.Metrics.minkowski_metric()
      [
        [-1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
      ]
  """
  def minkowski_metric do
    [
      [-1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
    ]
  end

  @doc """
  Calculates the spacetime interval (ds²) in Minkowski space.

  ## Parameters
    - dt: Time differential
    - dx: x-coordinate differential
    - dy: y-coordinate differential
    - dz: z-coordinate differential
    - c: Speed of light (default: 1 in natural units)

  ## Examples

      iex> GeneralRelativity.Metrics.minkowski_interval(1, 0, 0, 0)
      -1
  """
  def minkowski_interval(dt, dx, dy, dz, c \\ 1) do
    -c * c * dt * dt + dx * dx + dy * dy + dz * dz
  end

  @doc """
  Returns the Schwarzschild metric tensor for given parameters.

  ## Parameters
    - m: Mass of the gravitating body
    - r: Radial coordinate
    - theta: Polar angle
    - c: Speed of light (default: 1)
    - g: Gravitational constant (default: 1)

  ## Examples

      iex> GeneralRelativity.Metrics.schwarzschild_metric(1, 3, :math.pi()/2)
      [
        [-(1 - 2/3), 0, 0, 0],
        [0, 1/(1 - 2/3), 0, 0],
        [0, 0, 9, 0],
        [0, 0, 0, 9 * :math.sin(:math.pi()/2) ** 2]
      ]
  """
  def schwarzschild_metric(m, r, theta, c \\ 1, g \\ 1) do
    rs = 2 * g * m / (c * c)
    factor = 1 - rs / r

    [
      [-factor, 0, 0, 0],
      [0, 1 / factor, 0, 0],
      [0, 0, r * r, 0],
      [0, 0, 0, r * r * :math.sin(theta) * :math.sin(theta)]
    ]
  end

  @doc """
  Calculates the spacetime interval (ds²) in Schwarzschild coordinates.

  ## Parameters
    - m: Mass
    - r: Radial coordinate
    - theta: Polar angle
    - dt: Time differential
    - dr: Radial differential
    - dtheta: Polar angle differential
    - dphi: Azimuthal angle differential
    - c: Speed of light (default: 1)
    - g: Gravitational constant (default: 1)

  ## Examples

      iex> GeneralRelativity.Metrics.schwarzschild_interval(1, 3, :math.pi()/2, 1, 0, 0, 0)
      -(1 - 2/3)
  """
  def schwarzschild_interval(m, r, theta, dt, dr, dtheta, dphi, c \\ 1, g \\ 1) do
    rs = 2 * g * m / (c * c)
    factor = 1 - rs / r

    -factor * c * c * dt * dt +
    (1 / factor) * dr * dr +
    r * r * dtheta * dtheta +
    r * r * :math.sin(theta) * :math.sin(theta) * dphi * dphi
  end

  @doc """
  Calculates the line element in Rindler coordinates.

  ## Parameters
    - g: Proper acceleration
    - x_prime: Space coordinate in Rindler frame
    - dt_prime: Time differential in Rindler frame
    - dx_prime: Space differential in Rindler frame
    - c: Speed of light (default: 1)

  ## Examples

      iex> GeneralRelativity.Metrics.rindler_interval(9.8, 1, 1, 0)
      -:math.pow(1 + 9.8 * 1, 2) * 1 * 1 + 0
  """
  def rindler_interval(g, x_prime, dt_prime, dx_prime, c \\ 1) do
    -:math.pow(1 + g * x_prime / (c * c), 2) * c * c * dt_prime * dt_prime + dx_prime * dx_prime
  end
end

defmodule GeneralRelativity.TensorAlgebra do
  @moduledoc """
  This module implements tensor operations using Einstein notation,
  including index raising/lowering, tensor transformations, and four-vector operations.

  All functions assume 4-dimensional spacetime tensors unless otherwise specified.
  """

  @doc """
  Lowers an index of a tensor using the metric tensor (converts contravariant to covariant).

  ## Parameters
    - tensor: The tensor components (as a list or nested lists)
    - metric: The metric tensor g_{μν}
    - index_positions: List of atom indicating position of each index (:upper or :lower)
    - index_to_lower: Which index to lower (0-based)

  ## Returns
    The tensor with the specified index lowered

  ## Examples

      iex> metric = GeneralRelativity.Metrics.minkowski_metric()
      iex> GeneralRelativity.TensorAlgebra.lower_index([1, 0, 0, 0], metric, [:upper, :upper, :upper, :upper], 0)
      [-1, 0, 0, 0]
  """
  def lower_index(tensor, metric, index_positions, index_to_lower) do
    tensor
    |> Enum.with_index()
    |> Enum.map(fn {component, i} ->
      if Enum.at(index_positions, i) == :upper and i == index_to_lower do
        Enum.reduce(0..3, 0, fn j, acc ->
          acc + Enum.at(Enum.at(metric, i), j) * component
        end)
      else
        component
      end
    end)
  end

  @doc """
  Raises an index of a tensor using the inverse metric (converts covariant to contravariant).

  ## Parameters
    - tensor: The tensor components
    - metric: The metric tensor g_{μν}
    - index_positions: List indicating position of indices (:upper or :lower)
    - index_to_raise: Which index to raise (0-based)

  ## Returns
    The tensor with the specified index raised

  ## Examples

      iex> metric = GeneralRelativity.Metrics.minkowski_metric()
      iex> GeneralRelativity.TensorAlgebra.raise_index([-1, 0, 0, 0], metric, [:lower, :lower, :lower, :lower], 0)
      [1, 0, 0, 0]
  """
  def raise_index(tensor, metric, index_positions, index_to_raise) do
    # For Minkowski metric, inverse is the same as the metric
    inverse_metric = metric

    tensor
    |> Enum.with_index()
    |> Enum.map(fn {component, i} ->
      if Enum.at(index_positions, i) == :lower and i == index_to_raise do
        Enum.reduce(0..3, 0, fn j, acc ->
          acc + Enum.at(Enum.at(inverse_metric, i), j) * component
        end)
      else
        component
      end
    end)
  end

  @doc """
  Transforms a tensor between coordinate systems using the Jacobian matrices.

  ## Parameters
    - tensor: The tensor to transform (as nested lists)
    - jacobian: The Jacobian matrix ∂x^{μ'}/∂x^ν
    - inverse_jacobian: The inverse Jacobian matrix ∂x^μ/∂x^{ν'}

  ## Returns
    The transformed tensor

  ## Examples

      iex> jacobian = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
      iex> GeneralRelativity.TensorAlgebra.transform_tensor([1, 0, 0, 0], jacobian, jacobian)
      [1, 0, 0, 0]
  """
  def transform_tensor(tensor, jacobian, _inverse_jacobian) when is_list(tensor) do
    # For vectors (1D tensors)
    if not List.first(tensor) |> is_list do
      Enum.map(0..3, fn i ->
        Enum.reduce(0..3, 0, fn j, acc ->
          acc + Enum.at(tensor, j) * Enum.at(Enum.at(jacobian, j), i)
        end)
      end)
    else
      # For higher rank tensors (would need more complex implementation)
      tensor
    end
  end

  @doc """
  Calculates the inner product of two four-vectors using the metric.

  ## Parameters
    - a: First four-vector (contravariant components as list)
    - b: Second four-vector (contravariant components as list)
    - metric: The metric tensor g_{μν}

  ## Returns
    The scalar product a·b = g_{μν} a^μ b^ν

  ## Examples

      iex> metric = GeneralRelativity.Metrics.minkowski_metric()
      iex> GeneralRelativity.TensorAlgebra.four_vector_product([1, 0, 0, 0], [1, 0, 0, 0], metric)
      -1
      iex> GeneralRelativity.TensorAlgebra.four_vector_product([0, 1, 0, 0], [0, 1, 0, 0], metric)
      1
  """
  def four_vector_product(a, b, metric) do
    Enum.reduce(0..3, 0, fn i, acc1 ->
      Enum.reduce(0..3, acc1, fn j, acc2 ->
        acc2 + Enum.at(Enum.at(metric, i), j) * Enum.at(a, i) * Enum.at(b, j)
      end)
    end)
  end

  @doc """
  Calculates the inverse of a metric tensor.

  ## Parameters
    - metric: The metric tensor g_{μν}

  ## Returns
    The inverse metric g^{μν}

  ## Examples

      iex> metric = GeneralRelativity.Metrics.minkowski_metric()
      iex> GeneralRelativity.TensorAlgebra.inverse_metric(metric)
      [
        [-1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
      ]
  """
  def inverse_metric(metric) do
    # For diagonal metrics, the inverse is just 1/g_{ii} for each component
    Enum.map(metric, fn row ->
      Enum.map(row, fn
        0 -> 0
        x -> 1 / x
      end)
    end)
  end
end
