defmodule AstroEquations.Physics.GeneralRelativityTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.GeneralRelativity

  setup do
    minkowski_metric = [
      [-1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
    ]

    schwarzschild_metric = [
      [-0.5, 0, 0, 0],
      [0, 2.0, 0, 0],
      [0, 0, 9, 0],
      [0, 0, 0, 9]
    ]

    {:ok, mink: minkowski_metric, schw: schwarzschild_metric}
  end

  test "four_vector_product with Minkowski metric", %{mink: metric} do
    # Time-like vector with itself
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([1, 0, 0, 0], [1, 0, 0, 0], metric) == -1

    # Space-like vector with itself
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([0, 1, 0, 0], [0, 1, 0, 0], metric) == 1

    # Light-like vector
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([1, 1, 0, 0], [1, 1, 0, 0], metric) == 0

    # Different vectors
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([1, 1, 0, 0], [1, -1, 0, 0], metric) == -2
  end

  test "four_vector_product with Schwarzschild metric", %{schw: metric} do
    # Time component only
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([1, 0, 0, 0], [1, 0, 0, 0], metric) == -0.5

    # Radial component only
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([0, 1, 0, 0], [0, 1, 0, 0], metric) == 2.0

    # Angular component
    assert AstroEquations.Physics.GeneralRelativity.four_vector_product([0, 0, 0, 1], [0, 0, 0, 1], metric) == 9
  end

  test "lower_index with Minkowski metric", %{mink: metric} do
    # Lower time component of contravariant vector
    assert AstroEquations.Physics.GeneralRelativity.lower_index([1, 0, 0, 0], metric, [:upper, :upper, :upper, :upper], 0) ==
           [-1, 0, 0, 0]

    # Lower spatial component
    assert AstroEquations.Physics.GeneralRelativity.lower_index([0, 1, 0, 0], metric, [:upper, :upper, :upper, :upper], 1) ==
           [0, 1, 0, 0]  # No change for spatial in Minkowski

    # Try to lower already lowered index
    assert AstroEquations.Physics.GeneralRelativity.lower_index([-1, 0, 0, 0], metric, [:lower, :upper, :upper, :upper], 0) ==
           [-1, 0, 0, 0]  # Shouldn't change
  end

  test "raise_index with Minkowski metric", %{mink: metric} do
    # Raise time component of covariant vector
    assert AstroEquations.Physics.GeneralRelativity.raise_index([-1, 0, 0, 0], metric, [:lower, :lower, :lower, :lower], 0) ==
           [1, 0, 0, 0]

    # Raise spatial component
    assert AstroEquations.Physics.GeneralRelativity.raise_index([0, 1, 0, 0], metric, [:lower, :lower, :lower, :lower], 1) ==
           [0, 1, 0, 0]  # No change for spatial in Minkowski

    # Try to raise already raised index
    assert AstroEquations.Physics.GeneralRelativity.raise_index([1, 0, 0, 0], metric, [:upper, :lower, :lower, :lower], 0) ==
           [1, 0, 0, 0]  # Shouldn't change
  end

  test "transform_tensor for vectors" do
    # Identity transformation
    jacobian = [
      [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
    ]

    assert AstroEquations.Physics.GeneralRelativity.transform_tensor([1, 2, 3, 4], jacobian, jacobian) == [1, 2, 3, 4]

    # Simple scaling transformation
    scale_jacobian = [
      [2, 0, 0, 0],
      [0, 2, 0, 0],
      [0, 0, 2, 0],
      [0, 0, 0, 2]
    ]

    scale_inverse = [
      [0.5, 0, 0, 0],
      [0, 0.5, 0, 0],
      [0, 0, 0.5, 0],
      [0, 0, 0, 0.5]
    ]

    assert AstroEquations.Physics.GeneralRelativity.transform_tensor([1, 0, 0, 0], scale_jacobian, scale_inverse) == [2, 0, 0, 0]
  end

  test "inverse_metric calculation" do
    minkowski = [
      [-1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1]
    ]

    assert AstroEquations.Physics.GeneralRelativity.inverse_metric(minkowski) == minkowski

    schwarzschild = [
      [-2.0, 0, 0, 0],
      [0, 0.5, 0, 0],
      [0, 0, 1/9, 0],
      [0, 0, 0, 1/9]
    ]

    original = [
      [-0.5, 0, 0, 0],
      [0, 2.0, 0, 0],
      [0, 0, 9, 0],
      [0, 0, 0, 9]
    ]

    assert AstroEquations.Physics.GeneralRelativity.inverse_metric(original) == schwarzschild
  end

  test "index operations with Schwarzschild metric", %{schw: metric} do
    # Test lowering time component
    assert AstroEquations.Physics.GeneralRelativity.lower_index([1, 0, 0, 0], metric, [:upper, :upper, :upper, :upper], 0) ==
           [-0.5, 0, 0, 0]

    # Test raising radial component
    inverse_metric = AstroEquations.Physics.GeneralRelativity.inverse_metric(metric)
    assert AstroEquations.Physics.GeneralRelativity.raise_index([0, 2.0, 0, 0], inverse_metric, [:lower, :lower, :lower, :lower], 1) ==
           [0, 1, 0, 0]
  end
end
