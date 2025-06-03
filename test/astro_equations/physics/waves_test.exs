defmodule AstroEquations.Physics.WavesTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.Waves

  describe "work/2" do
    test "calculates work for parallel vectors" do
      assert AstroEquations.Physics.Waves.work([2, 0, 0], [3, 0, 0]) == 6.0
    end

    test "calculates work for non-parallel vectors" do
      assert AstroEquations.Physics.Waves.work([1, 2, 3], [4, 5, 6]) == 32.0
    end
  end
end
