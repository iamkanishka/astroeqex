defmodule AstroEquations.Physics.ElectromagnetismTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.Electromagnetism

  describe "gauss_law/2" do
    test "calculates electric flux for given enclosed charge" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law(1.0) == 1 / 8.8541878128e-12
    end
  end

  describe "gauss_law_differential/2" do
    test "calculates divergence of E field for given charge density" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law_differential(1.0) ==
               1 / 8.8541878128e-12
    end
  end

  describe "gauss_law_magnetism/0" do
    test "always returns zero" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law_magnetism() == 0
    end
  end

  describe "faraday_law/1" do
    test "calculates curl of E field from changing B field" do
      assert AstroEquations.Physics.Electromagnetism.faraday_law(2.5) == -2.5
    end
  end

  describe "ampere_law/4" do
    test "calculates curl of B field from current and changing E field" do
      current = 1.0
      dE_dt = 1.0
      expected = 1.25663706212e-6 * (current + 8.8541878128e-12 * dE_dt)
      assert AstroEquations.Physics.Electromagnetism.ampere_law(current, dE_dt) == expected
    end
  end

  describe "lorentz_force_point/4" do
    test "calculates force on point charge" do
      # Test case: E field in x, velocity in y, B field in z
      assert AstroEquations.Physics.Electromagnetism.lorentz_force_point(
               1.0,
               {1.0, 0.0, 0.0},
               {0.0, 1.0, 0.0},
               {0.0, 0.0, 1.0}
             ) == {1.0, 0.0, -1.0}
    end
  end

  describe "electric_field_point/3" do
    test "calculates electric field from point charge" do
      q = 1.0
      r = 1.0
      expected = 1 / (4 * :math.pi() * 8.8541878128e-12) * q / :math.pow(r, 2)
      assert AstroEquations.Physics.Electromagnetism.electric_field_point(q, r) == expected
    end
  end

  describe "drift_velocity/2" do
    test "calculates electron drift velocity" do
      assert AstroEquations.Physics.Electromagnetism.drift_velocity(0.001, 100) == 0.1
    end
  end

  describe "current_from_properties/5" do
    test "calculates current from charge carrier properties" do
      assert AstroEquations.Physics.Electromagnetism.current_from_properties(
               1.0e28,
               1.0e-6,
               1.6e-19,
               0.001,
               100
             ) ==
               1.6
    end
  end

  describe "electrical_power/2" do
    test "calculates power from current and voltage" do
      assert AstroEquations.Physics.Electromagnetism.electrical_power(2, 10) == 20
    end
  end

  describe "electrical_power_from_resistance/2" do
    test "calculates power from current and resistance" do
      assert AstroEquations.Physics.Electromagnetism.electrical_power_from_resistance(2, 5) == 20
    end
  end

  describe "ohms_law/2" do
    test "calculates voltage from current and resistance" do
      assert AstroEquations.Physics.Electromagnetism.ohms_law(2, 5) == 10
    end
  end

  describe "resistance/3" do
    test "calculates resistance from resistivity" do
      assert AstroEquations.Physics.Electromagnetism.resistance(1.68e-8, 1, 1.0e-6) == 0.0168
    end
  end

  describe "series_resistance/1" do
    test "calculates total series resistance" do
      assert AstroEquations.Physics.Electromagnetism.series_resistance([10, 20, 30]) == 60
    end
  end

  describe "parallel_resistance/1" do
    test "calculates total parallel resistance" do
      assert AstroEquations.Physics.Electromagnetism.parallel_resistance([10, 10]) == 5.0
    end
  end

  describe "capacitance/2" do
    test "calculates capacitance from charge and voltage" do
      assert AstroEquations.Physics.Electromagnetism.capacitance(1.0e-6, 10) == 1.0e-7
    end
  end

  describe "capacitance_from_geometry/3" do
    test "calculates capacitance from physical properties" do
      assert AstroEquations.Physics.Electromagnetism.capacitance_from_geometry(
               8.854e-12,
               0.01,
               1.0e-3
             ) ==
               8.854e-11
    end
  end

  describe "capacitor_energy/2" do
    test "calculates energy stored in capacitor" do
      assert AstroEquations.Physics.Electromagnetism.capacitor_energy(1.0e-6, 10) == 5.0e-5
    end
  end

  describe "capacitor_field/3" do
    test "calculates electric field in capacitor" do
      charge = 1.0e-6
      area = 0.01
      expected = charge / (8.8541878128e-12 * area)
      assert AstroEquations.Physics.Electromagnetism.capacitor_field(charge, area) == expected
    end
  end

  describe "biot_savart/4" do
    test "calculates magnetic field from current element" do
      # Current along z-axis, field point at x-axis
      result = Physics.Electromagnetism.biot_savart(1.0, {0.0, 0.0, 1.0e-3}, {1.0, 0.0, 0.0}, 1.0)
      assert elem(result, 1) == 1.0e-10
      assert elem(result, 0) == 0.0
      assert elem(result, 2) == 0.0
    end
  end

  describe "moving_charge_field/4" do
    test "calculates magnetic field from moving charge" do
      # Charge moving along y-axis, field point at x-axis
      q = 1.6e-19
      expected_z = -q / :math.pow(1.0, 2) * 1.25663706212e-6 / (4 * :math.pi())

      result =
        AstroEquations.Physics.Electromagnetism.moving_charge_field(
          q,
          {0.0, 1.0e6, 0.0},
          {1.0, 0.0, 0.0},
          1.0
        )

      assert elem(result, 2) == expected_z
    end
  end

  describe "wire_magnetic_field/3" do
    test "calculates magnetic field around wire" do
      assert AstroEquations.Physics.Electromagnetism.wire_magnetic_field(1.0, 0.1) == 2.0e-6
    end
  end

  describe "inductor_emf/2" do
    test "calculates induced EMF in inductor" do
      assert AstroEquations.Physics.Electromagnetism.inductor_emf(0.1, 10.0) == -1.0
    end
  end

  describe "inductor_energy/2" do
    test "calculates energy stored in inductor" do
      assert AstroEquations.Physics.Electromagnetism.inductor_energy(0.1, 2.0) == 0.2
    end
  end

  describe "vector_potential/3" do
    test "calculates magnetic vector potential" do
      assert AstroEquations.Physics.Electromagnetism.vector_potential({1.0, 0.0, 0.0}, 1.0) ==
               {1.25663706212e-7, 0.0, 0.0}
    end
  end

  describe "gauss_law_polarization/2" do
    test "calculates bound charge for perpendicular vectors" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law_polarization([1, 0, 0], [1, 0, 0]) ==
               -1.0

      assert AstroEquations.Physics.Electromagnetism.gauss_law_polarization([0, 2, 0], [0, 3, 0]) ==
               -6.0
    end

    test "returns zero for orthogonal vectors" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law_polarization([1, 0, 0], [0, 1, 0]) ==
               0.0
    end
  end

  describe "gauss_law_displacement/2" do
    test "calculates free charge for perpendicular vectors" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law_displacement([1, 0, 0], [1, 0, 0]) ==
               1.0

      assert AstroEquations.Physics.Electromagnetism.gauss_law_displacement([0, 2, 0], [0, 3, 0]) ==
               6.0
    end

    test "returns zero for orthogonal vectors" do
      assert AstroEquations.Physics.Electromagnetism.gauss_law_displacement([1, 0, 0], [0, 1, 0]) ==
               0.0
    end
  end

  describe "relative_permittivity/2" do
    test "calculates relative permittivity" do
      assert AstroEquations.Physics.Electromagnetism.relative_permittivity(1.77e-11, 8.854e-12) ==
               2.0
    end
  end

  describe "electric_susceptibility/1" do
    test "calculates electric susceptibility" do
      assert AstroEquations.Physics.Electromagnetism.electric_susceptibility(2.0) == 1.0
      assert AstroEquations.Physics.Electromagnetism.electric_susceptibility(1.0) == 0.0
    end
  end

  describe "absolute_permittivity/2" do
    test "calculates absolute permittivity" do
      assert_in_delta AstroEquations.Physics.Electromagnetism.absolute_permittivity(
                        2.0,
                        8.854e-12
                      ),
                      1.7708e-11,
                      1.0e-15
    end
  end

  describe "polarization/4" do
    test "calculates from susceptibility and electric field" do
      assert Maxwell.Materials.polarization(1.0, [1, 0, 0]) ==
               [8.8541878128e-12, 0.0, 0.0]
    end

    test "calculates from dipole moment density" do
      assert Maxwell.Materials.polarization(0, [0, 0, 0], 1.0e28, [1.0e-29, 0.0, 0.0]) ==
               [1.0e-1, 0.0, 0.0]
    end
  end

  describe "surface_bound_charge/2" do
    test "calculates surface bound charge" do
      assert Maxwell.Materials.surface_bound_charge([1, 2, 3], [0, 1, 0]) == 2.0
      assert Maxwell.Materials.surface_bound_charge([1, 1, 1], [1, 1, 1]) == 3.0
    end
  end

  describe "volume_bound_charge/2" do
    test "approximates volume bound charge" do
      assert Maxwell.Materials.volume_bound_charge([[1, 0, 0], [1.1, 0, 0]], 1.0e-3) == -100.0
    end
  end

  describe "electric_displacement/3" do
    test "calculates from permittivity" do
      assert Maxwell.Materials.electric_displacement(2.0, [1, 0, 0]) == [2.0, 0.0, 0.0]
    end

    test "calculates with polarization" do
      assert Maxwell.Materials.electric_displacement(1.0, [1, 0, 0], [0.5, 0, 0]) == [
               1.5,
               0.0,
               0.0
             ]
    end
  end

  describe "magnetic_field_strength/2" do
    test "calculates without magnetization" do
      h = Maxwell.Materials.magnetic_field_strength([1.0, 0, 0])
      assert_in_delta h |> hd(), 795_774.715, 0.001
    end

    test "calculates with magnetization" do
      h = Maxwell.Materials.magnetic_field_strength([1.0, 0, 0], [1000, 0, 0])
      assert_in_delta h |> hd(), 794_774.715, 0.001
    end
  end

  describe "bound_surface_current/2" do
    test "calculates surface current" do
      assert Maxwell.Materials.bound_surface_current([0, 0, 1], [1, 0, 0]) == [0.0, 1.0, 0.0]
    end
  end
end
