defmodule AstroEquations.Physics.ThermodynamicsTest do
  use ExUnit.Case
  doctest AstroEquations.Physics.Thermodynamics

  test "ideal gas law calculates pressure correctly" do
    result =
      AstroEquations.Physics.Thermodynamics.ideal_gas_law(
        p: nil,
        v: 0.0224,
        n: 6.022e23,
        t: 273.15
      )

    assert_in_delta result.p, 101_325, 100
  end

  test "heat energy calculation" do
    assert AstroEquations.Physics.Thermodynamics.heat_energy(1.0, 4186, 1.0) == 4186.0
  end

  test "heat capacity calculation" do
    assert AstroEquations.Physics.Thermodynamics.heat_capacity(4186, 1.0) == 4186.0
  end

  test "specific heat capacity calculation" do
    assert AstroEquations.Physics.Thermodynamics.specific_heat_capacity(4186, 1.0) == 4186.0
  end

  test "microstates calculation" do
    assert AstroEquations.Physics.Thermodynamics.microstates(3, 3) == 10
  end

  test "entropy calculation" do
    assert_in_delta AstroEquations.Physics.Thermodynamics.entropy(10), 3.179e-23, 1.0e-25
  end

  test "photon energy calculation" do
    assert AstroEquations.Physics.Thermodynamics.photon_energy(1.0e15) == 6.62607015e-19
  end

  test "Wien's displacement law" do
    assert_in_delta AstroEquations.Physics.Thermodynamics.wiens_displacement(5000),
                    5.7955458e-7,
                    1.0e-10
  end

  test "Stefan-Boltzmann law" do
    assert_in_delta AstroEquations.Physics.Thermodynamics.stefan_boltzmann(5000), 3.54375e7, 1.0e3
  end

  test "Planck's law (wavelength form)" do
    result = AstroEquations.Physics.Thermodynamics.planck_wavelength(5.0e-7, 5000)
    assert_in_delta result, 1.3714e13, 1.0e9
  end

  test "Planck's law (frequency form)" do
    result = AstroEquations.Physics.Thermodynamics.planck_frequency(1.0e14, 5000)
    assert_in_delta result, 1.1144e-16, 1.0e-18
  end
end
