defmodule AstroEquations.Physics.Thermodynamics do
  @moduledoc """
  A module for calculating various thermodynamics formulas and relationships.

  This module provides functions for ideal gases, microstates, entropy,
  and black body radiation calculations.
  """

  @doc """
  Calculates the ideal gas law: pV = NkBT

  ## Parameters
    - p: Pressure (Pa)
    - v: Volume (m³)
    - n: Number of particles (unitless)
    - k_b: Boltzmann constant (default: 1.380649e-23 J/K)
    - t: Temperature (K)

  ## Returns
    The missing value in the equation (pressure, volume, number of particles, or temperature)
    based on which parameter is set to nil.

  ## Examples
      iex> Thermodynamics.ideal_gas_law(p: 101325, v: 0.0224, n: 6.022e23, t: 273.15)
      %{k_b: 1.380649e-23}  # Verifies the Boltzmann constant
  """
  def ideal_gas_law(p: nil, v: v, n: n, k_b: k_b, t: t), do: %{p: n * k_b * t / v}
  def ideal_gas_law(p: p, v: nil, n: n, k_b: k_b, t: t), do: %{v: n * k_b * t / p}
  def ideal_gas_law(p: p, v: v, n: nil, k_b: k_b, t: t), do: %{n: p * v / (k_b * t)}
  def ideal_gas_law(p: p, v: v, n: n, k_b: k_b, t: nil), do: %{t: p * v / (n * k_b)}
  def ideal_gas_law(p: p, v: v, n: n, k_b: nil, t: t), do: %{k_b: p * v / (n * t)}
  def ideal_gas_law(p: p, v: v, n: n, t: t), do: ideal_gas_law(p: p, v: v, n: n, k_b: 1.380649e-23, t: t)

  @doc """
  Calculates heat/thermal energy: Q = mcΔT

  ## Parameters
    - m: Mass (kg)
    - c: Specific heat capacity (J/(kg·K))
    - delta_t: Temperature change (K)

  ## Returns
    The thermal energy in joules

  ## Examples
      iex> Thermodynamics.heat_energy(1.0, 4186, 1.0)
      4186.0
  """
  def heat_energy(m, c, delta_t), do: m * c * delta_t

  @doc """
  Calculates heat capacity: C = dQ/dT

  ## Parameters
    - dq: Change in heat (J)
    - dt: Change in temperature (K)

  ## Returns
    The heat capacity in J/K

  ## Examples
      iex> Thermodynamics.heat_capacity(4186, 1.0)
      4186.0
  """
  def heat_capacity(dq, dt), do: dq / dt

  @doc """
  Calculates specific heat capacity: c = C/m

  ## Parameters
    - c_heat: Heat capacity (J/K)
    - m: Mass (kg)

  ## Returns
    The specific heat capacity in J/(kg·K)

  ## Examples
      iex> Thermodynamics.specific_heat_capacity(4186, 1.0)
      4186.0
  """
  def specific_heat_capacity(c_heat, m), do: c_heat / m

  @doc """
  Calculates the number of microstates: Ω = (q + N - 1)! / (q!(N-1)!)

  ## Parameters
    - q: Number of energy quanta
    - n: Number of particles

  ## Returns
    The number of microstates

  ## Examples
      iex> Thermodynamics.microstates(3, 3)
      10
  """
  def microstates(q, n) do
    numerator = factorial(q + n - 1)
    denominator = factorial(q) * factorial(n - 1)
    div(numerator, denominator)
  end

  @doc """
  Calculates entropy: S = kB ln Ω

  ## Parameters
    - omega: Number of microstates
    - k_b: Boltzmann constant (default: 1.380649e-23 J/K)

  ## Returns
    The entropy in J/K

  ## Examples
      iex> Thermodynamics.entropy(10)
      3.179134965721792e-23
  """
  def entropy(omega, k_b \\ 1.380649e-23), do: k_b * :math.log(omega)

  @doc """
  Calculates photon energy: E = hf

  ## Parameters
    - f: Frequency (Hz)
    - h: Planck's constant (default: 6.62607015e-34 J·s)

  ## Returns
    The photon energy in joules

  ## Examples
      iex> Thermodynamics.photon_energy(1.0e15)
      6.62607015e-19
  """
  def photon_energy(f, h \\ 6.62607015e-34), do: h * f

  @doc """
  Calculates Wien's displacement law: λ_max = b/T

  ## Parameters
    - t: Temperature (K)
    - b: Wien's displacement constant (default: 2.8977729e-3 m·K)

  ## Returns
    The wavelength of maximum emission in meters

  ## Examples
      iex> Thermodynamics.wiens_displacement(5000)
      5.7955458e-7
  """
  def wiens_displacement(t, b \\ 2.8977729e-3), do: b / t

  @doc """
  Calculates Stefan-Boltzmann law: I = σT⁴

  ## Parameters
    - t: Temperature (K)
    - sigma: Stefan-Boltzmann constant (default: 5.670374419e-8 W/(m²·K⁴))

  ## Returns
    The radiant exitance in W/m²

  ## Examples
      iex> Thermodynamics.stefan_boltzmann(5000)
      3.54375e7
  """
  def stefan_boltzmann(t, sigma \\ 5.670374419e-8), do: sigma * :math.pow(t, 4)

  @doc """
  Calculates Planck's law for spectral radiance (wavelength form)

  ## Parameters
    - lambda: Wavelength (m)
    - t: Temperature (K)
    - h: Planck's constant (default: 6.62607015e-34 J·s)
    - c: Speed of light (default: 299792458 m/s)
    - k_b: Boltzmann constant (default: 1.380649e-23 J/K)

  ## Returns
    The spectral radiance in W/(sr·m³)

  ## Examples
      iex> Thermodynamics.planck_wavelength(5.0e-7, 5000) |> Float.round(8)
      1.3714e13
  """
  def planck_wavelength(lambda, t, h \\ 6.62607015e-34, c \\ 299792458, k_b \\ 1.380649e-23) do
    numerator = 2 * h * :math.pow(c, 2)
    denominator = :math.pow(lambda, 5)
    exponent = h * c / (lambda * k_b * t)
    (numerator / denominator) * 1 / (:math.exp(exponent) - 1)
  end

  @doc """
  Calculates Planck's law for spectral radiance (frequency form)

  ## Parameters
    - nu: Frequency (Hz)
    - t: Temperature (K)
    - h: Planck's constant (default: 6.62607015e-34 J·s)
    - c: Speed of light (default: 299792458 m/s)
    - k_b: Boltzmann constant (default: 1.380649e-23 J/K)

  ## Returns
    The spectral radiance in W/(sr·m²·Hz)

  ## Examples
      iex> Thermodynamics.planck_frequency(1.0e14, 5000) |> Float.round(8)
      1.1144e-16
  """
  def planck_frequency(nu, t, h \\ 6.62607015e-34, c \\ 299792458, k_b \\ 1.380649e-23) do
    numerator = 2 * h * :math.pow(nu, 3)
    denominator = :math.pow(c, 2)
    exponent = h * nu / (k_b * t)
    (numerator / denominator) * 1 / (:math.exp(exponent) - 1)
  end

  # Helper function for factorial calculation
  defp factorial(0), do: 1
  defp factorial(n) when n > 0, do: n * factorial(n - 1)
end
