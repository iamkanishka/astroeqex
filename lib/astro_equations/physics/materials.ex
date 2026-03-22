defmodule AstroEquations.Physics.Materials do
  @moduledoc """
  Material properties relevant to physics and astrophysics.

  Covers:
  - Density (discrete and continuous)
  - Stress and strain (Young's modulus, shear, bulk modulus)
  - Poisson ratio
  - Thermal expansion
  - Viscosity and Reynolds number
  - Mach number and compressibility
  - Electrical conductivity and resistivity
  - Thermal conductivity (Fourier's law) and thermal diffusivity
  - Speed of sound in solids, liquids, and gases
  - Elastic energy density
  - Equation of state (polytropic)
  - Stokes number (particle dynamics)
  """
  use AstroEquations.Guards

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Density in kg/m³. Must be positive."
  @type density :: float()

  @typedoc "Stress in pascals (Pa)."
  @type stress :: float()

  @typedoc "Strain (dimensionless, ΔL/L)."
  @type strain :: float()

  @typedoc "Elastic modulus in pascals (Pa). Must be positive."
  @type modulus :: float()

  @typedoc "Dynamic viscosity in Pa·s. Must be positive."
  @type viscosity :: float()

  # ---------------------------------------------------------------------------
  # Density
  # ---------------------------------------------------------------------------

  @doc """
  Density: ρ = m / V

  ## Parameters
    - mass:   Total mass (kg)
    - volume: Total volume (m³)

  ## Examples
      iex> Materials.density(10, 2)
      5.0
  """

  @spec density(number, number) :: float
  def density(mass, volume) when is_positive(mass) and is_positive(volume), do: mass / volume

  @doc """
  Continuous density as the derivative dm/dV (numerical central difference).

  ## Parameters
    - mass_function: m(V) function
    - volume:        V at which to evaluate (m³)
    - delta:         Step size for finite difference (default: 1.0e-6)

  ## Examples
      iex> Materials.continuous_density(fn v -> 2.5 * v end, 2.0) |> Float.round(4)
      2.5
  """

  @spec continuous_density((number -> number), number, number) :: float
  def continuous_density(mass_function, volume, delta \\ 1.0e-6) do
    (mass_function.(volume + delta) - mass_function.(volume - delta)) / (2 * delta)
  end

  @doc """
  Number density: n = N / V

  ## Examples
      iex> Materials.number_density(1.0e23, 0.001)
      1.0e26
  """

  @spec number_density(number, number) :: float
  def number_density(count, volume), do: count / volume

  # ---------------------------------------------------------------------------
  # Stress & Strain
  # ---------------------------------------------------------------------------

  @doc """
  Normal stress: σ = F / A

  ## Examples
      iex> Materials.normal_stress(1000, 0.01) |> Float.round(2)
      100000.0
  """

  @spec normal_stress(number, number) :: float
  def normal_stress(force, area), do: force / area

  @doc """
  Normal (engineering) strain: ε = ΔL / L₀

  ## Examples
      iex> Materials.normal_strain(0.001, 1.0)
      0.001
  """

  @spec normal_strain(number, number) :: float
  def normal_strain(delta_l, l0), do: delta_l / l0

  @doc """
  Young's modulus: E = σ / ε

  ## Examples
      iex> Materials.youngs_modulus(200.0e6, 0.001)
      2.0e11
  """

  @spec youngs_modulus(number, number) :: float
  def youngs_modulus(stress, strain) when is_number(stress) and is_positive(strain),
    do: stress / strain

  @doc """
  Extension from Young's modulus: ΔL = F L₀ / (E A)

  ## Examples
      iex> Materials.extension(1000, 1.0, 200.0e9, 1.0e-4) |> Float.round(8)
      5.0e-5
  """

  @spec extension(number, number, number, number) :: float
  def extension(force, original_length, youngs_modulus, area) do
    force * original_length / (youngs_modulus * area)
  end

  @doc """
  Shear stress: τ = F_shear / A

  ## Examples
      iex> Materials.shear_stress(500, 0.01)
      50000.0
  """

  @spec shear_stress(number, number) :: float
  def shear_stress(shear_force, area), do: shear_force / area

  @doc """
  Shear modulus (modulus of rigidity): G = τ / γ

  ## Examples
      iex> Materials.shear_modulus(80.0e9, 0.001) > 0
      true
  """

  @spec shear_modulus(number, number) :: float
  def shear_modulus(shear_stress, shear_strain), do: shear_stress / shear_strain

  @doc """
  Bulk modulus: K = -V (dP/dV) ≈ -ΔP / (ΔV/V)

  ## Examples
      iex> Materials.bulk_modulus(1.0e6, -1.0e-6, 1.0e-3) > 0
      true
  """

  @spec bulk_modulus(number, number, number) :: float
  def bulk_modulus(delta_p, delta_v, original_volume) do
    -delta_p / (delta_v / original_volume)
  end

  @doc """
  Poisson's ratio: ν = -ε_lateral / ε_axial

  Ratio of lateral contraction to axial extension under uniaxial stress.
  Typical values: metals ≈ 0.25–0.35, rubber ≈ 0.5, cork ≈ 0.

  ## Parameters
    - lateral_strain: ε_lateral (negative for contraction, positive sign here)
    - axial_strain:   ε_axial (positive for extension)

  ## Examples
      iex> Materials.poisson_ratio(0.001, 0.003) |> Float.round(4)
      0.3333
  """

  @spec poisson_ratio(number, number) :: float
  def poisson_ratio(lateral_strain, axial_strain) do
    lateral_strain / axial_strain
  end

  @doc """
  Compressibility (inverse bulk modulus): κ = 1 / K  (Pa⁻¹)

  ## Examples
      iex> Materials.compressibility(2.0e9) |> Float.round(12)
      5.0e-10
  """

  @spec compressibility(number) :: float
  def compressibility(bulk_modulus), do: 1.0 / bulk_modulus

  @doc """
  Elastic energy density: u = ½ E ε²

  ## Examples
      iex> Materials.elastic_energy_density(200.0e9, 0.001) > 0
      true
  """

  @spec elastic_energy_density(number, number) :: float
  def elastic_energy_density(youngs_modulus, strain) do
    0.5 * youngs_modulus * strain * strain
  end

  # ---------------------------------------------------------------------------
  # Thermal Properties
  # ---------------------------------------------------------------------------

  @doc """
  Linear thermal expansion: ΔL = α L₀ ΔT

  ## Examples
      iex> Materials.linear_expansion(12.0e-6, 1.0, 100) |> Float.round(6)
      0.0012
  """

  @spec linear_expansion(number, number, number) :: number()
  def linear_expansion(alpha, original_length, delta_t) do
    alpha * original_length * delta_t
  end

  @doc """
  Volumetric thermal expansion: ΔV = β V₀ ΔT  (β ≈ 3α for isotropic solids)

  ## Examples
      iex> Materials.volumetric_expansion(36.0e-6, 1.0, 100) |> Float.round(6)
      0.0036
  """

  @spec volumetric_expansion(number, number, number) :: number()
  def volumetric_expansion(beta, original_volume, delta_t) do
    beta * original_volume * delta_t
  end

  @doc """
  Fourier's Law of heat conduction: Q/t = k A ΔT / d

  ## Returns
    Heat flow rate (W)

  ## Examples
      iex> Materials.heat_conduction(401, 0.01, 100, 0.01) |> Float.round(1)
      40100.0
  """

  @spec heat_conduction(number, number, number, number) :: float
  def heat_conduction(k, area, delta_t, thickness) do
    k * area * delta_t / thickness
  end

  @doc """
  Thermal resistance: R_th = d / (k A)

  ## Examples
      iex> Materials.thermal_resistance(0.01, 401, 0.01) > 0
      true
  """

  @spec thermal_resistance(number, number, number) :: float
  def thermal_resistance(thickness, k, area), do: thickness / (k * area)

  @doc """
  Thermal diffusivity: α = k / (ρ c_p)

  Governs how quickly temperature disturbances spread through a material.

  ## Parameters
    - k:           Thermal conductivity (W/(m·K))
    - density:     ρ (kg/m³)
    - specific_cp: Specific heat at constant pressure (J/(kg·K))

  ## Examples
      iex> Materials.thermal_diffusivity(401, 8960, 385) > 0
      true
  """

  @spec thermal_diffusivity(number, number, number) :: float
  def thermal_diffusivity(k, density, specific_cp) when is_positive(density) do
    k / (density * specific_cp)
  end

  # ---------------------------------------------------------------------------
  # Fluid Properties
  # ---------------------------------------------------------------------------

  @doc """
  Dynamic viscosity shear stress: τ = η (dv/dy)

  ## Examples
      iex> Materials.viscous_shear_stress(1.0e-3, 100)
      0.1
  """

  @spec viscous_shear_stress(number, number) :: number()
  def viscous_shear_stress(eta, dv_dy), do: eta * dv_dy

  @doc """
  Kinematic viscosity: ν = η / ρ

  ## Examples
      iex> Materials.kinematic_viscosity(1.81e-5, 1.225) > 0
      true
  """

  @spec kinematic_viscosity(number, number) :: float
  def kinematic_viscosity(eta, rho), do: eta / rho

  @doc """
  Reynolds number: Re = ρ v L / η

  ## Examples
      iex> Materials.reynolds_number(1.225, 30, 0.1, 1.81e-5) > 4000
      true
  """

  @spec reynolds_number(number, number, number, number) :: float
  def reynolds_number(rho, v, l, eta)
      when is_positive(rho) and is_positive(l) and is_positive(eta),
      do: rho * v * l / eta

  @doc """
  Mach number: Ma = v / v_sound

  ## Parameters
    - velocity:    Flow speed (m/s)
    - sound_speed: Speed of sound in the medium (m/s)

  ## Examples
      iex> Materials.mach_number(343.0, 343.0) |> Float.round(4)
      1.0
  """

  @spec mach_number(number, number) :: float
  def mach_number(velocity, sound_speed), do: velocity / sound_speed

  @doc """
  Stokes number: Stk = τ_p v₀ / L

  Ratio of particle relaxation time to the characteristic flow time scale.
  Stk ≪ 1 → particles follow fluid; Stk ≫ 1 → particles resist fluid.

  ## Parameters
    - relaxation_time: τ_p = ρ_p d² / (18 η) (s)
    - flow_velocity:   v₀ (m/s)
    - length_scale:    L (m)

  ## Examples
      iex> Materials.stokes_number(1.0e-3, 10.0, 0.1) |> Float.round(4)
      0.1
  """

  @spec stokes_number(number, number, number) :: float
  def stokes_number(relaxation_time, flow_velocity, length_scale) do
    relaxation_time * flow_velocity / length_scale
  end

  @doc """
  Particle relaxation time (Stokes drag): τ_p = ρ_p d² / (18 η)

  ## Parameters
    - particle_density: ρ_p (kg/m³)
    - diameter:         d (m)
    - viscosity:        η (Pa·s)

  ## Examples
      iex> Materials.particle_relaxation_time(1000.0, 1.0e-5, 1.8e-5) > 0
      true
  """

  @spec particle_relaxation_time(number, number, number) :: float
  def particle_relaxation_time(particle_density, diameter, viscosity)
      when is_positive(diameter) do
    particle_density * diameter * diameter / (18 * viscosity)
  end

  # ---------------------------------------------------------------------------
  # Sound Speed
  # ---------------------------------------------------------------------------

  @doc """
  Speed of sound in a solid: v_s = √(E / ρ)

  ## Examples
      iex> Materials.speed_of_sound_solid(200.0e9, 7900.0) |> Float.round(0)
      5032.0
  """

  @spec speed_of_sound_solid(number, number) :: float
  def speed_of_sound_solid(youngs_modulus, density)
      when is_positive(youngs_modulus) and is_positive(density) do
    :math.sqrt(youngs_modulus / density)
  end

  @doc """
  Speed of sound in a liquid: v_s = √(K / ρ)

  Uses the bulk modulus K (appropriate for liquids that cannot sustain shear).

  ## Parameters
    - bulk_modulus: K (Pa)
    - density:      ρ (kg/m³)

  ## Examples
      iex> Materials.speed_of_sound_liquid(2.2e9, 1000.0) |> Float.round(0)
      1483.0
  """

  @spec speed_of_sound_liquid(number, number) :: float
  def speed_of_sound_liquid(bulk_modulus, density) when is_positive(density) do
    :math.sqrt(bulk_modulus / density)
  end

  @doc """
  Speed of sound in an ideal gas: v_s = √(γ R T / M)

  ## Examples
      iex> Materials.speed_of_sound_gas(1.4, 293.0, 0.029) |> Float.round(0)
      343.0
  """

  @spec speed_of_sound_gas(number, number, number, number) :: float
  def speed_of_sound_gas(gamma, temperature, molar_mass, r \\ 8.31446)
      when is_positive(gamma) and is_positive(temperature) and is_positive(molar_mass) do
    :math.sqrt(gamma * r * temperature / molar_mass)
  end

  # ---------------------------------------------------------------------------
  # Equation of State
  # ---------------------------------------------------------------------------

  @doc """
  Polytropic equation of state: P = K ρ^Γ

  Common in stellar interior models; Γ = 5/3 for adiabatic monatomic gas,
  Γ = 4/3 for a radiation-dominated or ultra-relativistic gas.

  ## Examples
      iex> Materials.polytropic_pressure(1.0e5, 1000.0, 1.4) > 0
      true
  """

  @spec polytropic_pressure(number, number, number) :: float

  def polytropic_pressure(k, rho, gamma) when is_positive(k) and is_positive(rho),
    do: k * :math.pow(rho, gamma)
end
