defmodule AstroEquations.Physics.Forces do
  @doc """
  Calculates buoyant force using Archimedes' Principle.

  ## Parameters
    - displaced_mass: Mass of displaced fluid (kg)
    - gravity: Acceleration due to gravity (m/s²), defaults to 9.81

  ## Examples
      iex> Physics.buoyancy(5)
      49.05
  """
  @spec buoyancy(number, number) :: float
  def buoyancy(displaced_mass, gravity \\ 9.81) do
    displaced_mass * gravity
  end

  @doc """
  Calculates buoyant force using fluid density and displaced volume.

  ## Parameters
    - density: Density of the fluid (kg/m³)
    - volume: Volume displaced (m³)
    - gravity: Acceleration due to gravity (m/s²), defaults to 9.81

  ## Examples
      iex> Physics.buoyancy_from_density(1000, 0.005)
      49.05
  """
  @spec buoyancy_from_density(number, number, number) :: float
  def buoyancy_from_density(density, volume, gravity \\ 9.81) do
    density * volume * gravity
  end

  @doc """
  Calculates kinetic friction force.

  ## Parameters
    - coefficient: Coefficient of kinetic friction (μₖ)
    - normal_force: Normal force (N)

  ## Examples
      iex> Physics.kinetic_friction(0.3, 10)
      3.0
  """
  @spec kinetic_friction(number, number) :: float
  def kinetic_friction(coefficient, normal_force) do
    coefficient * normal_force
  end

  @doc """
  Calculates maximum static friction force.

  ## Parameters
    - coefficient: Coefficient of static friction (μₛ)
    - normal_force: Normal force (N)

  ## Examples
      iex> Physics.static_friction(0.4, 10)
      4.0
  """
  @spec static_friction(number, number) :: float
  def static_friction(coefficient, normal_force) do
    coefficient * normal_force
  end
end
