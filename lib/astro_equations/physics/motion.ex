defmodule AstroEquations.Physics.Motion do
  @moduledoc """
  This module contains functions related to motion physics concepts including:
  - Velocity
  - Acceleration
  - Newton's Laws
  - Momentum
  - Centripetal Force
  - Kinetic Energy
  """

  @doc """
  Calculates velocity given displacement and time.

  ## Parameters
    - displacement: change in position (Δx) in meters
    - time: time interval (Δt) in seconds

  ## Returns
    Velocity in m/s

  ## Examples
      iex> Physics.Motion.velocity(10, 2)
      5.0
  """
  def velocity(displacement, time) do
    displacement / time
  end

  @doc """
  Calculates acceleration given change in velocity and time.

  ## Parameters
    - delta_v: change in velocity (Δv) in m/s
    - time: time interval (Δt) in seconds

  ## Returns
    Acceleration in m/s²

  ## Examples
      iex> Physics.Motion.acceleration(20, 5)
      4.0
  """
  def acceleration(delta_v, time) do
    delta_v / time
  end

  @doc """
  Newton's Second Law: calculates force given mass and acceleration.

  ## Parameters
    - mass: in kilograms
    - acceleration: in m/s²

  ## Returns
    Force in Newtons

  ## Examples
      iex> Physics.Motion.newtons_second_law(5, 2)
      10.0
  """
  def newtons_second_law(mass, acceleration) do
    mass * acceleration
  end

  @doc """
  Calculates momentum given mass and velocity.

  ## Parameters
    - mass: in kilograms
    - velocity: in m/s

  ## Returns
    Momentum in kg·m/s

  ## Examples
      iex> Physics.Motion.momentum(10, 5)
      50.0
  """
  def momentum(mass, velocity) do
    mass * velocity
  end

  @doc """
  Calculates change in momentum (impulse) given force and time interval.

  ## Parameters
    - force: in Newtons
    - time: time interval in seconds

  ## Returns
    Change in momentum in kg·m/s

  ## Examples
      iex> Physics.Motion.impulse(10, 5)
      50.0
  """
  def impulse(force, time) do
    force * time
  end

  @doc """
  Calculates centripetal force given mass, velocity, and radius.

  ## Parameters
    - mass: in kilograms
    - velocity: in m/s
    - radius: in meters

  ## Returns
    Centripetal force in Newtons

  ## Examples
      iex> Physics.Motion.centripetal_force(2, 5, 10)
      5.0
  """
  def centripetal_force(mass, velocity, radius) do
    mass * :math.pow(velocity, 2) / radius
  end

  @doc """
  Calculates kinetic energy given mass and velocity.

  ## Parameters
    - mass: in kilograms
    - velocity: in m/s

  ## Returns
    Kinetic energy in Joules

  ## Examples
      iex> Physics.Motion.kinetic_energy(4, 5)
      50.0
  """
  def kinetic_energy(mass, velocity) do
    0.5 * mass * :math.pow(velocity, 2)
  end
end
