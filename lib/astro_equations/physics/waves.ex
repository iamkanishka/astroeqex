defmodule AstroEquations.Physics.Waves do
  @moduledoc """
  A collection of fundamental physics formulas and calculations.

  This module provides functions for calculating wave number, wave velocity,
  wave properties, and related physics concepts.
  """

  @doc """
  Calculates wave number (k) from wavelength.

  ## Parameters
    - wavelength: Wavelength (m)

  ## Examples
      iex> Physics.wave_number(2)
      3.141592653589793
  """
  @spec wave_number(number) :: float
  def wave_number(wavelength) do
    2 * :math.pi() / wavelength
  end

  @doc """
  Calculates wave velocity from frequency and wavelength.

  ## Parameters
    - frequency: Frequency (Hz)
    - wavelength: Wavelength (m)

  ## Examples
      iex> Physics.wave_velocity(440, 0.78)
      343.2
  """
  @spec wave_velocity(number, number) :: float
  def wave_velocity(frequency, wavelength) do
    frequency * wavelength
  end

  @doc """
  Calculates angular frequency from period or frequency.

  ## Options
    - period: Period (s) - either period or frequency must be provided
    - frequency: Frequency (Hz)

  ## Examples
      iex> Physics.angular_frequency(period: 2)
      3.141592653589793

      iex> Physics.angular_frequency(frequency: 50)
      314.1592653589793
  """
  @spec angular_frequency(Keyword.t()) :: float
  def angular_frequency(opts) do
    cond do
      Keyword.has_key?(opts, :period) ->
        2 * :math.pi() / Keyword.get(opts, :period)

      Keyword.has_key?(opts, :frequency) ->
        2 * :math.pi() * Keyword.get(opts, :frequency)

      true ->
        raise ArgumentError, "Must provide either :period or :frequency"
    end
  end

  @doc """
  Generates a wave function at a specific time and position.

  ## Parameters
    - amplitude: Amplitude of the wave
    - angular_frequency: Angular frequency (ω)
    - wave_number: Wave number (k)
    - time: Time (s)
    - position: Position (m)
    - phase_constant: Phase constant (φ), defaults to 0

  ## Examples
      iex> Physics.wave_function(1, 2, 3, 4, 5)
      -0.5365729180004349
  """
  @spec wave_function(number, number, number, number, number, number) :: float
  def wave_function(
        amplitude,
        angular_frequency,
        wave_number,
        time,
        position,
        phase_constant \\ 0
      ) do
    amplitude * :math.sin(angular_frequency * time - wave_number * position + phase_constant)
  end
end
