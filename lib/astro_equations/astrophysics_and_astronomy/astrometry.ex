defmodule AstroEquations.AstrophysicsAndAstronomy.Astrometry do
  @moduledoc """
  A collection of astronomical formulas for redshift, magnitude, flux, and metallicity calculations.

  This module provides functions to perform common astrometric calculations including:
  - Redshift calculations
  - Apparent and absolute magnitude conversions
  - Flux-magnitude relationships
  - Color index calculations
  - Metallicity determinations
  """

  @doc """
  Calculates the redshift (z) from observed and emitted wavelengths.

  ## Parameters
    - λ_obs: Observed wavelength (in any consistent units)
    - λ_emit: Emitted wavelength (same units as λ_obs)

  ## Returns
    The redshift value (z)

  ## Examples
      iex> Astrometry.redshift(700, 656)
      0.06707317073170732
  """
  @spec redshift(number, number) :: float
  def redshift(λ_obs, λ_emit) do
    (λ_obs - λ_emit) / λ_emit
  end

  @doc """
  Calculates the redshift (z) using the ratio method (1 + z = λ_obs/λ_emit).

  ## Parameters
    - λ_obs: Observed wavelength (in any consistent units)
    - λ_emit: Emitted wavelength (same units as λ_obs)

  ## Returns
    The redshift value (z)

  ## Examples
      iex> Astrometry.redshift_ratio(700, 656)
      0.06707317073170732
  """
  @spec redshift_ratio(number, number) :: float
  def redshift_ratio(λ_obs, λ_emit) do
    λ_obs / λ_emit - 1
  end

  @doc """
  Calculates apparent magnitude difference from flux ratio.

  ## Parameters
    - m: Apparent magnitude of object
    - m0: Reference apparent magnitude
    - F: Flux of object
    - F0: Reference flux

  ## Returns
    The difference in apparent magnitudes (m - m0)

  ## Examples
      iex> Astrometry.apparent_magnitude_diff(12, 10, 1.0, 6.31)
      -2.500000000000001
  """
  @spec apparent_magnitude_diff(number, number, number, number) :: float
  def apparent_magnitude_diff(m, m0, F, F0) do
    -2.5 * :math.log10(F / F0)
  end

  @doc """
  Calculates absolute magnitude from apparent magnitude and distance.

  ## Parameters
    - m: Apparent magnitude
    - d: Distance in parsecs

  ## Returns
    The absolute magnitude (M)

  ## Examples
      iex> Astrometry.absolute_magnitude(5, 100)
      0.0
  """
  @spec absolute_magnitude(number, number) :: float
  def absolute_magnitude(m, d) do
    m - 5 * :math.log10(d / 10)
  end

  @doc """
  Calculates the flux ratio from relative magnitudes.

  ## Parameters
    - ma: Magnitude of object A
    - mb: Magnitude of object B

  ## Returns
    The flux ratio (Ia/Ib)

  ## Examples
      iex> Astrometry.flux_ratio_from_magnitudes(10, 15)
      100.0
  """
  @spec flux_ratio_from_magnitudes(number, number) :: float
  def flux_ratio_from_magnitudes(ma, mb) do
    10 ** (0.4 * (mb - ma))
  end

  @doc """
  Calculates flux from magnitude using a reference flux.

  ## Parameters
    - m: Magnitude
    - F0: Reference flux for zero magnitude

  ## Returns
    The flux (F)

  ## Examples
      iex> Astrometry.flux_from_magnitude(0, 1.0)
      1.0
      iex> Astrometry.flux_from_magnitude(1, 1.0) |> Float.round(4)
      0.3981
  """
  @spec flux_from_magnitude(number, number) :: float
  def flux_from_magnitude(m, F0) do
    F0 * 10 ** (-0.4 * m)
  end

  @doc """
  Calculates color index from fluxes in two different filters.

  ## Parameters
    - F_f1: Flux in filter 1
    - F_f2: Flux in filter 2

  ## Returns
    The color index

  ## Examples
      iex> Astrometry.color_index(1.0, 2.0) |> Float.round(4)
      0.7526
  """
  @spec color_index(number, number) :: float
  def color_index(F_f1, F_f2) do
    -2.5 * :math.log(F_f1 / F_f2)
  end

  @doc """
  Calculates metallicity [Fe/H] relative to solar.

  ## Parameters
    - fe_h_star: [Fe/H] ratio for the star
    - fe_h_solar: [Fe/H] ratio for the Sun (solar value)

  ## Returns
    The metallicity (Z)

  ## Examples
      iex> Astrometry.metallicity(-0.5, 0.0)
      -0.5
  """
  @spec metallicity(number, number) :: float
  def metallicity(fe_h_star, fe_h_solar) do
    :math.log10(fe_h_star / fe_h_solar)
  end
end
