defmodule AstroEquations.AstrophysicsAndAstronomy.Astrometry do
  @moduledoc """
  A collection of astronomical formulas for redshift, magnitude, flux, metallicity,
  parallax, proper motion, and distance calculations.

  Provides:
  - Redshift calculations (classical, ratio, relativistic)
  - Apparent and absolute magnitude conversions
  - Flux-magnitude relationships and zero-point conversions
  - Color index, color excess, and interstellar extinction
  - Metallicity determinations
  - Parallax and distance modulus
  - Proper motion and space velocity
  - Bolometric corrections
  - Angular diameter, luminosity, and angular diameter distances
  - Surface brightness
  """

  # arcseconds per radian
  @arcsec_per_rad 206_265.0

  # ---------------------------------------------------------------------------
  # Redshift
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the redshift (z) from observed and emitted wavelengths.

  ## Parameters
    - λ_obs:   Observed wavelength (any consistent unit)
    - λ_emit:  Emitted (rest-frame) wavelength (same unit)

  ## Returns
    Redshift value z

  ## Examples
      iex> Astrometry.redshift(700, 656)
      0.06707317073170732
  """
  @spec redshift(number, number) :: float
  def redshift(λ_obs, λ_emit), do: (λ_obs - λ_emit) / λ_emit

  @doc """
  Calculates redshift using the ratio method (equivalent to `redshift/2`).

  ## Examples
      iex> Astrometry.redshift_ratio(700, 656)
      0.06707317073170732
  """
  @spec redshift_ratio(number, number) :: float
  def redshift_ratio(λ_obs, λ_emit), do: λ_obs / λ_emit - 1

  @doc """
  Calculates the recession velocity from redshift (non-relativistic, z ≪ 1).

  ## Parameters
    - z: Redshift
    - c: Speed of light in km/s (default: 299_792.458)

  ## Returns
    Recession velocity in km/s

  ## Examples
      iex> Astrometry.recession_velocity(0.1)
      29979.2458
  """
  @spec recession_velocity(number, number) :: float
  def recession_velocity(z, c \\ 299_792.458), do: z * c

  @doc """
  Calculates recession velocity from redshift using the relativistic Doppler formula.

  v/c = ((1+z)² − 1) / ((1+z)² + 1)

  ## Parameters
    - z: Redshift
    - c: Speed of light in km/s (default: 299_792.458)

  ## Returns
    Velocity in km/s

  ## Examples
      iex> Astrometry.relativistic_velocity(0.5) |> Float.round(0)
      115305.0
  """
  @spec relativistic_velocity(number, number) :: float
  def relativistic_velocity(z, c \\ 299_792.458) do
    z2 = (1 + z) * (1 + z)
    c * (z2 - 1) / (z2 + 1)
  end

  # ---------------------------------------------------------------------------
  # Magnitude & Flux
  # ---------------------------------------------------------------------------

  @doc """
  Calculates apparent magnitude difference from flux ratio.

  ## Parameters
    - f:  Flux of the object
    - f0: Reference flux

  ## Returns
    Magnitude difference (m - m0)

  ## Examples
      iex> Astrometry.apparent_magnitude_diff(1.0, 6.31) |> Float.round(4)
      -2.0
  """
  @spec apparent_magnitude_diff(number, number) :: float
  def apparent_magnitude_diff(f, f0), do: -2.5 * :math.log10(f / f0)

  @doc """
  Calculates absolute magnitude from apparent magnitude and distance.

  ## Parameters
    - m: Apparent magnitude
    - d: Distance in parsecs

  ## Returns
    Absolute magnitude M

  ## Examples
      iex> Astrometry.absolute_magnitude(5, 100)
      0.0
  """
  @spec absolute_magnitude(number, number) :: float
  def absolute_magnitude(m, d), do: m - 5 * :math.log10(d / 10)

  @doc """
  Calculates absolute magnitude with interstellar extinction correction.

  ## Parameters
    - m: Apparent magnitude
    - d: Distance in parsecs
    - a: Extinction in magnitudes

  ## Returns
    Corrected absolute magnitude M

  ## Examples
      iex> Astrometry.absolute_magnitude_extinction(5, 100, 0.3)
      -0.3
  """
  @spec absolute_magnitude_extinction(number, number, number) :: float
  def absolute_magnitude_extinction(m, d, a), do: m - 5 * :math.log10(d / 10) - a

  @doc """
  Calculates the distance modulus: μ = 5 log₁₀(d / 10 pc)

  ## Returns
    Distance modulus μ

  ## Examples
      iex> Astrometry.distance_modulus(100)
      0.0
      iex> Astrometry.distance_modulus(1000) |> Float.round(4)
      5.0
  """
  @spec distance_modulus(number) :: float
  def distance_modulus(d), do: 5 * :math.log10(d / 10)

  @doc """
  Calculates distance in parsecs from the distance modulus: d = 10^(1 + μ/5)

  ## Examples
      iex> Astrometry.distance_from_modulus(0.0)
      10.0
  """
  @spec distance_from_modulus(number) :: float
  def distance_from_modulus(mu), do: 10 * :math.pow(10, mu / 5)

  @doc """
  Calculates the flux ratio from magnitude difference.

  ## Examples
      iex> Astrometry.flux_ratio_from_magnitudes(10, 15)
      100.0
  """
  @spec flux_ratio_from_magnitudes(number, number) :: float
  def flux_ratio_from_magnitudes(ma, mb), do: :math.pow(10, 0.4 * (mb - ma))

  @doc """
  Calculates flux from magnitude and a reference (zero-point) flux: F = F₀ × 10^(-0.4m)

  ## Examples
      iex> Astrometry.flux_from_magnitude(0, 1.0)
      1.0
  """
  @spec flux_from_magnitude(number, number) :: float
  def flux_from_magnitude(m, f0), do: f0 * :math.pow(10, -0.4 * m)

  @doc """
  Converts apparent magnitude to flux density in milli-Janskys (AB system).

  F_ν [mJy] = F₀ × 10^(-0.4 m) × 1000  where F₀ = 3631 Jy for AB magnitudes.

  ## Parameters
    - m:       AB magnitude
    - zero_pt: Zero-point flux density in Jy (default: 3631 Jy for AB)

  ## Returns
    Flux density in milli-Janskys (mJy)

  ## Examples
      iex> Astrometry.magnitude_to_flux_density(8.9) |> Float.round(0)
      3631.0
  """
  @spec magnitude_to_flux_density(number, number) :: float
  def magnitude_to_flux_density(m, zero_pt \\ 3631.0) do
    zero_pt * :math.pow(10, -0.4 * m) * 1000.0
  end

  @doc """
  Calculates color index (magnitude difference between two filters).

  C = −2.5 log₁₀(F₁/F₂)

  ## Examples
      iex> Astrometry.color_index(1.0, 2.0) |> Float.round(4)
      0.7526
  """
  @spec color_index(number, number) :: float
  def color_index(f_f1, f_f2), do: -2.5 * :math.log10(f_f1 / f_f2)

  @doc """
  Color excess from reddening: E(B-V) = (B-V)_obs − (B-V)_intrinsic

  Quantifies the additional reddening due to interstellar dust along the line of sight.

  ## Parameters
    - bv_observed:  Observed B-V color index
    - bv_intrinsic: Intrinsic (unreddened) B-V color index

  ## Returns
    Color excess E(B-V) in magnitudes

  ## Examples
      iex> Astrometry.color_excess(0.5, 0.3) |> Float.round(4)
      0.2
  """
  @spec color_excess(number, number) :: float
  def color_excess(bv_observed, bv_intrinsic), do: bv_observed - bv_intrinsic

  @doc """
  Visual extinction from color excess: A_V = R_V × E(B-V)

  R_V ≈ 3.1 for the standard diffuse ISM. Higher values (up to ~5.5) occur
  in dense molecular clouds.

  ## Parameters
    - ebv: Color excess E(B-V)
    - r_v: Total-to-selective extinction ratio (default: 3.1)

  ## Returns
    Visual extinction A_V in magnitudes

  ## Examples
      iex> Astrometry.extinction_from_color_excess(0.2) |> Float.round(4)
      0.62
  """
  @spec extinction_from_color_excess(number, number) :: float
  def extinction_from_color_excess(ebv, r_v \\ 3.1), do: r_v * ebv

  @doc """
  Calculates bolometric magnitude from absolute magnitude and bolometric correction.

  M_bol = M_V + BC

  ## Examples
      iex> Astrometry.bolometric_magnitude(4.83, -0.07)
      4.76
  """
  @spec bolometric_magnitude(number, number) :: float
  def bolometric_magnitude(m_v, bc), do: m_v + bc

  @doc """
  Calculates luminosity from bolometric magnitude (relative to Sun).

  L/L☉ = 10^(0.4 × (M_bol,☉ − M_bol))

  ## Parameters
    - m_bol:      Bolometric magnitude
    - m_bol_sun:  Solar bolometric magnitude (default: 4.74)

  ## Returns
    Luminosity in solar units L/L☉

  ## Examples
      iex> Astrometry.luminosity_from_bolometric(4.74) |> Float.round(4)
      1.0
  """
  @spec luminosity_from_bolometric(number, number) :: float
  def luminosity_from_bolometric(m_bol, m_bol_sun \\ 4.74) do
    :math.pow(10, 0.4 * (m_bol_sun - m_bol))
  end

  @doc """
  Surface brightness in magnitudes per arcsec².

  μ = −2.5 log₁₀(F_total / Ω)  where Ω is the aperture area in arcsec².

  ## Parameters
    - total_flux:      Total flux within the aperture
    - area_arcsec_sq:  Aperture area in arcsec²

  ## Returns
    Surface brightness in magnitudes per arcsec² (larger = fainter)

  ## Examples
      iex> Astrometry.surface_brightness(1000.0, 100.0) < Astrometry.surface_brightness(100.0, 100.0)
      true
  """
  @spec surface_brightness(number, number) :: float
  def surface_brightness(total_flux, area_arcsec_sq) do
    -2.5 * :math.log10(total_flux / area_arcsec_sq)
  end

  # ---------------------------------------------------------------------------
  # Parallax & Proper Motion
  # ---------------------------------------------------------------------------

  @doc """
  Calculates distance in parsecs from parallax angle: d = 1/p

  ## Examples
      iex> Astrometry.distance_from_parallax(0.1)
      10.0
  """
  @spec distance_from_parallax(number) :: float
  def distance_from_parallax(p), do: 1.0 / p

  @doc """
  Calculates parallax angle in arcseconds from distance.

  ## Examples
      iex> Astrometry.parallax_from_distance(10)
      0.1
  """
  @spec parallax_from_distance(number) :: float
  def parallax_from_distance(d), do: 1.0 / d

  @doc """
  Calculates total proper motion from RA and Dec components.

  μ_total = √(μ_α*² + μ_δ²)

  ## Examples
      iex> Astrometry.proper_motion_total(3.0, 4.0)
      5.0
  """
  @spec proper_motion_total(number, number) :: float
  def proper_motion_total(mu_ra, mu_dec), do: :math.sqrt(mu_ra ** 2 + mu_dec ** 2)

  @doc """
  Position angle of proper motion vector (east of north).

  ## Parameters
    - mu_ra:  Proper motion component in RA (arcsec/yr, including cos δ)
    - mu_dec: Proper motion component in Dec (arcsec/yr)

  ## Returns
    Position angle in radians [0, 2π)

  ## Examples
      iex> Astrometry.proper_motion_angle(0.0, 1.0) |> Float.round(4)
      0.0
  """
  @spec proper_motion_angle(number, number) :: float
  def proper_motion_angle(mu_ra, mu_dec) do
    pa = :math.atan2(mu_ra, mu_dec)
    :math.fmod(pa + 2 * :math.pi(), 2 * :math.pi())
  end

  @doc """
  Calculates the transverse (tangential) velocity from proper motion and distance.

  v_T = 4.74047 × μ × d  (km/s, with μ in arcsec/yr, d in pc)

  ## Examples
      iex> Astrometry.transverse_velocity(0.1, 10) |> Float.round(3)
      4.74
  """
  @spec transverse_velocity(number, number) :: float
  def transverse_velocity(mu, d), do: 4.74 * mu * d

  @doc """
  Calculates total space velocity from transverse and radial components.

  ## Examples
      iex> Astrometry.space_velocity(3.0, 4.0)
      5.0
  """
  @spec space_velocity(number, number) :: float
  def space_velocity(v_t, v_r), do: :math.sqrt(v_t ** 2 + v_r ** 2)

  # ---------------------------------------------------------------------------
  # Angular Sizes & Cosmological Distances
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the angular diameter (in arcseconds) of an object.

  ## Examples
      iex> Astrometry.angular_diameter(1.0, 206_265.0) |> Float.round(6)
      1.0
  """
  @spec angular_diameter(number, number) :: float
  def angular_diameter(physical_size, distance), do: @arcsec_per_rad * physical_size / distance

  @doc """
  Calculates physical size from angular diameter and distance.

  ## Examples
      iex> Astrometry.physical_size(1.0, 206_265.0) |> Float.round(4)
      1.0
  """
  @spec physical_size(number, number) :: float
  def physical_size(angle_arcsec, distance), do: angle_arcsec * distance / @arcsec_per_rad

  @doc """
  Comoving distance (low-z approximation): D_c ≈ cz/H₀

  ## Returns
    Comoving distance in Mpc

  ## Examples
      iex> Astrometry.hubble_distance(0.1) |> Float.round(2)
      428.28
  """
  @spec hubble_distance(number, number) :: float
  def hubble_distance(z, h0 \\ 70.0), do: z * 299_792.458 / h0

  @doc """
  Luminosity distance: D_L = (1 + z) × D_c

  Used to convert observed flux to luminosity: L = 4π D_L² F.

  ## Returns
    Luminosity distance in Mpc

  ## Examples
      iex> Astrometry.luminosity_distance(0.1) |> Float.round(2)
      471.11
  """
  @spec luminosity_distance(number, number) :: float
  def luminosity_distance(z, h0 \\ 70.0) do
    hubble_distance(z, h0) * (1 + z)
  end

  @doc """
  Angular diameter distance: D_A = D_c / (1 + z)

  Relates the observed angular size θ to the physical transverse size: θ = d / D_A.

  ## Returns
    Angular diameter distance in Mpc

  ## Examples
      iex> Astrometry.angular_diameter_distance(0.1) |> Float.round(2)
      389.35
  """
  @spec angular_diameter_distance(number, number) :: float
  def angular_diameter_distance(z, h0 \\ 70.0) do
    hubble_distance(z, h0) / (1 + z)
  end

  @doc """
  Distance modulus from luminosity distance in Mpc.

  μ = 5 log₁₀(D_L / 10 pc)

  ## Examples
      iex> Astrometry.distance_modulus_from_dl(0.01) |> Float.round(4)
      0.0
  """
  @spec distance_modulus_from_dl(number) :: float
  def distance_modulus_from_dl(d_l_mpc) do
    d_l_pc = d_l_mpc * 1.0e6
    5 * :math.log10(d_l_pc / 10)
  end

  # ---------------------------------------------------------------------------
  # Metallicity
  # ---------------------------------------------------------------------------

  @doc """
  Calculates metallicity [Fe/H] relative to solar values.

  [Fe/H] = log₁₀((Fe/H)_star / (Fe/H)_☉)

  ## Examples
      iex> Astrometry.metallicity(0.001, 0.001)
      0.0
  """
  @spec metallicity(number, number) :: float
  def metallicity(fe_h_star, fe_h_solar), do: :math.log10(fe_h_star / fe_h_solar)

  @doc """
  Converts [Fe/H] to metallicity fraction Z (approximate solar-scaled).

  Z ≈ 0.017 × 10^[Fe/H]

  ## Examples
      iex> Astrometry.feh_to_z(0.0) |> Float.round(4)
      0.017
  """
  @spec feh_to_z(number) :: float
  def feh_to_z(fe_h), do: 0.017 * :math.pow(10, fe_h)
end
