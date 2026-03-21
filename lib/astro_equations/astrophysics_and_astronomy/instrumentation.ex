defmodule AstroEquations.AstrophysicsAndAstronomy.Instrumentation do
  @moduledoc """
  Optical instrumentation parameters for astronomical telescopes and detectors.

  Covers:
  - Lensmaker's equation and focal ratio
  - Fields of view, plate scales, and image scale
  - Diffraction and seeing (atmospheric) resolution limits
  - Nyquist sampling criterion
  - Adaptive optics (fitting error, total AO error, Strehl ratio)
  - Signal-to-noise ratio and limiting magnitude
  - CCD characterisation (dynamic range, saturation time)
  - Spectral resolving power and grating dispersion
  - Atmospheric extinction and airmass
  - Tsiolkovsky rocket equation and specific impulse
  - Telescope collecting area and étendue
  - Photon count rate from source magnitude
  - Optimal extraction aperture radius
  """

  @arcsec_per_rad 206_265.0

  # ---------------------------------------------------------------------------
  # Optics
  # ---------------------------------------------------------------------------

  @doc """
  Calculates the focal length using the lensmaker's equation.

  Supports both thin-lens approximation and full thick-lens formula.

  Thin:  1/f = (n−1)(1/R₁ − 1/R₂)
  Thick: 1/f = (n−1)[1/R₁ − 1/R₂ + (n−1)d/(nR₁R₂)]

  ## Parameters
    - n:         Refractive index
    - r1:        Radius of curvature of first surface (m)
    - r2:        Radius of curvature of second surface (m)
    - d:         Lens thickness (m)
    - thin_lens: Use thin-lens approximation (default: false)

  ## Returns
    Focal length in meters

  ## Examples
      iex> Instrumentation.lensmakers_equation(1.5, 0.1, -0.1, 0.01)
      0.10101010101010101
      iex> Instrumentation.lensmakers_equation(1.5, 0.1, -0.1, 0.01, true)
      0.1
  """
  @spec lensmakers_equation(number, number, number, number, boolean) :: float
  def lensmakers_equation(n, r1, r2, d, thin_lens \\ false) do
    if thin_lens do
      1 / ((n - 1) * (1 / r1 - 1 / r2))
    else
      term = 1 / r1 - 1 / r2 + (n - 1) * d / (n * r1 * r2)
      1 / ((n - 1) * term)
    end
  end

  @doc """
  Calculates the focal ratio (f-number): N = f / D.

  ## Examples
      iex> Instrumentation.focal_ratio(0.5, 0.1)
      5.0
  """
  @spec focal_ratio(number, number) :: float
  def focal_ratio(f, d), do: f / d

  @doc """
  Calculates the field of view in radians: FoV = detector_width / focal_length

  ## Parameters
    - wd:   Detector width (m)
    - dt:   Telescope diameter (m)
    - n:    Focal ratio
    - fsys: System focal length (m); if nil, uses dt × n

  ## Returns
    Field of view in radians

  ## Examples
      iex> Instrumentation.field_of_view(0.01, 0.2, 10)
      0.005
  """
  @spec field_of_view(number, number, number, number | nil) :: float
  def field_of_view(wd, dt, n, fsys \\ nil) do
    fl = fsys || dt * n
    wd / fl
  end

  @doc """
  Converts field of view from radians to arcseconds.

  ## Examples
      iex> Instrumentation.fov_to_arcsec(0.001) |> Float.round(3)
      206.265
  """
  @spec fov_to_arcsec(number) :: float
  def fov_to_arcsec(fov_rad), do: fov_rad * @arcsec_per_rad

  @doc """
  Telescope collecting area: A = π (D/2)²

  ## Examples
      iex> Instrumentation.collecting_area(0.5) |> Float.round(6)
      0.196350
  """
  @spec collecting_area(number) :: float
  def collecting_area(diameter) do
    :math.pi() * :math.pow(diameter / 2, 2)
  end

  @doc """
  Étendue (throughput) of an optical system: A Ω

  A fundamental invariant of any optical system; limits the maximum
  flux that can be collected and focused.

  ## Parameters
    - area:       Collecting area (m²)
    - solid_angle: Field solid angle (steradians)

  ## Returns
    Étendue in m² sr

  ## Examples
      iex> Instrumentation.etendue(0.2, 1.0e-6) |> Float.round(10)
      2.0e-7
  """
  @spec etendue(number, number) :: float
  def etendue(area, solid_angle), do: area * solid_angle

  # ---------------------------------------------------------------------------
  # Resolution
  # ---------------------------------------------------------------------------

  @doc """
  Rayleigh diffraction limit in radians: θ = 1.22 λ / D

  ## Examples
      iex> Instrumentation.diffraction_limit(500.0e-9, 0.1)
      6.1e-6
  """
  @spec diffraction_limit(number, number) :: float
  def diffraction_limit(wavelength, diameter), do: 1.22 * wavelength / diameter

  @doc """
  Seeing (atmospheric) resolution limit in radians: θ_see = 0.98 λ / r₀

  ## Examples
      iex> Instrumentation.seeing_limit(500.0e-9, 0.15) > 0
      true
  """
  @spec seeing_limit(number, number) :: float
  def seeing_limit(wavelength, r0), do: 0.98 * wavelength / r0

  @doc """
  Combines diffraction and seeing limits in quadrature: θ_total = √(θ_diff² + θ_see²)

  ## Examples
      iex> Instrumentation.total_resolution_limit(500.0e-9, 0.4, 0.2) > 0
      true
  """
  @spec total_resolution_limit(number, number, number) :: float
  def total_resolution_limit(wavelength, diameter, r0) do
    d = diffraction_limit(wavelength, diameter)
    s = seeing_limit(wavelength, r0)
    :math.sqrt(d * d + s * s)
  end

  @doc """
  Diffraction-limited angular resolution in arcseconds.

  θ_arcsec = 0.2516 × λ_μm / D_m  (Dawes limit approximation)

  ## Parameters
    - wavelength_um: Wavelength in microns
    - diameter_m:    Aperture diameter in meters

  ## Returns
    Resolution in arcseconds

  ## Examples
      iex> Instrumentation.resolution_arcsec(0.5, 0.1) |> Float.round(4)
      1.258
  """
  @spec resolution_arcsec(number, number) :: float
  def resolution_arcsec(wavelength_um, diameter_m) do
    0.2516 * wavelength_um / diameter_m
  end

  # ---------------------------------------------------------------------------
  # Sampling & Plate Scale
  # ---------------------------------------------------------------------------

  @doc """
  Nyquist sampling: number of pixels per resolution element.

  Returns the number of pixels that span one resolution element (should be ≥ 2
  to satisfy the Nyquist criterion). Values < 2 indicate under-sampling.

  n_pix = 1.22 λ f_sys / (D × pixel_size)

  ## Parameters
    - p:          Pixel size (m)
    - fsys:       System focal length (m)
    - wavelength: Wavelength (m)
    - dt:         Telescope diameter (m)

  ## Returns
    Pixels per resolution element

  ## Examples
      iex> Instrumentation.nyquist_sampling(6.0e-6, 2.0, 500.0e-9, 0.1) > 0
      true
  """
  @spec nyquist_sampling(number, number, number, number) :: float
  def nyquist_sampling(p, fsys, wavelength, dt) do
    resolution_element = 1.22 * wavelength * fsys / dt
    resolution_element / p
  end

  @doc """
  Plate scale in radians/m and arcsec/m: {1/f, 206265/f}

  ## Examples
      iex> Instrumentation.plate_scale(1.0)
      {1.0, 206265.0}
  """
  @spec plate_scale(number) :: {float, float}
  def plate_scale(f), do: {1 / f, @arcsec_per_rad / f}

  @doc """
  Image scale in arcsec/pixel: scale = (arcsec_per_rad × pixel_size) / focal_length

  ## Examples
      iex> Instrumentation.image_scale(6.0e-6, 2.0) > 0
      true
  """
  @spec image_scale(number, number) :: float
  def image_scale(pixel_size, f), do: @arcsec_per_rad * pixel_size / f

  @doc """
  Photon count rate from an apparent magnitude source.

  N_phot = F₀ × 10^(-0.4 m) × A × QE × Δλ/λ

  where F₀ is the zero-point photon flux in photons/s/m²/Δλ/λ.

  ## Parameters
    - magnitude:  Apparent AB magnitude
    - area_m2:    Collecting area (m²)
    - qe:         Detector quantum efficiency (0–1, default: 1.0)
    - f0_photons: Zero-point photon flux (photons/s/m², default: 1.0e7 for visible)

  ## Returns
    Photon count rate (photons/s)

  ## Examples
      iex> Instrumentation.photon_count_rate(0.0, 1.0) > 0
      true
  """
  @spec photon_count_rate(number, number, number, number) :: float
  def photon_count_rate(magnitude, area_m2, qe \\ 1.0, f0_photons \\ 1.0e7) do
    f0_photons * :math.pow(10, -0.4 * magnitude) * area_m2 * qe
  end

  # ---------------------------------------------------------------------------
  # Adaptive Optics
  # ---------------------------------------------------------------------------

  @doc """
  Wavefront fitting error variance: σ_fit² = 0.26 × (d_sub/r₀)^(5/3)

  ## Examples
      iex> Instrumentation.fitting_error(0.1, 0.2) |> Float.round(4)
      0.103
  """
  @spec fitting_error(number, number) :: float
  def fitting_error(d_sub, r0), do: 0.26 * :math.pow(d_sub / r0, 5 / 3)

  @doc """
  Total AO error variance from fitting, anisoplanatism, temporal, and WFS noise terms.

  σ² = σ_fit² + σ_aniso² + σ_temp² + σ_wfs²

  ## Parameters
    - d_sub:   Subaperture diameter (m)
    - r0:      Fried parameter (m)
    - theta:   Angular offset from guide star (rad)
    - theta_0: Isoplanatic angle (rad)
    - tau:     Time delay (s)
    - tau_0:   Coherence time (s)
    - c_wfs:   Wavefront sensor noise coefficient
    - lambda:  Wavelength (m)
    - f_t:     Telescope focal length (m)

  ## Returns
    Total AO error variance in rad²

  ## Examples
      iex> Instrumentation.adaptive_optics_error(0.1, 0.15, 0.01, 0.2, 0.5, 0.02, 10, 0.3, 0.001) > 0
      true
  """
  @spec adaptive_optics_error(
          number,
          number,
          number,
          number,
          number,
          number,
          number,
          number,
          number
        ) ::
          float
  def adaptive_optics_error(d_sub, r0, theta, theta_0, tau, tau_0, c_wfs, lambda, f_t) do
    fitting_term = 0.3 * :math.pow(d_sub / r0, 5 / 3)
    angular_term = :math.pow(theta / theta_0, 5 / 3)
    time_term = 28.4 * :math.pow(tau / tau_0, 5 / 3)
    wfs_term = c_wfs * :math.pow(lambda / (f_t * d_sub), 2)

    fitting_term + angular_term + time_term + wfs_term
  end

  @doc """
  Strehl ratio (Maréchal approximation): S ≈ exp(−σ²)

  Valid for σ² < 1 rad²; S = 1 for a diffraction-limited system.

  ## Examples
      iex> Instrumentation.strehl_ratio(0.0) |> Float.round(4)
      1.0
  """
  @spec strehl_ratio(number) :: float
  def strehl_ratio(sigma_sq), do: :math.exp(-sigma_sq)

  # ---------------------------------------------------------------------------
  # Detector / SNR
  # ---------------------------------------------------------------------------

  @doc """
  Signal-to-noise ratio for a CCD detector.

  SNR = S·t / √(S·t + B_s·N_p·t + D·N_p·t + R²·N_p)

  ## Parameters
    - f:   Source flux (photons/s)
    - t:   Integration time (s)
    - b_s: Sky background flux (photons/s/pixel)
    - n_p: Number of pixels
    - d:   Dark current (electrons/s/pixel)
    - r:   Read noise (electrons/pixel/read)

  ## Returns
    SNR

  ## Examples
      iex> Instrumentation.signal_to_noise(100, 10, 5, 1000, 0.1, 2.0) |> Float.round(3)
      27.386
  """
  @spec signal_to_noise(number, number, number, number, number, number) :: float
  def signal_to_noise(f, t, b_s, n_p, d, r) do
    signal = f * t
    noise = :math.sqrt(signal + b_s * n_p * t + d * n_p * t + :math.pow(r, 2) * n_p)
    signal / noise
  end

  @doc """
  Required exposure time to reach a target SNR (sky-background-limited case).

  t = (SNR)² × (B_s × N_p + D × N_p) / (f - SNR² × ???)
  Simplified sky-limit: t ≈ SNR² × N_noise / f²

  Returns the quadratic-formula solution for t in the full noise model.

  ## Parameters
    - target_snr: Desired SNR
    - f:          Source flux (photons/s)
    - b_s:        Sky background (photons/s/pixel)
    - n_p:        Number of pixels
    - d:          Dark current (electrons/s/pixel)
    - r:          Read noise (electrons/pixel)

  ## Returns
    Required exposure time in seconds

  ## Examples
      iex> Instrumentation.exposure_time_from_snr(10.0, 100, 5, 4, 0.1, 2.0) > 0
      true
  """
  @spec exposure_time_from_snr(number, number, number, number, number, number) :: float
  def exposure_time_from_snr(target_snr, f, b_s, n_p, d, r) do
    snr2 = target_snr * target_snr
    # Quadratic: f² t² - snr²(f + b_s*n_p + d*n_p) t - snr² r²*n_p = 0
    a = f * f
    b = -snr2 * (f + b_s * n_p + d * n_p)
    c_coeff = -snr2 * :math.pow(r, 2) * n_p
    (-b + :math.sqrt(b * b - 4 * a * c_coeff)) / (2 * a)
  end

  @doc """
  Limiting magnitude for a given SNR threshold (sky-background-limited, simplified).

  m_lim = m₀ − 2.5 log₁₀(SNR × √(B_sky Ω t) / (QE × A × t))

  ## Parameters
    - snr_threshold: Desired minimum SNR
    - f0:            Photon flux for m = 0 (photons/s/m²)
    - area:          Telescope collecting area (m²)
    - t:             Exposure time (s)
    - b_sky:         Sky background (photons/s/arcsec²)
    - omega:         Solid angle per resolution element (arcsec²)
    - qe:            Quantum efficiency (0–1, default: 1.0)

  ## Returns
    Limiting magnitude

  ## Examples
      iex> Instrumentation.limiting_magnitude(10, 5, 100, 0.1, 2.0, 5.0) > 0
      true
  """
  @spec limiting_magnitude(number, number, number, number, number, number, number) :: float
  def limiting_magnitude(snr_threshold, f0, area, t, b_sky, omega, qe \\ 1.0) do
    signal_needed = snr_threshold * :math.sqrt(b_sky * omega * t)
    actual_flux = signal_needed / (qe * area * t)
    -2.5 * :math.log10(actual_flux / f0)
  end

  @doc """
  CCD dynamic range in decibels: DR = 20 log₁₀(FWC / σ_read)

  ## Examples
      iex> Instrumentation.dynamic_range_db(65_000, 5) > 0
      true
  """
  @spec dynamic_range_db(number, number) :: float
  def dynamic_range_db(full_well_capacity, read_noise) do
    20 * :math.log10(full_well_capacity / read_noise)
  end

  @doc """
  CCD saturation time — time until the brightest source fills the full well.

  t_sat = FWC / (source_count_rate + sky_rate × n_pix)

  ## Parameters
    - full_well_capacity: Full-well capacity (electrons)
    - source_rate:        Source count rate (electrons/s)
    - sky_rate:           Sky background rate per pixel (electrons/s/pixel)
    - n_pix:              Number of pixels in extraction aperture

  ## Returns
    Saturation time in seconds

  ## Examples
      iex> Instrumentation.ccd_saturation_time(65_000, 10_000, 50, 9) > 0
      true
  """
  @spec ccd_saturation_time(number, number, number, number) :: float
  def ccd_saturation_time(full_well_capacity, source_rate, sky_rate, n_pix) do
    full_well_capacity / (source_rate + sky_rate * n_pix)
  end

  @doc """
  Optimal circular extraction aperture radius for maximum SNR.

  r_opt ≈ 1.4 × FWHM  (empirical approximation for Gaussian PSF and sky-limited case)

  ## Parameters
    - fwhm_pixels: PSF full-width at half-maximum in pixels

  ## Returns
    Optimal aperture radius in pixels

  ## Examples
      iex> Instrumentation.optimal_extraction_aperture(3.0) |> Float.round(4)
      4.2
  """
  @spec optimal_extraction_aperture(number) :: float
  def optimal_extraction_aperture(fwhm_pixels), do: 1.4 * fwhm_pixels

  # ---------------------------------------------------------------------------
  # Spectroscopy
  # ---------------------------------------------------------------------------

  @doc """
  Spectral resolving power: R = λ / Δλ

  ## Examples
      iex> Instrumentation.resolving_power(550.0e-9, 0.05e-9)
      11000.0
  """
  @spec resolving_power(number, number) :: float
  def resolving_power(lambda, delta_lambda), do: lambda / delta_lambda

  @doc """
  Linear (reciprocal) dispersion of a diffraction grating in the image plane.

  dλ/dx = d cos(β) / (m × f)

  ## Parameters
    - grating_spacing: Grating period d (m)
    - beta:            Diffraction angle (rad)
    - order:           Diffraction order m
    - focal_length:    Camera focal length (m)

  ## Returns
    Reciprocal linear dispersion in m/m (multiply by 10⁶ for nm/mm)

  ## Examples
      iex> Instrumentation.grating_dispersion(1.0e-6, 0.0, 1, 0.5) > 0
      true
  """
  @spec grating_dispersion(number, number, number, number) :: float
  def grating_dispersion(grating_spacing, beta, order, focal_length) do
    grating_spacing * :math.cos(beta) / (order * focal_length)
  end

  @doc """
  Blaze wavelength of a diffraction grating: λ_blaze = 2 d sin(θ_blaze) / m

  ## Examples
      iex> Instrumentation.blaze_wavelength(1000, 0.5236) > 0
      true
  """
  @spec blaze_wavelength(number, number, number) :: float
  def blaze_wavelength(grating_spacing, blaze_angle, order \\ 1) do
    2 * grating_spacing * :math.sin(blaze_angle) / order
  end

  # ---------------------------------------------------------------------------
  # Atmospheric Extinction
  # ---------------------------------------------------------------------------

  @doc """
  Apparent magnitude corrected for atmospheric extinction.

  m_obs = m_zenith + A_λ × sec(z)

  At higher airmass (larger z), stars appear fainter (larger magnitude).

  ## Parameters
    - m_lambda_z: Zenithal magnitude (at z = 0)
    - a_lambda:   Extinction coefficient (mag/airmass)
    - z:          Zenith angle (radians)

  ## Returns
    Observed (extincted) magnitude at zenith angle z

  ## Examples
      iex> Instrumentation.atmospheric_extinction(10.0, 0.2, :math.pi/3)
      10.4
  """
  @spec atmospheric_extinction(number, number, number) :: float
  def atmospheric_extinction(m_lambda_z, a_lambda, z) do
    m_lambda_z + a_lambda * (1 / :math.cos(z))
  end

  @doc """
  Airmass (plane-parallel atmosphere): X = sec(z) = 1 / cos(z)

  ## Examples
      iex> Instrumentation.airmass(0) |> Float.round(4)
      1.0
  """
  @spec airmass(number) :: float
  def airmass(z), do: 1 / :math.cos(z)

  # ---------------------------------------------------------------------------
  # Tsiolkovsky Rocket Equation
  # ---------------------------------------------------------------------------

  @doc """
  Delta-v from Tsiolkovsky's rocket equation: Δv = v_e × ln(m₀/m_f)

  ## Examples
      iex> Instrumentation.tsiolkovsky_rocket_equation(2500, 1000, 500) |> Float.round(3)
      1732.868
  """
  @spec tsiolkovsky_rocket_equation(number, number, number) :: float
  def tsiolkovsky_rocket_equation(v_e, m_0, m_f) do
    v_e * :math.log(m_0 / m_f)
  end

  @doc """
  Specific impulse from exhaust velocity: Isp = v_e / g₀

  ## Examples
      iex> Instrumentation.specific_impulse(4400, 9.80665) |> Float.round(1)
      448.7
  """
  @spec specific_impulse(number, number) :: float
  def specific_impulse(v_e, g0 \\ 9.80665), do: v_e / g0
end
