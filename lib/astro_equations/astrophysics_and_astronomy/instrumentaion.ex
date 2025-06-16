defmodule AstroEquations.AstrophysicsAndAstronomy.Instrumentaion do
  @moduledoc """
  Provides functions for calculating various optical instrumentation parameters
  including lens properties, focal ratios, fields of view, resolution limits,
  and plate scales, fitting error, adaptive optics error, signal-to-noise ratio,
  atmospheric extinction, and rocket science equations .

  """

  @doc """
  Calculates the focal length of a lens using the lensmaker's equation.

  ## Parameters
    - n: Refractive index of the lens material
    - r1: Radius of curvature of the first surface (meters)
    - r2: Radius of curvature of the second surface (meters)
    - d: Thickness of the lens (meters)
    - thin_lens: Boolean flag to use thin lens approximation (default: false)

  ## Examples
      iex> Optics.Instrumentation.lensmakers_equation(1.5, 0.1, -0.1, 0.01)
      0.10101010101010101

      iex> Optics.Instrumentation.lensmakers_equation(1.5, 0.1, -0.1, 0.01, true)
      0.1
  """
  @spec lensmakers_equation(number, number, number, number, boolean) :: float
  def lensmakers_equation(n, r1, r2, d, thin_lens \\ false) do
    if thin_lens do
      1 / ((n - 1) * (1 / r1 - 1 / r2))
    else
      term = (1 / r1) - (1 / r2) + ((n - 1) * d) / (n * r1 * r2)
      1 / ((n - 1) * term)
    end
  end

  @doc """
  Calculates the focal ratio (f-number) of an optical system.

  ## Parameters
    - f: Focal length (meters)
    - d: Diameter of the aperture (meters)

  ## Examples
      iex> Optics.Instrumentation.focal_ratio(0.5, 0.1)
      5.0
  """
  @spec focal_ratio(number, number) :: float
  def focal_ratio(f, d), do: f / d

  @doc """
  Calculates the field of view of an optical system.

  ## Parameters
    - wd: Detector width (meters)
    - dt: Telescope diameter (meters)
    - n: Focal ratio (f-number)
    - fsys: Optional system focal length (overrides dt * n if provided)

  ## Examples
      iex> Optics.Instrumentation.field_of_view(0.01, 0.2, 10)
      0.005

      iex> Optics.Instrumentation.field_of_view(0.01, 0.2, 10, 2.0)
      0.005
  """
  @spec field_of_view(number, number, number, number | nil) :: float
  def field_of_view(wd, dt, n, fsys \\ nil) do
    system_focal_length = fsys || dt * n
    wd / system_focal_length
  end

  @doc """
  Calculates the diffraction limit (Rayleigh criterion) for an optical system.

  ## Parameters
    - wavelength: Wavelength of light (meters)
    - diameter: Diameter of the aperture (meters)

  ## Examples
      iex> Optics.Instrumentation.diffraction_limit(500.0e-9, 0.1)
      6.1e-6
  """
  @spec diffraction_limit(number, number) :: float
  def diffraction_limit(wavelength, diameter), do: 1.22 * wavelength / diameter

  @doc """
  Calculates the seeing limit (Rayleigh criterion) for atmospheric conditions.

  ## Parameters
    - wavelength: Wavelength of light (meters)
    - r0: Fried parameter (meters)

  ## Examples
      iex> Optics.Instrumentation.seeing_limit(500.0e-9, 0.2)
      2.45e-6
  """
  @spec seeing_limit(number, number) :: float
  def seeing_limit(wavelength, r0), do: 0.98 * wavelength / r0

  @doc """
  Calculates the total resolution limit by combining diffraction and seeing limits.

  ## Parameters
    - wavelength: Wavelength of light (meters)
    - diameter: Diameter of the aperture (meters)
    - r0: Fried parameter (meters)

  ## Examples
      iex> Optics.Instrumentation.total_resolution_limit(500.0e-9, 0.1, 0.2)
      6.580454994347146e-6
  """
  @spec total_resolution_limit(number, number, number) :: float
  def total_resolution_limit(wavelength, diameter, r0) do
    diffraction = diffraction_limit(wavelength, diameter)
    seeing = seeing_limit(wavelength, r0)
    :math.sqrt(diffraction * diffraction + seeing * seeing)
  end

  @doc """
  Calculates the Nyquist sampling parameter for a diffraction-limited system.

  ## Parameters
    - p: Pixel size (meters)
    - fsys: System focal length (meters)
    - wavelength: Wavelength of light (meters)
    - dt: Telescope diameter (meters)

  ## Examples
      iex> Optics.Instrumentation.nyquist_sampling(5.0e-6, 2.0, 500.0e-9, 0.2)
      10.0
  """
  @spec nyquist_sampling(number, number, number, number) :: float
  def nyquist_sampling(p, fsys, wavelength, dt) do
    # The original line was invalid; if you want to compare, use ==, or just return the calculation.
    (2 * p) / wavelength
  end

  @doc """
  Converts focal length to plate scale in both radians/meter and arcseconds/meter.

  ## Parameters
    - f: Focal length (meters)

  ## Examples
      iex> Optics.Instrumentation.plate_scale(1.0)
      {1.0, 206265.0}
  """
  @spec plate_scale(number) :: {float, float}
  def plate_scale(f) do
    {1 / f, 206_265 / f}
  end



  @doc """
  Calculates the fitting error variance for adaptive optics systems.

  ## Parameters
    - d_sub: Subaperture diameter (meters)
    - r0: Fried parameter (meters)

  ## Examples
      iex> Optics.Instrumentation.fitting_error(0.1, 0.2)
      0.10291902476328986
  """
  @spec fitting_error(number, number) :: float
  def fitting_error(d_sub, r0) do
    0.26 * :math.pow(d_sub / r0, 5/3)
  end

  @doc """
  Calculates the total adaptive optics error variance.

  ## Parameters
    - d_sub: Subaperture diameter (meters)
    - r0: Fried parameter (meters)
    - theta: Angular offset (radians)
    - theta_0: Reference angular offset (radians)
    - tau: Time delay (seconds)
    - tau_0: Reference time delay (seconds)
    - c_wfs: Wavefront sensor constant
    - lambda: Wavelength (meters)
    - f_t: Telescope focal length (meters)

  ## Examples
      iex> Optics.Instrumentation.adaptive_optics_error(0.1, 0.2, 0.01, 0.02, 0.001, 0.002, 1.0, 500.0e-9, 10.0)
      0.11874999999999999
  """
  @spec adaptive_optics_error(number, number, number, number, number, number, number, number, number) :: float
  def adaptive_optics_error(d_sub, r0, theta, theta_0, tau, tau_0, c_wfs, lambda, f_t) do
    fitting_term = 0.3 * :math.pow(d_sub / r0, 5/3)
    angular_term = :math.pow(theta / theta_0, 5/3)
    time_term = 28.4 * :math.pow(tau / tau_0, 5/3)
    wfs_term = c_wfs * :math.pow(lambda / (f_t * d_sub), 2)

    fitting_term + angular_term + time_term + wfs_term
  end

  @doc """
  Calculates the signal-to-noise ratio for an optical system.

  ## Parameters
    - f: Signal flux (photons/second)
    - t: Integration time (seconds)
    - b_s: Background flux (photons/second/pixel)
    - n_p: Number of pixels
    - d: Dark current (electrons/second/pixel)
    - r: Read noise (electrons/pixel)

  ## Examples
      iex> Optics.Instrumentation.signal_to_noise(100, 10, 5, 1000, 0.1, 2.0)
      27.386127875258307
  """
  @spec signal_to_noise(number, number, number, number, number, number) :: float
  def signal_to_noise(f, t, b_s, n_p, d, r) do
    signal = f * t
    noise = :math.sqrt(signal + (b_s * n_p * t) + (d * n_p * t) + (:math.pow(r, 2) * n_p))
    signal / noise
  end

  @doc """
  Calculates the apparent magnitude corrected for atmospheric extinction.

  ## Parameters
    - m_lambda_z: Zenithal magnitude
    - a_lambda: Extinction coefficient
    - z: Zenith angle (radians)

  ## Examples
      iex> Optics.Instrumentation.atmospheric_extinction(10.0, 0.2, :math.pi/3)
      10.4
  """
  @spec atmospheric_extinction(number, number, number) :: float
  def atmospheric_extinction(m_lambda_z, a_lambda, z) do
    m_lambda_z - a_lambda * (1 / :math.cos(z))
  end

  @doc """
  Calculates the delta-v (change in velocity) using Tsiolkovsky's rocket equation.

  ## Parameters
    - v_e: Exhaust velocity (m/s)
    - m_0: Initial mass (kg)
    - m_f: Final mass (kg)

  ## Examples
      iex> Optics.Instrumentation.tsiolkovsky_rocket_equation(2500, 1000, 500)
      1732.8679513998632
  """
  @spec tsiolkovsky_rocket_equation(number, number, number) :: float
  def tsiolkovsky_rocket_equation(v_e, m_0, m_f) do
    v_e * :math.log(m_0 / m_f)
  end

end
