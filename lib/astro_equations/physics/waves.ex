defmodule AstroEquations.Physics.Waves do
  @moduledoc """
  Wave mechanics and acoustics.

  Covers:
  - Wave fundamentals (number, velocity, angular frequency, period)
  - Wave function (transverse, standing)
  - Superposition and interference (constructive/destructive, path difference)
  - Diffraction (single slit, double slit, grating)
  - Doppler effect (acoustic and relativistic)
  - Sound intensity and decibels
  - Standing waves on a string and in open/closed pipes
  - Group and phase velocity
  - Wave energy and power
  - Snell's Law (refraction), critical angle, Malus's Law
  - Acoustic impedance
  - Thin lens and mirror equations
  """

  # m/s at 20 °C
  @speed_of_sound 343.0
  # m/s
  @speed_of_light 2.99792458e8

  # ---------------------------------------------------------------------------
  # Fundamentals
  # ---------------------------------------------------------------------------

  @doc """
  Wave number: k = 2π / λ

  ## Examples
      iex> Waves.wave_number(2) |> Float.round(6)
      3.141593
  """
  @spec wave_number(number) :: float
  def wave_number(wavelength), do: 2 * :math.pi() / wavelength

  @doc """
  Wave speed from frequency and wavelength: v = f λ

  ## Examples
      iex> Waves.wave_velocity(440, 0.780)
      343.2
  """
  @spec wave_velocity(number, number) :: float
  def wave_velocity(frequency, wavelength), do: frequency * wavelength * 1.0

  @doc "Frequency of a wave from its speed and wavelength: f = v/λ."
  @spec frequency(number, number) :: float
  def frequency(wave_speed, wavelength), do: wave_speed / wavelength

  @doc "Wavelength of a wave from its speed and frequency: λ = v/f."
  @spec wavelength(number, number) :: float
  def wavelength(wave_speed, freq), do: wave_speed / freq

  @doc "Period of oscillation from frequency: T = 1/f."
  @spec period(number) :: float
  def period(frequency), do: 1.0 / frequency

  @doc "Angular frequency from a keyword-list argument: either `period:` or `frequency:`."
  @spec angular_frequency(Keyword.t()) :: float | no_return()
  def angular_frequency(opts) do
    cond do
      Keyword.has_key?(opts, :period) -> 2 * :math.pi() / Keyword.get(opts, :period)
      Keyword.has_key?(opts, :frequency) -> 2 * :math.pi() * Keyword.get(opts, :frequency)
      true -> raise ArgumentError, "Must provide :period or :frequency"
    end
  end

  @doc "Angular frequency from ordinary frequency: omega = 2πf."
  @spec omega_from_frequency(number) :: float
  def omega_from_frequency(f), do: 2 * :math.pi() * f

  # ---------------------------------------------------------------------------
  # Wave Function
  # ---------------------------------------------------------------------------

  @doc """
  Transverse wave displacement: y(x, t) = A sin(omega t - k x + φ)

  ## Parameters
    - amplitude:         A (m)
    - angular_frequency: omega (rad/s)
    - wave_number:       k (rad/m)
    - time:              t (s)
    - position:          x (m)
    - phase_constant:    φ (rad, default: 0)

  ## Examples
      iex> Waves.wave_function(1, 2, 3, 4, 5) |> Float.round(4)
      -0.5366
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

  @doc "Standing wave formed by two counter-propagating waves: y = 2A sin(kx) cos(omegat)."
  @spec standing_wave(number, number, number, number, number) :: float
  def standing_wave(amplitude, k, omega, x, t) do
    2 * amplitude * :math.sin(k * x) * :math.cos(omega * t)
  end

  # ---------------------------------------------------------------------------
  # Superposition & Interference
  # ---------------------------------------------------------------------------

  @doc "Checks the constructive interference condition: Δd = mλ."
  @spec constructive_interference?(number, number, number) :: boolean
  def constructive_interference?(path_diff, wavelength, m \\ 0) do
    abs(path_diff - m * wavelength) < wavelength * 1.0e-9
  end

  @doc "Checks the destructive interference condition: Δd = (m + ½)λ."
  @spec destructive_interference?(number, number, number) :: boolean
  def destructive_interference?(path_diff, wavelength, m \\ 0) do
    abs(path_diff - (m + 0.5) * wavelength) < wavelength * 1.0e-9
  end

  @doc "Beat frequency between two oscillators: f_beat = |f₁ − f₂|."
  @spec beat_frequency(number, number) :: float
  def beat_frequency(f1, f2), do: abs(f1 - f2) * 1.0

  # ---------------------------------------------------------------------------
  # Diffraction
  # ---------------------------------------------------------------------------

  @doc """
  Single-slit diffraction minimum angle: sin θ = m λ / a

  ## Parameters
    - m:          Diffraction order (integer ≠ 0)
    - wavelength: λ (m)
    - slit_width: a (m)

  ## Returns
    Angle θ in radians

  ## Examples
      iex> Waves.single_slit_minimum(1, 500.0e-9, 0.001) > 0
      true
  """
  @spec single_slit_minimum(number, number, number) :: float
  def single_slit_minimum(m, wavelength, slit_width) do
    :math.asin(m * wavelength / slit_width)
  end

  @doc """
  Double-slit bright fringe angle: sin θ = m λ / d

  ## Examples
      iex> Waves.double_slit_bright(0, 500.0e-9, 5.0e-4) |> Float.round(4)
      0.0
  """
  @spec double_slit_bright(number, number, number) :: float
  def double_slit_bright(m, wavelength, slit_sep) do
    :math.asin(m * wavelength / slit_sep)
  end

  @doc """
  Diffraction grating angle: m λ = d sin θ

  ## Returns
    Diffraction angle in radians

  ## Examples
      iex> Waves.grating_angle(0, 500.0e-9, 5.0e-4) |> Float.round(4)
      0.0
  """
  @spec grating_angle(number, number, number) :: float
  def grating_angle(m, wavelength, grating_sep) do
    :math.asin(m * wavelength / grating_sep)
  end

  @doc """
  Resolving power of a diffraction grating: R = m N

  ## Examples
      iex> Waves.grating_resolving_power(1, 500)
      500
  """
  @spec grating_resolving_power(number, number) :: number
  def grating_resolving_power(m, n), do: m * n

  # ---------------------------------------------------------------------------
  # Doppler Effect
  # ---------------------------------------------------------------------------

  @doc """
  Acoustic Doppler shift (observer and source both potentially moving):

  f_obs = f_src * (v_sound + v_obs) / (v_sound - v_src)

  Positive v_obs = observer moving toward source.
  Positive v_src = source moving toward observer.

  ## Parameters
    - f_src:   Source frequency (Hz)
    - v_obs:   Observer velocity (m/s, + toward source)
    - v_src:   Source velocity (m/s, + toward observer)
    - v_sound: Speed of sound (m/s, default: 343)

  ## Examples
      iex> Waves.doppler(440.0, 0.0, 0.0) |> Float.round(4)
      440.0
  """
  @spec doppler(number, number, number, number) :: float
  def doppler(f_src, v_obs, v_src, v_sound \\ @speed_of_sound) do
    f_src * (v_sound + v_obs) / (v_sound - v_src)
  end

  @doc """
  Relativistic Doppler shift: f_obs = f_src √((1 + β) / (1 - β)).

  Positive β = source approaching observer.

  ## Examples
      iex> Waves.relativistic_doppler(1.0e14, 0) |> Float.round(2)
      1.0e14
  """
  @spec relativistic_doppler(number, number, number) :: float
  def relativistic_doppler(f_src, v, c \\ @speed_of_light) do
    b = v / c
    f_src * :math.sqrt((1 + b) / (1 - b))
  end

  # ---------------------------------------------------------------------------
  # Sound Intensity & Decibels
  # ---------------------------------------------------------------------------

  @doc """
  Sound intensity level: L = 10 log₁₀(I / I₀) dB

  ## Parameters
    - intensity: I (W/m²)
    - i0:        Reference intensity (default: 1.0e-12 W/m²)

  ## Examples
      iex> Waves.intensity_decibels(1.0e-12) |> Float.round(4)
      0.0
  """
  @spec intensity_decibels(number, number) :: float
  def intensity_decibels(intensity, i0 \\ 1.0e-12) do
    10 * :math.log10(intensity / i0)
  end

  @doc "Sound intensity from a decibel level: I = I₀ × 10^(L/10)."
  @spec intensity_from_decibels(number, number) :: float
  def intensity_from_decibels(level_db, i0 \\ 1.0e-12) do
    i0 * :math.pow(10, level_db / 10)
  end

  @doc "Intensity ratio from amplitude ratio: I₂/I₁ = (A₂/A₁)²."
  @spec intensity_ratio_from_amplitude(number, number) :: float
  def intensity_ratio_from_amplitude(a2, a1), do: (a2 / a1) ** 2 * 1.0

  @doc """
  Specific acoustic impedance of a medium: Z = ρ v

  ## Parameters
    - density:    Medium density ρ (kg/m³)
    - wave_speed: Speed of sound in medium (m/s)

  ## Returns
    Acoustic impedance Z (Pa·s/m = Rayleigh)

  ## Examples
      iex> Waves.acoustic_impedance(1.225, 343.0) |> Float.round(0)
      420.0
  """
  @spec acoustic_impedance(number, number) :: float
  def acoustic_impedance(density, wave_speed), do: density * wave_speed

  @doc """
  Sound power from intensity and area: P = I × A

  ## Examples
      iex> Waves.sound_power(1.0e-3, 0.1) |> Float.round(6)
      0.0001
  """
  @spec sound_power(number, number) :: float
  def sound_power(intensity, area), do: intensity * area

  # ---------------------------------------------------------------------------
  # Standing Waves
  # ---------------------------------------------------------------------------

  @doc """
  Harmonic frequencies on a string fixed at both ends: f_n = n v / (2 L)

  ## Examples
      iex> Waves.string_harmonic(1, 340, 0.68) |> Float.round(1)
      250.0
  """
  @spec string_harmonic(number, number, number) :: float
  def string_harmonic(n, wave_speed, length) do
    n * wave_speed / (2 * length)
  end

  @doc "Harmonic frequencies of an open cylindrical pipe: fₙ = nv/(2L)."
  @spec open_pipe_harmonic(number, number, number) :: float
  def open_pipe_harmonic(n, wave_speed, length) do
    n * wave_speed / (2 * length)
  end

  @doc """
  Closed pipe (one open end) harmonics: f_n = (2n-1) v / (4 L).

  Only odd harmonics are present.

  ## Examples
      iex> Waves.closed_pipe_harmonic(1, 340, 0.25)
      340.0
  """
  @spec closed_pipe_harmonic(number, number, number) :: float
  def closed_pipe_harmonic(n, wave_speed, length) do
    (2 * n - 1) * wave_speed / (4 * length)
  end

  # ---------------------------------------------------------------------------
  # Phase & Group Velocity
  # ---------------------------------------------------------------------------

  @doc "Phase velocity of a wave: v_ph = omega/k."
  @spec phase_velocity(number, number) :: float
  def phase_velocity(omega, k), do: omega / k * 1.0

  @doc "Group velocity by finite difference: v_g = Δomega/Δk."
  @spec group_velocity(number, number, number, number) :: float
  def group_velocity(omega1, omega2, k1, k2) do
    (omega2 - omega1) / (k2 - k1)
  end

  @doc "Checks whether phase and group velocities are equal (non-dispersive medium)."
  @spec is_dispersionless?(number, number, number) :: boolean
  def is_dispersionless?(v_phase, v_group, tol \\ 1.0e-9) do
    abs(v_phase - v_group) < tol
  end

  # ---------------------------------------------------------------------------
  # Wave Power & Energy
  # ---------------------------------------------------------------------------

  @doc """
  Power of a transverse wave on a string: P = ½ μ omega² A² v

  ## Examples
      iex> Waves.wave_power(0.01, :math.pi(), 0.1, 340) > 0
      true
  """
  @spec wave_power(number, number, number, number) :: float
  def wave_power(mu, omega, amplitude, wave_speed) do
    0.5 * mu * omega ** 2 * amplitude ** 2 * wave_speed
  end

  @doc "Wave intensity (power per unit area): I = P/A."
  @spec wave_intensity(number, number) :: float
  def wave_intensity(power, area), do: power / area * 1.0

  @doc "Spherical wave intensity (inverse-square law): I = P/(4πr²)."
  @spec spherical_wave_intensity(number, number) :: float
  def spherical_wave_intensity(power, radius) do
    power / (4 * :math.pi() * radius ** 2)
  end

  # ---------------------------------------------------------------------------
  # Optics
  # ---------------------------------------------------------------------------

  @doc """
  Snell's Law: n₁ sin θ₁ = n₂ sin θ₂ — solves for θ₂.

  ## Examples
      iex> Waves.snells_law(1.0, 0.0, 1.5) |> Float.round(4)
      0.0
  """
  @spec snells_law(number, number, number) :: float
  def snells_law(n1, theta1, n2) do
    :math.asin(n1 * :math.sin(theta1) / n2)
  end

  @doc "Critical angle for total internal reflection: θ_c = arcsin(n₂/n₁)."
  @spec critical_angle(number, number) :: float
  def critical_angle(n1, n2), do: :math.asin(n2 / n1)

  @doc "Refractive index of a medium: n = c/v."
  @spec refractive_index(number, number) :: float
  def refractive_index(wave_speed, c \\ @speed_of_light), do: c / wave_speed

  @doc "Malus's Law for polarised light transmitted through a polariser: I = I₀ cos²θ."
  @spec malus_law(number, number) :: float
  def malus_law(i0, theta), do: i0 * :math.cos(theta) ** 2

  @doc """
  Thin lens equation: 1/f = 1/d_o + 1/d_i — solves for image distance d_i.

  ## Parameters
    - focal_length:    f (m), negative for diverging lens
    - object_distance: d_o (m, positive)

  ## Returns
    Image distance d_i (m); positive = real image, negative = virtual

  ## Examples
      iex> Waves.thin_lens_image_distance(0.1, 0.2) |> Float.round(4)
      0.2
  """
  @spec thin_lens_image_distance(number, number) :: float
  def thin_lens_image_distance(focal_length, object_distance) do
    1.0 / (1.0 / focal_length - 1.0 / object_distance)
  end

  @doc """
  Lateral magnification of a lens: m = -d_i / d_o.

  Negative value indicates an inverted real image.

  ## Examples
      iex> Waves.lens_magnification(0.2, 0.2)
      -1.0
  """
  @spec lens_magnification(number, number) :: float
  def lens_magnification(image_distance, object_distance) do
    -image_distance / object_distance
  end

  @doc """
  Mirror equation image distance: 1/f = 1/d_o + 1/d_i → d_i = f d_o / (d_o - f).

  Convention: f > 0 for concave (converging) mirrors.

  ## Examples
      iex> Waves.mirror_image_distance(0.5, 2.0) |> Float.round(4)
      0.6667
  """
  @spec mirror_image_distance(number, number) :: float
  def mirror_image_distance(focal_length, object_distance) do
    focal_length * object_distance / (object_distance - focal_length)
  end
end
