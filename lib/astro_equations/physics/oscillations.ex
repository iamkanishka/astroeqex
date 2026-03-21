defmodule AstroEquations.Physics.Oscillations do
  @moduledoc """
  Oscillatory and wave mechanics.

  Covers:
  - Hooke's Law (spring force and potential energy)
  - Simple harmonic motion (SHM) — displacement, velocity, acceleration
  - Angular frequency, period, frequency
  - Energy in SHM
  - Damped oscillations (under-, over-, critically-damped)
  - Logarithmic decrement and settling time
  - Driven oscillations, resonance, and Q-factor
  - Pendula (simple and physical)
  - LC circuit oscillations and energy
  - Coupled oscillators (beat frequency)
  - Springs in series and parallel
  """

  # ---------------------------------------------------------------------------
  # Spring Force (Hooke's Law)
  # ---------------------------------------------------------------------------

  @doc """
  Spring restoring force: F = -k x

  ## Examples
      iex> Oscillations.force(10, 0.5)
      -5.0
  """
  @spec force(number, number) :: float
  def force(k_s, x), do: -k_s * x * 1.0

  @doc """
  Elastic potential energy: U = ½ k x²

  ## Examples
      iex> Oscillations.potential_energy(10, 0.5)
      1.25
  """
  @spec potential_energy(number, number) :: float
  def potential_energy(k_s, x), do: 0.5 * k_s * :math.pow(x, 2)

  @doc """
  Effective spring constant for two springs in series: 1/k_eff = 1/k₁ + 1/k₂

  ## Examples
      iex> Oscillations.springs_series(10.0, 10.0) |> Float.round(4)
      5.0
  """
  @spec springs_series(number, number) :: float
  def springs_series(k1, k2), do: k1 * k2 / (k1 + k2)

  @doc """
  Effective spring constant for two springs in parallel: k_eff = k₁ + k₂

  ## Examples
      iex> Oscillations.springs_parallel(10.0, 5.0) |> Float.round(4)
      15.0
  """
  @spec springs_parallel(number, number) :: float
  def springs_parallel(k1, k2), do: k1 + k2

  # ---------------------------------------------------------------------------
  # SHM Frequencies
  # ---------------------------------------------------------------------------

  @doc """
  Angular frequency of spring-mass system: omega = √(k / m)

  ## Examples
      iex> Oscillations.angular_frequency(10, 2.5)
      2.0
  """
  @spec angular_frequency(number, number) :: float
  def angular_frequency(k_s, m), do: :math.sqrt(k_s / m)

  @doc """
  Angular frequency of a simple pendulum (small angle): omega = √(g / L)

  ## Examples
      iex> Oscillations.pendulum_angular_frequency(1.0) |> Float.round(4)
      3.1321
  """
  @spec pendulum_angular_frequency(number, number) :: float
  def pendulum_angular_frequency(length, g \\ 9.81), do: :math.sqrt(g / length)

  @doc """
  Angular frequency of a physical (compound) pendulum: omega = √(m g d / I)

  ## Examples
      iex> Oscillations.physical_pendulum_frequency(1, 0.5, 0.5) |> Float.round(4)
      3.1321
  """
  @spec physical_pendulum_frequency(number, number, number, number) :: float
  def physical_pendulum_frequency(mass, d, moment_of_inertia, g \\ 9.81) do
    :math.sqrt(mass * g * d / moment_of_inertia)
  end

  @doc "Oscillation period from angular frequency: T = 2π/omega."
  @spec period_from_omega(number) :: float
  def period_from_omega(omega), do: 2 * :math.pi() / omega

  @doc "Period of a spring-mass system: T = 2π√(m/k)."
  @spec spring_period(number, number) :: float
  def spring_period(mass, k), do: 2 * :math.pi() * :math.sqrt(mass / k)

  @doc "Period of a simple pendulum (small angle): T = 2π√(L/g)."
  @spec pendulum_period(number, number) :: float
  def pendulum_period(length, g \\ 9.81), do: 2 * :math.pi() * :math.sqrt(length / g)

  @doc "Ordinary frequency from period: f = 1/T."
  @spec frequency_from_period(number) :: float
  def frequency_from_period(period), do: 1.0 / period

  @doc "Angular frequency from ordinary frequency: omega = 2πf."
  @spec omega_from_frequency(number) :: float
  def omega_from_frequency(f), do: 2 * :math.pi() * f

  # ---------------------------------------------------------------------------
  # SHM Kinematics
  # ---------------------------------------------------------------------------

  @doc "Simple harmonic displacement: x(t) = A cos(omegat + φ)."
  @spec shm_displacement(number, number, number, number) :: float
  def shm_displacement(amplitude, omega, t, phi \\ 0.0) do
    amplitude * :math.cos(omega * t + phi)
  end

  @doc "Simple harmonic velocity: v(t) = −Aomega sin(omegat + φ)."
  @spec shm_velocity(number, number, number, number) :: float
  def shm_velocity(amplitude, omega, t, phi \\ 0.0) do
    -amplitude * omega * :math.sin(omega * t + phi)
  end

  @doc "Simple harmonic acceleration: a(t) = -Aomega² cos(omegat + φ)."
  @spec shm_acceleration(number, number, number, number) :: float
  def shm_acceleration(amplitude, omega, t, phi \\ 0.0) do
    -amplitude * :math.pow(omega, 2) * :math.cos(omega * t + phi)
  end

  @doc "Instantaneous SHM speed from position: v = omega√(A² − x²)."
  @spec shm_velocity_from_position(number, number, number) :: float
  def shm_velocity_from_position(amplitude, omega, x) do
    omega * :math.sqrt(max(:math.pow(amplitude, 2) - :math.pow(x, 2), 0.0))
  end

  # ---------------------------------------------------------------------------
  # SHM Energy
  # ---------------------------------------------------------------------------

  @doc "Total mechanical energy in SHM (constant): E = ½kA²."
  @spec shm_total_energy(number, number) :: float
  def shm_total_energy(k, amplitude), do: 0.5 * k * :math.pow(amplitude, 2)

  @doc "Instantaneous kinetic energy in SHM: KE = ½momega²(A² − x²)."
  @spec shm_kinetic_energy(number, number, number, number) :: float
  def shm_kinetic_energy(mass, omega, amplitude, x) do
    0.5 * mass * :math.pow(omega, 2) * (:math.pow(amplitude, 2) - :math.pow(x, 2))
  end

  @doc """
  SHM energy at time t for a damped oscillator: E(t) = ½kA² e^(-bt/m)

  ## Parameters
    - k:         Spring constant (N/m)
    - amplitude: Initial amplitude A (m)
    - b:         Damping coefficient (kg/s)
    - mass:      Mass (kg)
    - t:         Time (s)

  ## Examples
      iex> Oscillations.shm_energy_at_time(10, 0.1, 0, 1.0, 5.0) |> Float.round(4)
      0.05
  """
  @spec shm_energy_at_time(number, number, number, number, number) :: float
  def shm_energy_at_time(k, amplitude, b, mass, t) do
    0.5 * k * :math.pow(amplitude, 2) * :math.exp(-b * t / mass)
  end

  # ---------------------------------------------------------------------------
  # Damped Oscillations
  # ---------------------------------------------------------------------------

  @doc """
  Damping ratio: ζ = b / (2 √(k m))

  ## Examples
      iex> Oscillations.damping_ratio(2, 10, 1.0) |> Float.round(4)
      0.3162
  """
  @spec damping_ratio(number, number, number) :: float
  def damping_ratio(b, k, mass), do: b / (2 * :math.sqrt(k * mass))

  @doc """
  Damped natural frequency: omega_d = omega₀ √(1 - ζ²)  (valid for ζ < 1).
  """
  @spec damped_frequency(number, number) :: float
  def damped_frequency(omega_0, zeta) when zeta < 1.0 do
    omega_0 * :math.sqrt(1 - zeta * zeta)
  end

  @doc """
  Underdamped displacement: x(t) = A e^(-ζomega₀t) cos(omega_d t + φ)

  ## Examples
      iex> Oscillations.underdamped_displacement(0.1, 0, 10, 10, 0) |> Float.round(4)
      0.1
  """
  @spec underdamped_displacement(number, number, number, number, number, number) :: float
  def underdamped_displacement(amplitude, zeta, omega_0, omega_d, t, phi \\ 0.0) do
    amplitude *
      :math.exp(-zeta * omega_0 * t) *
      :math.cos(omega_d * t + phi)
  end

  @doc """
  Quality factor (Q-factor): Q = omega₀ m / b

  ## Examples
      iex> Oscillations.quality_factor(10, 1.0, 0.5) |> Float.round(4)
      20.0
  """
  @spec quality_factor(number, number, number) :: float
  def quality_factor(omega_0, mass, b), do: omega_0 * mass / b

  @doc """
  Logarithmic decrement: δ = π b / (m omega_d) = 2π ζ / √(1 - ζ²)

  The natural log of the ratio of successive amplitude peaks in an underdamped oscillator.
  δ = ln(x_n / x_{n+1}).

  ## Parameters
    - zeta:    Damping ratio
    - omega_0: Natural angular frequency (rad/s) — used for normalisation if needed
    - omega_d: Damped angular frequency (rad/s)

  ## Examples
      iex> Oscillations.log_decrement(0.1, 10.0, 9.95) |> Float.round(4)
      0.6315
  """
  @spec log_decrement(number, number, number) :: float
  def log_decrement(zeta, omega_0, omega_d) do
    :math.pi() * zeta * omega_0 / omega_d
  end

  @doc """
  Approximate settling time to within 2% of final value:
  t_s ≈ 4 / (ζ omega₀)

  Valid for underdamped systems (ζ < 1).

  ## Parameters
    - zeta:    Damping ratio
    - omega_0: Natural angular frequency (rad/s)

  ## Examples
      iex> Oscillations.settling_time(0.1, 10.0) |> Float.round(4)
      4.0
  """
  @spec settling_time(number, number) :: float
  def settling_time(zeta, omega_0), do: 4.0 / (zeta * omega_0)

  # ---------------------------------------------------------------------------
  # Driven Oscillations
  # ---------------------------------------------------------------------------

  @doc """
  Resonance amplitude of a driven damped oscillator:
  A(omega) = F₀/m / √((omega₀² - omega²)² + (b omega/m)²)

  ## Examples
      iex> Oscillations.resonance_amplitude(1, 1, 10, 0, 0.1) > 0
      true
  """
  @spec resonance_amplitude(number, number, number, number, number) :: float
  def resonance_amplitude(f0, mass, omega_0, omega, b) do
    f0 / mass /
      :math.sqrt(
        :math.pow(omega_0 * omega_0 - omega * omega, 2) +
          :math.pow(b * omega / mass, 2)
      )
  end

  @doc """
  Resonant frequency of a driven oscillator (maximum amplitude):
  omega_res = √(omega₀² - b²/(2m²))

  ## Examples
      iex> Oscillations.resonant_frequency(10, 0, 1) |> Float.round(4)
      10.0
  """
  @spec resonant_frequency(number, number, number) :: float
  def resonant_frequency(omega_0, b, mass) do
    arg = omega_0 * omega_0 - b * b / (2 * mass * mass)
    if arg > 0, do: :math.sqrt(arg), else: 0.0
  end

  @doc """
  Half-power (3 dB) bandwidth of a resonance peak: BW = omega₀ / Q = b / m

  ## Parameters
    - omega_0: Natural angular frequency (rad/s)
    - q:       Q-factor

  ## Examples
      iex> Oscillations.half_power_bandwidth(10.0, 20.0) |> Float.round(4)
      0.5
  """
  @spec half_power_bandwidth(number, number) :: float
  def half_power_bandwidth(omega_0, q), do: omega_0 / q

  # ---------------------------------------------------------------------------
  # LC Circuit
  # ---------------------------------------------------------------------------

  @doc """
  LC circuit oscillation angular frequency: omega = 1 / √(LC)

  ## Examples
      iex> Oscillations.lc_angular_frequency(1.0e-3, 10.0e-6) |> Float.round(0)
      10000.0
  """
  @spec lc_angular_frequency(number, number) :: float
  def lc_angular_frequency(inductance, capacitance) do
    1.0 / :math.sqrt(inductance * capacitance)
  end

  @doc """
  Total energy stored in an LC circuit at peak current: E = ½ L I_max²

  Equivalently E = ½ C V_max² at peak voltage.

  ## Parameters
    - inductance:    L (H)
    - peak_current:  I_max (A)

  ## Examples
      iex> Oscillations.lc_energy(1.0e-3, 2.0) |> Float.round(6)
      0.002
  """
  @spec lc_energy(number, number) :: float
  def lc_energy(inductance, peak_current) do
    0.5 * inductance * peak_current ** 2
  end

  # ---------------------------------------------------------------------------
  # Coupled Oscillators / Beats
  # ---------------------------------------------------------------------------

  @doc """
  Beat frequency between two oscillators: f_beat = |f₁ - f₂|

  ## Examples
      iex> Oscillations.beat_frequency(440, 444)
      4.0
  """
  @spec beat_frequency(number, number) :: float
  def beat_frequency(f1, f2), do: abs(f1 - f2) * 1.0
end
