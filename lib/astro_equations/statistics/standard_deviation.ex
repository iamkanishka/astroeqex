defmodule AstroEquations.Statistics.StandardDeviation do
  @moduledoc """
  Standard deviation and related statistics for astronomical data analysis.

  Provides:
  - Population and sample standard deviation
  - Weighted standard deviation
  - Running (online) standard deviation
  - Normalised (z-score) values
  - Propagation of uncertainty (error propagation)
  - Photometric error estimates (Poisson, read-noise limited)
  - Signal-to-noise ratio from count data
  """

  alias AstroEquations.Statistics.Variance

  # ---------------------------------------------------------------------------
  # Standard Deviation
  # ---------------------------------------------------------------------------

  @doc """
  Population standard deviation: σ = sqrt(σ²)

  ## Parameters
    - data: List of numbers

  ## Returns
    Population standard deviation

  ## Examples
      iex> StandardDeviation.population([2, 4, 4, 4, 5, 5, 7, 9])
      2.0
  """
  @spec population([number]) :: float | {:error, String.t()}
  def population(data) do
    case Variance.population(data) do
      {:error, _} = e -> e
      v -> :math.sqrt(v)
    end
  end

  @doc """
  Sample standard deviation: s = sqrt(s²).

  Uses Bessel's correction (N-1 denominator).

  ## Parameters
    - data: List of numbers (at least 2 elements)

  ## Returns
    Sample standard deviation

  ## Examples
      iex> StandardDeviation.sample([2, 4, 4, 4, 5, 5, 7, 9]) |> Float.round(4)
      2.1381
  """
  @spec sample([number]) :: float | {:error, String.t()}
  def sample(data) do
    case Variance.sample(data) do
      {:error, _} = e -> e
      v -> :math.sqrt(v)
    end
  end

  @doc """
  Weighted standard deviation.

  ## Parameters
    - data:    List of values
    - weights: List of non-negative weights

  ## Returns
    Weighted standard deviation

  ## Examples
      iex> StandardDeviation.weighted([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]) |> Float.round(4)
      0.8165
  """
  @spec weighted([number], [number]) :: float | {:error, String.t()}
  def weighted(data, weights) do
    case Variance.weighted(data, weights) do
      {:error, _} = e -> e
      v -> :math.sqrt(v)
    end
  end

  @doc """
  Standard error of the mean: SE = s / √N

  ## Parameters
    - data: List of numbers

  ## Returns
    Standard error

  ## Examples
      iex> StandardDeviation.standard_error([1.0, 2.0, 3.0]) > 0
      true
  """
  @spec standard_error([number]) :: float | {:error, String.t()}
  def standard_error(data), do: Variance.standard_error(data)

  # ---------------------------------------------------------------------------
  # Z-Score Normalisation
  # ---------------------------------------------------------------------------

  @doc """
  Z-score of a single value: z = (x - μ) / σ

  ## Parameters
    - x:     Value to normalise
    - mean:  Distribution mean μ
    - sigma: Standard deviation σ

  ## Returns
    Z-score

  ## Examples
      iex> StandardDeviation.z_score(5.0, 5.0, 1.0) |> Float.round(4)
      0.0
  """
  @spec z_score(number, number, number) :: float
  def z_score(x, mean, sigma), do: (x - mean) / sigma

  @doc """
  Normalises a list to zero mean and unit standard deviation (z-score normalisation).

  ## Parameters
    - data: List of numbers

  ## Returns
    List of z-scores, or {:error, reason}
  """
  @spec normalise([number]) :: [float] | {:error, String.t()}
  def normalise([]), do: {:error, "Empty dataset"}

  def normalise(data) do
    n = length(data)
    mean = Enum.sum(data) / n

    case sample(data) do
      {:error, _} = e -> e
      sigma -> Enum.map(data, fn x -> (x - mean) / sigma end)
    end
  end

  @doc """
  Coefficient of variation: CV = σ / |μ| × 100%

  Useful for comparing variability across datasets with different units or scales.

  ## Parameters
    - data: List of numbers

  ## Returns
    CV as a percentage, or {:error, reason}

  ## Examples
      iex> StandardDeviation.coefficient_of_variation([2, 4, 4, 4, 5, 5, 7, 9]) > 0
      true
  """
  @spec coefficient_of_variation([number]) :: float | {:error, String.t()}
  def coefficient_of_variation([]), do: {:error, "Empty dataset"}

  def coefficient_of_variation(data) do
    mean = Enum.sum(data) / length(data)

    if mean == 0 do
      {:error, "Mean is zero; CV undefined"}
    else
      case sample(data) do
        {:error, _} = e -> e
        sigma -> sigma / abs(mean) * 100.0
      end
    end
  end

  @doc """
  Pooled standard deviation for two groups of equal size:
  s_p = √((s₁² + s₂²) / 2)

  Used when combining uncertainty estimates from independent sub-samples.

  ## Parameters
    - sigma1: Standard deviation of group 1
    - sigma2: Standard deviation of group 2

  ## Returns
    Pooled standard deviation

  ## Examples
      iex> StandardDeviation.pooled_std(3.0, 4.0) |> Float.round(4)
      3.5355
  """
  @spec pooled_std(number, number) :: float
  def pooled_std(sigma1, sigma2) do
    :math.sqrt((sigma1 * sigma1 + sigma2 * sigma2) / 2)
  end

  # ---------------------------------------------------------------------------
  # Error Propagation
  # ---------------------------------------------------------------------------

  @doc """
  Propagated uncertainty for addition or subtraction: σ_c = sqrt(σ_a² + σ_b²)

  ## Parameters
    - sigma_a: Uncertainty in a
    - sigma_b: Uncertainty in b

  ## Returns
    Combined uncertainty

  ## Examples
      iex> StandardDeviation.propagate_addition(3.0, 4.0) |> Float.round(4)
      5.0
  """
  @spec propagate_addition(number, number) :: float
  def propagate_addition(sigma_a, sigma_b), do: :math.sqrt(sigma_a ** 2 + sigma_b ** 2)

  @doc """
  Propagated relative uncertainty for multiplication or division:
  (σ_c/c)² = (σ_a/a)² + (σ_b/b)²

  ## Parameters
    - a, sigma_a: Value and its uncertainty
    - b, sigma_b: Value and its uncertainty

  ## Returns
    Absolute uncertainty in c = a * b (or a / b)

  ## Examples
      iex> StandardDeviation.propagate_product(10, 1, 5, 0.5) |> Float.round(4)
      5.5902
  """
  @spec propagate_product(number, number, number, number) :: float
  def propagate_product(a, sigma_a, b, sigma_b) do
    abs(a * b) * :math.sqrt((sigma_a / a) ** 2 + (sigma_b / b) ** 2)
  end

  @doc """
  Propagated uncertainty for a power law: y = x^n → σ_y = |n| |x|^(n-1) σ_x

  ## Parameters
    - x:       Value
    - sigma_x: Uncertainty in x
    - n:       Exponent

  ## Returns
    Uncertainty in y = x^n

  ## Examples
      iex> StandardDeviation.propagate_power(3.0, 0.1, 2) |> Float.round(4)
      0.6
  """
  @spec propagate_power(number, number, number) :: float
  def propagate_power(x, sigma_x, n), do: abs(n) * abs(:math.pow(x, n - 1)) * sigma_x

  @doc """
  Propagated uncertainty through a natural logarithm: y = ln(x) → σ_y = σ_x / |x|

  ## Parameters
    - x:       Value (must be > 0)
    - sigma_x: Uncertainty in x

  ## Returns
    Uncertainty in y = ln(x)

  ## Examples
      iex> StandardDeviation.propagate_log(100.0, 5.0) |> Float.round(4)
      0.05
  """
  @spec propagate_log(number, number) :: float
  def propagate_log(x, sigma_x), do: sigma_x / abs(x)

  @doc """
  Propagated uncertainty for general function f(x): σ_y ≈ |df/dx| σ_x
  using numerical differentiation.

  ## Parameters
    - f:       Function f(x)
    - x:       Value at which to evaluate
    - sigma_x: Uncertainty in x
    - h:       Step size (default: 1.0e-7)

  ## Returns
    Propagated uncertainty σ_y

  ## Examples
      iex> StandardDeviation.propagate_function(&:math.sqrt/1, 4.0, 0.1) |> Float.round(4)
      0.025
  """
  @spec propagate_function((number -> number), number, number, number) :: float
  def propagate_function(f, x, sigma_x, h \\ 1.0e-7) do
    df_dx = (f.(x + h) - f.(x - h)) / (2 * h)
    abs(df_dx) * sigma_x
  end

  # ---------------------------------------------------------------------------
  # Astronomical / Photometric Errors
  # ---------------------------------------------------------------------------

  @doc """
  Photometric signal-to-noise ratio in the Poisson (shot-noise) limit.

  SNR = N_star / sqrt(N_star + N_sky + N_dark + N_read²)

  where all N values are in electrons.

  ## Parameters
    - n_star:  Source electrons
    - n_sky:   Sky background electrons (total in aperture)
    - n_dark:  Dark current electrons (total in aperture)
    - n_read:  Read noise in electrons (per pixel, square-root of n_pixels applied)

  ## Returns
    SNR

  ## Examples
      iex> StandardDeviation.photometric_snr(10_000, 0, 0, 0) |> Float.round(2)
      100.0
  """
  @spec photometric_snr(number, number, number, number) :: float
  def photometric_snr(n_star, n_sky, n_dark, n_read) do
    n_star / :math.sqrt(n_star + n_sky + n_dark + n_read ** 2)
  end

  @doc """
  Photometric magnitude uncertainty from SNR.

  σ_m ≈ 2.5 / (ln(10) × SNR) ≈ 1.0857 / SNR

  ## Parameters
    - snr: Signal-to-noise ratio

  ## Returns
    Magnitude uncertainty in magnitudes

  ## Examples
      iex> StandardDeviation.magnitude_uncertainty(50) |> Float.round(4)
      0.0217
  """
  @spec magnitude_uncertainty(number) :: float
  def magnitude_uncertainty(snr), do: 1.0857 / snr

  @doc """
  Minimum detectable source flux in the sky-background-limited case.

  f_lim = SNR_threshold × √(N_sky × t_exp) / t_exp

  ## Parameters
    - snr_threshold:  Required SNR
    - sky_electrons:  Sky electrons per second in aperture
    - exposure_time:  Integration time (s)

  ## Returns
    Minimum source flux in electrons/s

  ## Examples
      iex> StandardDeviation.limiting_flux(5, 100, 60) > 0
      true
  """
  @spec limiting_flux(number, number, number) :: float
  def limiting_flux(snr_threshold, sky_electrons, exposure_time) do
    snr_threshold * :math.sqrt(sky_electrons * exposure_time) / exposure_time
  end

  @doc """
  Poisson uncertainty on a count: σ = √N

  ## Parameters
    - n: Observed count

  ## Returns
    Uncertainty (√N)

  ## Examples
      iex> StandardDeviation.poisson_uncertainty(100) |> Float.round(4)
      10.0
  """
  @spec poisson_uncertainty(number) :: float
  def poisson_uncertainty(n), do: :math.sqrt(n)

  # ---------------------------------------------------------------------------
  # Running / Online Standard Deviation
  # ---------------------------------------------------------------------------

  @doc """
  Computes mean and standard deviation in a single streaming pass using
  Welford's online algorithm.

  ## Parameters
    - data: Enumerable of numbers

  ## Returns
    {:ok, mean, population_std, sample_std}, or {:error, reason}

  ## Examples
      iex> StandardDeviation.online([1.0, 2.0, 3.0]) |> elem(0)
      :ok
  """
  @spec online(Enumerable.t()) :: {:ok, float, float, float} | {:error, String.t()}
  def online(data) do
    acc = Variance.welford_accumulate(data)

    case Variance.welford_finalize(acc) do
      {:error, _} = e ->
        e

      {:ok, mean, pop_var, samp_var} ->
        {:ok, mean, :math.sqrt(pop_var), :math.sqrt(samp_var)}
    end
  end
end
