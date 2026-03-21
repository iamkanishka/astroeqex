defmodule AstroEquations.Statistics.Variance do
  @moduledoc """
  Variance and related dispersion statistics for astronomical data analysis.

  Provides:
  - Population and sample variance
  - Weighted variance
  - Biased/unbiased variance estimators
  - Variance of the mean (standard error of the mean)
  - Running / online variance (Welford's algorithm)
  - Covariance and correlation coefficient (Pearson, Spearman)
  - Chi-squared goodness-of-fit statistic
  - Reduced chi-squared
  - Median absolute deviation (MAD)
  - Interquartile range (IQR)
  - RMS value
  """

  # ---------------------------------------------------------------------------
  # Basic Variance
  # ---------------------------------------------------------------------------

  @doc """
  Population variance: σ² = Σ(xᵢ - μ)² / N

  ## Parameters
    - data: List of numbers

  ## Returns
    Population variance, or {:error, reason}

  ## Examples
      iex> Variance.population([2, 4, 4, 4, 5, 5, 7, 9])
      4.0
  """
  @spec population([number]) :: float | {:error, String.t()}
  def population([]), do: {:error, "Empty dataset"}

  def population(data) do
    n = length(data)
    mean = Enum.sum(data) / n
    Enum.sum(Enum.map(data, fn x -> :math.pow(x - mean, 2) end)) / n
  end

  @doc """
  Sample variance (unbiased): s² = Σ(xᵢ - x̄)² / (N - 1)

  Bessel's correction applied — use for samples from a population.

  ## Parameters
    - data: List of numbers (at least 2 elements)

  ## Returns
    Sample variance, or {:error, reason}

  ## Examples
      iex> Variance.sample([2, 4, 4, 4, 5, 5, 7, 9])
      4.571428571428571
  """
  @spec sample([number]) :: float | {:error, String.t()}
  def sample([]), do: {:error, "Empty dataset"}
  def sample([_]), do: {:error, "Need at least 2 data points for sample variance"}

  def sample(data) do
    n = length(data)
    mean = Enum.sum(data) / n
    Enum.sum(Enum.map(data, fn x -> :math.pow(x - mean, 2) end)) / (n - 1)
  end

  # ---------------------------------------------------------------------------
  # Weighted Variance
  # ---------------------------------------------------------------------------

  @doc """
  Weighted population variance: σ²_w = Σ wᵢ (xᵢ - μ_w)² / Σ wᵢ

  ## Parameters
    - data:    List of values
    - weights: List of non-negative weights (same length as data)

  ## Returns
    Weighted variance, or {:error, reason}

  ## Examples
      iex> Variance.weighted([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]) |> Float.round(4)
      0.6667
  """
  @spec weighted([number], [number]) :: float | {:error, String.t()}
  def weighted([], _), do: {:error, "Empty dataset"}
  def weighted(_, []), do: {:error, "Empty weights"}

  def weighted(data, weights) when length(data) != length(weights),
    do: {:error, "Data and weights must have the same length"}

  def weighted(data, weights) do
    sum_w = Enum.sum(weights)
    w_mean = Enum.zip(data, weights) |> Enum.reduce(0.0, fn {x, w}, acc -> acc + w * x end)
    mean = w_mean / sum_w

    w_sq_dev =
      Enum.zip(data, weights)
      |> Enum.reduce(0.0, fn {x, w}, acc -> acc + w * :math.pow(x - mean, 2) end)

    w_sq_dev / sum_w
  end

  # ---------------------------------------------------------------------------
  # Standard Error
  # ---------------------------------------------------------------------------

  @doc """
  Standard error of the mean: SE = s / √N

  ## Parameters
    - data: List of numbers

  ## Returns
    Standard error of the mean, or {:error, reason}
  """
  @spec standard_error([number]) :: float | {:error, String.t()}
  def standard_error([]), do: {:error, "Empty dataset"}

  def standard_error(data) do
    case sample(data) do
      {:error, _} = err -> err
      s2 -> :math.sqrt(s2) / :math.sqrt(length(data))
    end
  end

  # ---------------------------------------------------------------------------
  # Welford Online Algorithm
  # ---------------------------------------------------------------------------

  @doc """
  Online (streaming) variance using Welford's algorithm.

  Returns {count, mean, M2} accumulator; call `welford_finalize/1` to get variance.

  ## Parameters
    - data: Enumerable of numbers

  ## Returns
    {n, mean, M2} accumulator tuple

  ## Examples
      iex> Variance.welford_accumulate([1.0, 2.0, 3.0]) |> elem(0)
      3
  """
  @spec welford_accumulate(Enumerable.t()) :: {non_neg_integer, float, float}
  def welford_accumulate(data) do
    Enum.reduce(data, {0, 0.0, 0.0}, fn x, {n, mean, m2} ->
      n1 = n + 1
      delta = x - mean
      new_mean = mean + delta / n1
      delta2 = x - new_mean
      {n1, new_mean, m2 + delta * delta2}
    end)
  end

  @doc """
  Finalises Welford accumulator, returning {:ok, mean, population_variance, sample_variance}.

  ## Parameters
    - acc: {n, mean, M2} from `welford_accumulate/1`

  ## Examples
      iex> Variance.welford_finalize({3, 2.0, 2.0}) |> elem(0)
      :ok
  """
  @spec welford_finalize({non_neg_integer, float, float}) ::
          {:ok, float, float, float} | {:error, String.t()}
  def welford_finalize({0, _, _}), do: {:error, "No data"}
  def welford_finalize({1, mean, _}), do: {:ok, mean, 0.0, 0.0}

  def welford_finalize({n, mean, m2}) do
    {:ok, mean, m2 / n, m2 / (n - 1)}
  end

  # ---------------------------------------------------------------------------
  # Covariance & Correlation
  # ---------------------------------------------------------------------------

  @doc """
  Sample covariance between two variables: Cov(X, Y) = Σ(xᵢ - x̄)(yᵢ - ȳ) / (N - 1)

  ## Parameters
    - xs: List of x values
    - ys: List of y values (same length)

  ## Returns
    Sample covariance, or {:error, reason}
  """
  @spec covariance([number], [number]) :: float | {:error, String.t()}
  def covariance([], _), do: {:error, "Empty dataset"}

  def covariance(xs, ys) when length(xs) != length(ys),
    do: {:error, "Lists must have equal length"}

  def covariance(xs, ys) do
    n = length(xs)
    mx = Enum.sum(xs) / n
    my = Enum.sum(ys) / n

    Enum.zip(xs, ys)
    |> Enum.reduce(0.0, fn {x, y}, acc -> acc + (x - mx) * (y - my) end)
    |> Kernel./(n - 1)
  end

  @doc """
  Pearson correlation coefficient: r = Cov(X,Y) / (σ_X σ_Y)

  ## Parameters
    - xs: List of x values
    - ys: List of y values

  ## Returns
    Pearson r in [-1, 1], or {:error, reason}

  ## Examples
      iex> Variance.pearson_r([1, 2, 3], [1, 2, 3]) |> Float.round(4)
      1.0
  """
  @spec pearson_r([number], [number]) :: float | {:error, String.t()}
  def pearson_r(xs, ys) do
    with cov when is_float(cov) <- covariance(xs, ys),
         sx2 when is_float(sx2) <- sample(xs),
         sy2 when is_float(sy2) <- sample(ys) do
      cov / (:math.sqrt(sx2) * :math.sqrt(sy2))
    end
  end

  @doc """
  Spearman rank correlation coefficient.

  A non-parametric measure of rank-order association: rₛ = 1 - 6Σd²/(n(n²-1))

  ## Parameters
    - xs: List of x values
    - ys: List of y values (same length)

  ## Returns
    Spearman r in [-1, 1], or {:error, reason}

  ## Examples
      iex> Variance.spearman_r([1, 2, 3, 4], [1, 2, 3, 4]) |> Float.round(4)
      1.0
  """
  @spec spearman_r([number], [number]) :: float | {:error, String.t()}
  def spearman_r(xs, ys) when length(xs) != length(ys),
    do: {:error, "Lists must have equal length"}

  def spearman_r([], _), do: {:error, "Empty dataset"}

  def spearman_r(xs, ys) do
    rank = fn lst ->
      sorted = Enum.sort(lst)
      Enum.map(lst, fn x -> Enum.find_index(sorted, &(&1 == x)) + 1.0 end)
    end

    pearson_r(rank.(xs), rank.(ys))
  end

  # ---------------------------------------------------------------------------
  # Chi-Squared Statistics
  # ---------------------------------------------------------------------------

  @doc """
  Chi-squared goodness-of-fit statistic: χ² = Σ (O - E)² / E

  ## Parameters
    - observed: List of observed counts
    - expected: List of expected counts

  ## Returns
    χ² value, or {:error, reason}
  """
  @spec chi_squared([number], [number]) :: float | {:error, String.t()}
  def chi_squared([], _), do: {:error, "Empty dataset"}

  def chi_squared(observed, expected) when length(observed) != length(expected),
    do: {:error, "Lists must have equal length"}

  def chi_squared(observed, expected) do
    Enum.zip(observed, expected)
    |> Enum.reduce(0.0, fn {o, e}, acc -> acc + :math.pow(o - e, 2) / e end)
  end

  @doc """
  Reduced chi-squared: χ²_ν = χ² / ν  (ν = degrees of freedom).

  A value near 1.0 indicates a good fit.

  ## Parameters
    - chi2: Chi-squared statistic
    - dof:  Degrees of freedom ν

  ## Returns
    Reduced chi-squared

  ## Examples
      iex> Variance.reduced_chi_squared(6.0, 3) |> Float.round(4)
      2.0
  """
  @spec reduced_chi_squared(number, number) :: float
  def reduced_chi_squared(chi2, dof), do: chi2 / dof

  # ---------------------------------------------------------------------------
  # Robust Dispersion Estimators
  # ---------------------------------------------------------------------------

  @doc """
  Median absolute deviation (MAD): median |xᵢ - median(x)|

  Robust alternative to standard deviation for data with outliers.

  ## Parameters
    - data: List of numbers

  ## Returns
    MAD value, or {:error, reason}
  """
  @spec mad([number]) :: float | {:error, String.t()}
  def mad([]), do: {:error, "Empty dataset"}

  def mad(data) do
    med = median(data)
    deviations = Enum.map(data, fn x -> abs(x - med) end)
    median(deviations)
  end

  @doc """
  Interquartile range: IQR = Q3 - Q1

  ## Parameters
    - data: List of numbers

  ## Returns
    IQR value, or {:error, reason}
  """
  @spec iqr([number]) :: float | {:error, String.t()}
  def iqr([]), do: {:error, "Empty dataset"}

  def iqr(data) do
    sorted = Enum.sort(data)
    q1 = percentile(sorted, 25)
    q3 = percentile(sorted, 75)
    q3 - q1
  end

  @doc """
  Root mean square (RMS) of a dataset: RMS = √(Σxᵢ²/N)

  Useful in signal processing and photometry to quantify noise levels.

  ## Parameters
    - data: List of numbers

  ## Returns
    RMS value, or {:error, reason}

  ## Examples
      iex> Variance.rms([3.0, 4.0]) |> Float.round(4)
      3.5355
  """
  @spec rms([number]) :: float | {:error, String.t()}
  def rms([]), do: {:error, "Empty dataset"}

  def rms(data) do
    n = length(data)
    :math.sqrt(Enum.sum(Enum.map(data, fn x -> x * x end)) / n)
  end

  # ---------------------------------------------------------------------------
  # Private helpers
  # ---------------------------------------------------------------------------

  defp median(data) do
    sorted = Enum.sort(data)
    n = length(sorted)
    mid = div(n, 2)

    if rem(n, 2) == 0,
      do: (Enum.at(sorted, mid - 1) + Enum.at(sorted, mid)) / 2.0,
      else: Enum.at(sorted, mid) * 1.0
  end

  defp percentile(sorted, p) do
    n = length(sorted)
    idx = p / 100 * (n - 1)
    lo = trunc(idx)
    hi = lo + 1
    frac = idx - lo

    v_lo = Enum.at(sorted, lo)
    v_hi = Enum.at(sorted, min(hi, n - 1))
    v_lo + frac * (v_hi - v_lo)
  end
end
