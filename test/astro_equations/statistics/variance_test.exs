defmodule AstroEquations.Statistics.VarianceTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Statistics.Variance

  @data [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]

  describe "population/1" do
    test "classic Pearson example: σ² = 4" do
      assert_in_delta Variance.population(@data), 4.0, 1.0e-10
    end

    test "constant dataset has zero variance" do
      assert_in_delta Variance.population([3.0, 3.0, 3.0]), 0.0, 1.0e-10
    end

    test "single element has zero variance" do
      assert_in_delta Variance.population([5.0]), 0.0, 1.0e-10
    end

    test "empty list returns error" do
      assert {:error, _} = Variance.population([])
    end
  end

  describe "sample/1" do
    test "Bessel correction: s² > σ²" do
      s2 = Variance.sample(@data)
      sig2 = Variance.population(@data)
      assert s2 > sig2
    end

    test "classic example: s² ≈ 4.571" do
      assert_in_delta Variance.sample(@data), 4.571, 0.001
    end

    test "single element returns error (need ≥ 2)" do
      assert {:error, _} = Variance.sample([5.0])
    end
  end

  describe "weighted/2" do
    test "equal weights gives same result as population variance" do
      weights = List.duplicate(1.0, length(@data))
      assert_in_delta Variance.weighted(@data, weights), Variance.population(@data), 1.0e-10
    end

    test "error on mismatched lengths" do
      assert {:error, _} = Variance.weighted([1, 2, 3], [1, 2])
    end

    test "error on empty list" do
      assert {:error, _} = Variance.weighted([], [])
    end
  end

  describe "standard_error/1" do
    test "SE = σ/√N" do
      se = Variance.standard_error(@data)
      n = length(@data)
      sigma = :math.sqrt(Variance.population(@data))
      assert_in_delta se * :math.sqrt(n), sigma, 1.0e-6
    end

    test "larger dataset → smaller standard error" do
      small = Variance.standard_error([1.0, 2.0, 3.0])
      large = Variance.standard_error([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
      assert large < small
    end
  end

  describe "welford_accumulate/1 and welford_finalize/1" do
    test "streaming gives same mean as direct" do
      acc = Variance.welford_accumulate(@data)
      {:ok, mean, _, _} = Variance.welford_finalize(acc)
      expected_mean = Enum.sum(@data) / length(@data)
      assert_in_delta mean, expected_mean, 1.0e-10
    end

    test "streaming variance matches batch" do
      acc = Variance.welford_accumulate(@data)
      {:ok, _mean, pop_var, samp_var} = Variance.welford_finalize(acc)
      assert_in_delta pop_var, Variance.population(@data), 1.0e-10
      assert_in_delta samp_var, Variance.sample(@data), 1.0e-10
    end

    test "empty accumulation returns error" do
      acc = Variance.welford_accumulate([])
      assert {:error, _} = Variance.welford_finalize(acc)
    end

    test "single element returns ok with zero variance" do
      acc = Variance.welford_accumulate([42.0])
      {:ok, mean, pop_var, _} = Variance.welford_finalize(acc)
      assert_in_delta mean, 42.0, 1.0e-10
      assert_in_delta pop_var, 0.0, 1.0e-10
    end
  end

  describe "covariance/2" do
    test "perfectly correlated: Cov(x, x) = Var(x)" do
      cov = Variance.covariance(@data, @data)
      var = Variance.sample(@data)
      assert_in_delta cov, var, 1.0e-10
    end

    test "anti-correlated: Cov(x, -x) = -Var(x)" do
      neg = Enum.map(@data, fn x -> -x end)
      cov = Variance.covariance(@data, neg)
      var = Variance.sample(@data)
      assert_in_delta cov, -var, 1.0e-10
    end

    test "error on unequal lengths" do
      assert {:error, _} = Variance.covariance([1, 2, 3], [1, 2])
    end
  end

  describe "pearson_r/2" do
    test "perfectly positively correlated: r = 1" do
      x = [1.0, 2.0, 3.0, 4.0, 5.0]
      y = [2.0, 4.0, 6.0, 8.0, 10.0]
      assert_in_delta Variance.pearson_r(x, y), 1.0, 1.0e-10
    end

    test "perfectly negatively correlated: r = -1" do
      x = [1.0, 2.0, 3.0, 4.0, 5.0]
      y = [-1.0, -2.0, -3.0, -4.0, -5.0]
      assert_in_delta Variance.pearson_r(x, y), -1.0, 1.0e-10
    end

    test "r is in [-1, 1]" do
      x = Enum.to_list(1..20) |> Enum.map(&(&1 * 1.0))
      y = Enum.map(x, fn v -> v + :rand.uniform() end)
      r = Variance.pearson_r(x, y)
      assert r >= -1.0 and r <= 1.0
    end
  end

  describe "chi_squared/2" do
    test "perfect fit: O = E → χ² = 0" do
      e = [10.0, 20.0, 30.0]
      assert_in_delta Variance.chi_squared(e, e), 0.0, 1.0e-10
    end

    test "always non-negative" do
      assert Variance.chi_squared([10, 20, 30], [12, 18, 30]) >= 0
    end

    test "error on unequal lengths" do
      assert {:error, _} = Variance.chi_squared([1, 2], [1, 2, 3])
    end
  end

  describe "reduced_chi_squared/2" do
    test "χ²_ν = χ²/ν" do
      chi2 = 9.0
      dof = 3
      assert_in_delta Variance.reduced_chi_squared(chi2, dof), 3.0, 1.0e-10
    end
  end

  describe "mad/1" do
    test "classic example: median of [2,4,4,4,5,5,7,9] = 4.5, deviations = [2.5,0.5,0.5,0.5,0.5,0.5,2.5,4.5], MAD = 0.5" do
      assert_in_delta Variance.mad(@data), 0.5, 1.0e-10
    end

    test "constant dataset has MAD = 0" do
      assert_in_delta Variance.mad([3.0, 3.0, 3.0, 3.0]), 0.0, 1.0e-10
    end
  end

  describe "iqr/1" do
    test "IQR > 0 for non-constant data" do
      assert Variance.iqr(@data) > 0
    end

    test "IQR = 0 for constant data" do
      assert_in_delta Variance.iqr([5.0, 5.0, 5.0, 5.0]), 0.0, 1.0e-10
    end
  end
end

defmodule AstroEquations.Statistics.StandardDeviationTest do
  use ExUnit.Case, async: true

  alias AstroEquations.Statistics.StandardDeviation, as: SD

  @data [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]

  describe "population/1" do
    test "σ = 2 for the classic Pearson dataset" do
      assert_in_delta SD.population(@data), 2.0, 1.0e-10
    end

    test "σ = √σ²" do
      sigma = SD.population(@data)
      assert_in_delta sigma * sigma, 4.0, 1.0e-10
    end

    test "empty list returns error" do
      assert {:error, _} = SD.population([])
    end
  end

  describe "sample/1" do
    test "s ≈ 2.138 for the classic dataset" do
      assert_in_delta SD.sample(@data), 2.1381, 0.0001
    end

    test "s > σ (Bessel correction)" do
      assert SD.sample(@data) > SD.population(@data)
    end
  end

  describe "standard_error/1" do
    test "SE = s/√N" do
      n = length(@data)
      s = SD.sample(@data)
      se = SD.standard_error(@data)
      assert_in_delta se, s / :math.sqrt(n), 1.0e-10
    end
  end

  describe "z_score/3" do
    test "mean is at z = 0" do
      mean = Enum.sum(@data) / length(@data)
      sigma = SD.population(@data)
      assert_in_delta SD.z_score(mean, mean, sigma), 0.0, 1.0e-10
    end

    test "one sigma above mean is z = 1" do
      mean = 0.0
      sigma = 1.0
      assert_in_delta SD.z_score(1.0, mean, sigma), 1.0, 1.0e-10
    end
  end

  describe "normalise/1" do
    test "normalised mean is 0" do
      z_scores = SD.normalise(@data)
      mean_z = Enum.sum(z_scores) / length(z_scores)
      assert_in_delta mean_z, 0.0, 1.0e-10
    end

    test "normalised standard deviation is 1" do
      z_scores = SD.normalise(@data)
      s = SD.sample(z_scores)
      assert_in_delta s, 1.0, 1.0e-8
    end

    test "empty list returns error" do
      assert {:error, _} = SD.normalise([])
    end
  end

  describe "propagate_addition/2" do
    test "σ_c = √(σ_a² + σ_b²): 3, 4 → 5" do
      assert_in_delta SD.propagate_addition(3.0, 4.0), 5.0, 1.0e-10
    end

    test "larger than either individual uncertainty" do
      sa = 0.05
      sb = 0.03
      sc = SD.propagate_addition(sa, sb)
      assert sc > sa
      assert sc > sb
    end
  end

  describe "propagate_product/4" do
    test "result is always positive" do
      assert SD.propagate_product(10.0, 0.5, 5.0, 0.2) > 0
    end

    test "relative uncertainty adds in quadrature" do
      a = 10.0
      sa = 1.0
      b = 5.0
      sb = 0.5
      sc = SD.propagate_product(a, sa, b, sb)
      rel_c = sc / (a * b)
      rel_a = sa / a
      rel_b = sb / b
      assert_in_delta rel_c, :math.sqrt(rel_a ** 2 + rel_b ** 2), 1.0e-10
    end
  end

  describe "propagate_power/3" do
    test "y = x² → σ_y = 2x σ_x" do
      x = 3.0
      sx = 0.1
      sigma_y = SD.propagate_power(x, sx, 2)
      assert_in_delta sigma_y, 2 * x * sx, 1.0e-10
    end
  end

  describe "propagate_log/2" do
    test "σ_y = σ_x / |x|" do
      x = 100.0
      sx = 5.0
      assert_in_delta SD.propagate_log(x, sx), sx / x, 1.0e-10
    end
  end

  describe "photometric_snr/4" do
    test "pure source dominated: SNR = √N_star" do
      n_star = 10_000.0
      snr = SD.photometric_snr(n_star, 0, 0, 0)
      assert_in_delta snr, :math.sqrt(n_star), 0.001
    end

    test "sky noise reduces SNR" do
      n_star = 10_000.0
      snr_no_sky = SD.photometric_snr(n_star, 0, 0, 0)
      snr_sky = SD.photometric_snr(n_star, 100, 5, 0)
      assert snr_sky < snr_no_sky
    end
  end

  describe "magnitude_uncertainty/1" do
    test "σ_m ≈ 1.0857 / SNR" do
      snr = 50.0
      assert_in_delta SD.magnitude_uncertainty(snr), 1.0857 / snr, 1.0e-4
    end

    test "higher SNR → smaller magnitude uncertainty" do
      assert SD.magnitude_uncertainty(100) < SD.magnitude_uncertainty(10)
    end
  end

  describe "poisson_uncertainty/1" do
    test "σ = √N: N = 100 → σ = 10" do
      assert_in_delta SD.poisson_uncertainty(100), 10.0, 1.0e-10
    end

    test "σ = √N: N = 10000 → σ = 100" do
      assert_in_delta SD.poisson_uncertainty(10_000), 100.0, 1.0e-10
    end
  end

  describe "online/1" do
    test "streaming gives correct mean and std dev" do
      data = [1.0, 2.0, 3.0, 4.0, 5.0]
      {:ok, mean, _pop_std, samp_std} = SD.online(data)
      assert_in_delta mean, 3.0, 1.0e-10
      assert_in_delta samp_std, SD.sample(data), 1.0e-10
    end

    test "works on large ranges" do
      {:ok, mean, _, _} = SD.online(1..1000)
      assert_in_delta mean, 500.5, 0.001
    end
  end
end
