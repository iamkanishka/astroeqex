defmodule AstroEquations.Mathematics.Trigonometry do
  @moduledoc """
  Trigonometric utilities tailored for astronomical calculations.

  Provides:
  - Spherical law of cosines and sines
  - Hour angle, altitude, and azimuth conversions
  - Parallactic angle
  - Refraction correction (Bennett)
  - Ecliptic ↔ equatorial coordinate conversion
  - Galactic ↔ equatorial coordinate conversion
  - Sidereal time (GMST, LST)
  - Sunrise / sunset / twilight hour angles
  - Approximate solar declination
  """

  # ---------------------------------------------------------------------------
  # Types
  # ---------------------------------------------------------------------------

  @typedoc "Angle in radians."
  @type angle :: float()

  @typedoc "Right ascension in radians [0, 2π)."
  @type right_ascension :: float()

  @typedoc "Declination in radians (−π/2 to π/2)."
  @type declination :: float()

  @typedoc "Julian Date (continuous count of days from epoch)."
  @type julian_date :: float()

  @pi :math.pi()
  @deg_to_rad @pi / 180.0
  # IAU J2000 galactic coordinate pole constants
  # RA of north galactic pole (rad)
  @alpha_ngp 3.366_032_942_76
  # Dec of north galactic pole (rad)
  @delta_ngp 0.473_477_222_9
  # l of north celestial pole (rad)
  @l_ncp 2.145_568_156_1

  defp to_rad(deg), do: deg * @deg_to_rad

  # ---------------------------------------------------------------------------
  # Spherical Trigonometry
  # ---------------------------------------------------------------------------

  @doc """
  Spherical law of cosines: cos(c) = cos(a) cos(b) + sin(a) sin(b) cos(C).

  Solves for side c given sides a, b and included angle C (all in radians).

  ## Examples
      iex> Trigonometry.spherical_law_of_cosines(:math.pi()/3, :math.pi()/3, :math.pi()/3) |> Float.round(4)
      1.0472
  """

  @spec spherical_law_of_cosines(number, number, number) :: float
  def spherical_law_of_cosines(a, b, big_c) do
    :math.acos(
      :math.cos(a) * :math.cos(b) +
        :math.sin(a) * :math.sin(b) * :math.cos(big_c)
    )
  end

  @doc """
  Spherical law of sines: sin(A)/sin(a) = sin(B)/sin(b).

  Solves for angle A given angle B and the opposite sides a, b.

  ## Examples
      iex> Trigonometry.spherical_law_of_sines(:math.pi()/3, :math.pi()/3, :math.pi()/4) |> Float.round(4)
      0.6847
  """

  @spec spherical_law_of_sines(number, number, number) :: float
  def spherical_law_of_sines(big_b, b, a) do
    :math.asin(:math.sin(big_b) * :math.sin(a) / :math.sin(b))
  end

  # ---------------------------------------------------------------------------
  # Altitude & Azimuth
  # ---------------------------------------------------------------------------

  @doc """
  Altitude of a celestial object: sin(alt) = sin(δ) sin(φ) + cos(δ) cos(φ) cos(H)

  ## Parameters
    - dec:        Declination δ (radians)
    - latitude:   Observer latitude φ (radians)
    - hour_angle: Local hour angle H (radians)

  ## Returns
    Altitude in radians

  ## Examples
      iex> Trigonometry.altitude(0.0, 0.0, 0.0) |> Float.round(4)
      1.5708
  """

  @spec altitude(number, number, number) :: float
  def altitude(dec, latitude, hour_angle)
      when is_number(dec) and is_number(latitude) and is_number(hour_angle) do
    :math.asin(
      :math.sin(dec) * :math.sin(latitude) +
        :math.cos(dec) * :math.cos(latitude) * :math.cos(hour_angle)
    )
  end

  @doc """
  Azimuth of a celestial object (measured north through east).

  ## Returns
    Azimuth in radians [0, 2π)

  ## Examples
      iex> Trigonometry.azimuth(0.3, 0.8, 1.0) >= 0
      true
  """

  @spec azimuth(number, number, number) :: float
  def azimuth(dec, latitude, hour_angle) do
    az =
      :math.atan2(
        :math.sin(hour_angle),
        :math.cos(hour_angle) * :math.sin(latitude) - :math.tan(dec) * :math.cos(latitude)
      )

    :math.fmod(az + 2 * @pi, 2 * @pi)
  end

  @doc """
  Local hour angle: H = LST - α, normalised to (-π, π].

  ## Examples
      iex> Trigonometry.hour_angle(1.5, 1.5) |> Float.round(4)
      0.0
  """

  @spec hour_angle(number, number) :: float
  def hour_angle(lst, ra) do
    h = :math.fmod(lst - ra, 2 * @pi)
    if h > @pi, do: h - 2 * @pi, else: h
  end

  @doc """
  Parallactic angle: angle between the zenith and north celestial pole as seen from the object.

  tan(q) = sin(H) / (tan(φ) cos(δ) - sin(δ) cos(H))

  ## Examples
      iex> is_float(Trigonometry.parallactic_angle(0.5, 0.3, 0.5))
      true
  """

  @spec parallactic_angle(number, number, number) :: float
  def parallactic_angle(latitude, dec, hour_angle) do
    :math.atan2(
      :math.sin(hour_angle),
      :math.tan(latitude) * :math.cos(dec) - :math.sin(dec) * :math.cos(hour_angle)
    )
  end

  # ---------------------------------------------------------------------------
  # Sidereal Time
  # ---------------------------------------------------------------------------

  @doc """
  Greenwich Mean Sidereal Time (GMST) from Julian Date (IAU formula).

  ## Returns
    GMST in radians [0, 2π)
  """

  @spec gmst(number) :: float
  def gmst(jd) do
    t = (jd - 2_451_545.0) / 36_525.0

    gmst_deg =
      280.46061837 +
        360.98564736629 * (jd - 2_451_545.0) +
        0.000387933 * t * t -
        t * t * t / 38_710_000.0

    :math.fmod(gmst_deg * @deg_to_rad + 4 * @pi, 2 * @pi)
  end

  @doc """
  Local Sidereal Time from GMST and observer longitude.

  ## Parameters
    - jd:        Julian Date
    - longitude: Observer longitude (radians, east positive)

  ## Returns
    LST in radians [0, 2π)

  ## Examples
      iex> Trigonometry.local_sidereal_time(2_451_545.0, 0.0) >= 0
      true
  """

  @spec local_sidereal_time(number, number) :: float
  def local_sidereal_time(jd, longitude) do
    :math.fmod(gmst(jd) + longitude + 2 * @pi, 2 * @pi)
  end

  # ---------------------------------------------------------------------------
  # Sunrise / Sunset
  # ---------------------------------------------------------------------------

  @doc """
  Hour angle of sunrise or sunset.

  cos(H0) = (sin(h0) - sin(φ) sin(δ)) / (cos(φ) cos(δ))
  where h0 is the altitude at rising/setting (default: -0.833° standard refraction).

  ## Returns
    Hour angle H0 (radians), or nil if circumpolar / never rises
  """

  @spec sunrise_hour_angle(number, number, number) :: float | nil
  def sunrise_hour_angle(latitude, dec, h0 \\ -0.01454) do
    cos_h0 =
      (:math.sin(h0) - :math.sin(latitude) * :math.sin(dec)) /
        (:math.cos(latitude) * :math.cos(dec))

    cond do
      cos_h0 < -1.0 -> nil
      cos_h0 > 1.0 -> nil
      true -> :math.acos(cos_h0)
    end
  end

  # ---------------------------------------------------------------------------
  # Atmospheric Refraction
  # ---------------------------------------------------------------------------

  @doc """
  Atmospheric refraction correction (Bennett's formula):
  R ≈ 1.02 / tan(alt + 10.3 / (alt + 5.11))  [arcminutes]

  ## Parameters
    - altitude_deg: Apparent altitude in degrees (above horizon)

  ## Returns
    Refraction in arcminutes

  ## Examples
      iex> Trigonometry.atmospheric_refraction(45.0) |> Float.round(3)
      0.979
  """

  @spec atmospheric_refraction(number) :: float
  def atmospheric_refraction(altitude_deg) when is_number(altitude_deg) do
    1.02 / :math.tan(to_rad(altitude_deg + 10.3 / (altitude_deg + 5.11)))
  end

  # ---------------------------------------------------------------------------
  # Coordinate Conversions
  # ---------------------------------------------------------------------------

  @doc """
  Converts equatorial (RA, Dec) to ecliptic (λ, β) coordinates.

  ## Parameters
    - ra:      Right ascension α (radians)
    - dec:     Declination δ (radians)
    - epsilon: Obliquity of the ecliptic ε (radians, default: J2000 value ≈ 23.44°)

  ## Returns
    {lambda, beta} ecliptic longitude and latitude in radians

  ## Examples
      iex> Trigonometry.equatorial_to_ecliptic(1.5, 0.3) |> elem(0) |> is_float()
      true
  """

  @spec equatorial_to_ecliptic(number, number, number) :: {float, float}
  def equatorial_to_ecliptic(ra, dec, epsilon \\ 0.409_092_600_59) do
    {sin_beta, lambda_raw} = ecliptic_components(ra, dec, epsilon)
    beta = :math.asin(sin_beta)
    {:math.fmod(lambda_raw + 2 * @pi, 2 * @pi), beta}
  end

  defp ecliptic_components(ra, dec, epsilon) do
    cos_eps = :math.cos(epsilon)
    sin_eps = :math.sin(epsilon)
    sin_beta = :math.sin(dec) * cos_eps - :math.cos(dec) * sin_eps * :math.sin(ra)
    lambda = :math.atan2(:math.sin(ra) * cos_eps + :math.tan(dec) * sin_eps, :math.cos(ra))
    {sin_beta, lambda}
  end

  @doc """
  Converts equatorial (RA, Dec) to galactic (l, b) coordinates (J2000).

  ## Returns
    {l, b} galactic longitude and latitude in radians

  ## Examples
      iex> Trigonometry.equatorial_to_galactic(1.5, 0.3) |> elem(0) |> is_float()
      true
  """

  @spec equatorial_to_galactic(number, number) :: {float, float}
  def equatorial_to_galactic(ra, dec) do
    alpha_ngp = @alpha_ngp
    delta_ngp = @delta_ngp
    l_ncp = @l_ncp

    sin_b =
      :math.sin(dec) * :math.sin(delta_ngp) +
        :math.cos(dec) * :math.cos(delta_ngp) * :math.cos(ra - alpha_ngp)

    b = :math.asin(sin_b)

    l =
      :math.fmod(
        l_ncp -
          :math.atan2(
            :math.cos(dec) * :math.sin(ra - alpha_ngp),
            :math.sin(dec) * :math.cos(delta_ngp) -
              :math.cos(dec) * :math.sin(delta_ngp) * :math.cos(ra - alpha_ngp)
          ) +
          2 * @pi,
        2 * @pi
      )

    {l, b}
  end

  @doc """
  Converts galactic (l, b) back to equatorial (RA, Dec) in radians.

  ## Examples
      iex> Trigonometry.galactic_to_equatorial(0.0, 0.0) |> elem(0) |> is_float()
      true
  """

  @spec galactic_to_equatorial(number, number) :: {float, float}
  def galactic_to_equatorial(l, b) do
    alpha_ngp = @alpha_ngp
    delta_ngp = @delta_ngp
    l_ncp = @l_ncp

    sin_dec =
      :math.sin(b) * :math.sin(delta_ngp) +
        :math.cos(b) * :math.cos(delta_ngp) * :math.cos(l_ncp - l)

    dec = :math.asin(sin_dec)

    ra =
      :math.fmod(
        :math.atan2(
          :math.cos(b) * :math.sin(l_ncp - l),
          :math.sin(b) * :math.cos(delta_ngp) -
            :math.cos(b) * :math.sin(delta_ngp) * :math.cos(l_ncp - l)
        ) +
          alpha_ngp +
          2 * @pi,
        2 * @pi
      )

    {ra, dec}
  end

  @doc """
  Approximate solar declination for a given day of year (low-precision):
  δ ≈ -23.45° × cos(360°/365 × (d + 10))

  ## Parameters
    - day_of_year: d (1 = Jan 1)

  ## Returns
    Solar declination in radians

  ## Examples
      iex> Trigonometry.solar_declination(172) |> Float.round(3)
      0.409
  """

  @spec solar_declination(number) :: float

  def solar_declination(day_of_year) do
    to_rad(-23.45 * :math.cos(to_rad(360.0 / 365.0 * (day_of_year + 10))))
  end
end
