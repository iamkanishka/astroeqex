defmodule AstroEquations.Mathematics.Geometry do
  @moduledoc """
  Geometric utilities relevant to astronomy and astrophysics.

  Provides:
  - Spherical coordinate conversions (equatorial ↔ galactic ↔ ecliptic)
  - Great-circle angular separation
  - Solid angle calculations
  - Coordinate rotation matrices
  - Precession of equinoxes
  - Parallax geometry
  - Orbital geometry helpers (semi-latus rectum, eccentric/true anomaly)
  - Gravitational lensing geometry (Einstein radius)
  """

  @pi :math.pi()
  @deg_to_rad @pi / 180.0
  @rad_to_deg 180.0 / @pi
  # arcsec per radian
  @arcsec_per_rad 206_265.0

  # ---------------------------------------------------------------------------
  # Unit Conversion Helpers
  # ---------------------------------------------------------------------------

  @doc """
  Converts degrees to radians.

  ## Examples
      iex> Geometry.deg_to_rad(180) |> Float.round(4)
      3.1416
  """
  @spec deg_to_rad(number) :: float
  def deg_to_rad(deg), do: deg * @deg_to_rad

  @doc """
  Converts radians to degrees.

  ## Examples
      iex> Geometry.rad_to_deg(:math.pi()) |> Float.round(4)
      180.0
  """
  @spec rad_to_deg(number) :: float
  def rad_to_deg(rad), do: rad * @rad_to_deg

  @doc """
  Converts hours of right ascension to degrees.

  ## Examples
      iex> Geometry.hours_to_deg(12)
      180.0
  """
  @spec hours_to_deg(number) :: float
  def hours_to_deg(h), do: h * 15.0

  @doc """
  Converts degrees to hours of right ascension.

  ## Examples
      iex> Geometry.deg_to_hours(180)
      12.0
  """
  @spec deg_to_hours(number) :: float
  def deg_to_hours(deg), do: deg / 15.0

  @doc """
  Converts arcseconds to radians.

  ## Examples
      iex> Geometry.arcsec_to_rad(206_265.0) |> Float.round(4)
      1.0
  """
  @spec arcsec_to_rad(number) :: float
  def arcsec_to_rad(arcsec), do: arcsec / @arcsec_per_rad

  @doc """
  Converts radians to arcseconds.

  ## Examples
      iex> Geometry.rad_to_arcsec(1.0) |> Float.round(0)
      206265.0
  """
  @spec rad_to_arcsec(number) :: float
  def rad_to_arcsec(rad), do: rad * @arcsec_per_rad

  # ---------------------------------------------------------------------------
  # Angular Separation
  # ---------------------------------------------------------------------------

  @doc """
  Great-circle angular separation (Haversine formula) in radians.

  ## Parameters
    - ra1, dec1: RA and Dec of first point (radians)
    - ra2, dec2: RA and Dec of second point (radians)

  ## Returns
    Angular separation in radians

  ## Examples
      iex> Geometry.angular_separation(0, 0, :math.pi/2, 0) |> Float.round(6)
      1.570796
  """
  @spec angular_separation(number, number, number, number) :: float
  def angular_separation(ra1, dec1, ra2, dec2) do
    d_ra = (ra2 - ra1) / 2
    d_dec = (dec2 - dec1) / 2

    a =
      :math.pow(:math.sin(d_dec), 2) +
        :math.cos(dec1) * :math.cos(dec2) * :math.pow(:math.sin(d_ra), 2)

    2 * :math.asin(:math.sqrt(a))
  end

  @doc """
  Great-circle angular separation in degrees (inputs in degrees).

  ## Parameters
    - ra1_deg, dec1_deg: First position in degrees
    - ra2_deg, dec2_deg: Second position in degrees

  ## Returns
    Angular separation in degrees

  ## Examples
      iex> Geometry.angular_separation_deg(0.0, 0.0, 90.0, 0.0) |> Float.round(4)
      90.0
  """
  @spec angular_separation_deg(number, number, number, number) :: float
  def angular_separation_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg) do
    angular_separation(
      deg_to_rad(ra1_deg),
      deg_to_rad(dec1_deg),
      deg_to_rad(ra2_deg),
      deg_to_rad(dec2_deg)
    )
    |> rad_to_deg()
  end

  @doc """
  Position angle from point 1 toward point 2, measured east of north.

  ## Parameters
    - ra1, dec1: Reference point (radians)
    - ra2, dec2: Target point (radians)

  ## Returns
    Position angle in radians [0, 2π)

  ## Examples
      iex> Geometry.position_angle(0.0, 0.0, 0.1, 0.0) >= 0
      true
  """
  @spec position_angle(number, number, number, number) :: float
  def position_angle(ra1, dec1, ra2, dec2) do
    d_ra = ra2 - ra1

    pa =
      :math.atan2(
        :math.sin(d_ra),
        :math.cos(dec1) * :math.tan(dec2) - :math.sin(dec1) * :math.cos(d_ra)
      )

    :math.fmod(pa + 2 * @pi, 2 * @pi)
  end

  # ---------------------------------------------------------------------------
  # Solid Angle
  # ---------------------------------------------------------------------------

  @doc """
  Solid angle of a cone with half-angle θ: Ω = 2π(1 - cos θ).

  ## Parameters
    - theta: Half-angle in radians

  ## Returns
    Solid angle in steradians

  ## Examples
      iex> Geometry.solid_angle_cone(0) |> Float.round(4)
      0.0
  """
  @spec solid_angle_cone(number) :: float
  def solid_angle_cone(theta), do: 2 * @pi * (1 - :math.cos(theta))

  @doc """
  Solid angle of a small sky rectangle near declination δ: Ω ≈ Δα Δδ cos δ.

  ## Parameters
    - d_alpha: Width in RA (radians)
    - d_delta: Height in Dec (radians)
    - delta:   Central declination (radians)

  ## Returns
    Solid angle in steradians

  ## Examples
      iex> Geometry.solid_angle_rectangle(0.01, 0.01, 0) |> Float.round(6)
      0.0001
  """
  @spec solid_angle_rectangle(number, number, number) :: float
  def solid_angle_rectangle(d_alpha, d_delta, delta), do: d_alpha * d_delta * :math.cos(delta)

  @doc """
  Total solid angle of a sphere: 4π steradians.

  ## Examples
      iex> Geometry.sphere_solid_angle() |> Float.round(4)
      12.5664
  """
  @spec sphere_solid_angle() :: float
  def sphere_solid_angle(), do: 4 * @pi

  # ---------------------------------------------------------------------------
  # Orbital Geometry
  # ---------------------------------------------------------------------------

  @doc """
  Semi-latus rectum: l = a(1 - e²).

  ## Parameters
    - a: Semi-major axis
    - e: Eccentricity

  ## Examples
      iex> Geometry.semi_latus_rectum(1.0, 0.0)
      1.0
  """
  @spec semi_latus_rectum(number, number) :: float
  def semi_latus_rectum(a, e), do: a * (1 - e * e)

  @doc """
  Orbital radius at true anomaly ν: r = l/(1 + e cos ν).

  ## Parameters
    - a:   Semi-major axis
    - e:   Eccentricity
    - nu:  True anomaly (radians)

  ## Returns
    Radial distance r

  ## Examples
      iex> Geometry.orbital_radius(1.0, 0.5, 0) |> Float.round(4)
      0.5
  """
  @spec orbital_radius(number, number, number) :: float
  def orbital_radius(a, e, nu) do
    l = semi_latus_rectum(a, e)
    l / (1 + e * :math.cos(nu))
  end

  @doc """
  Periapsis distance: r_p = a(1 - e).

  ## Examples
      iex> Geometry.periapsis(1.0, 0.5)
      0.5
  """
  @spec periapsis(number, number) :: float
  def periapsis(a, e), do: a * (1 - e)

  @doc """
  Apoapsis distance: r_a = a(1 + e).

  ## Examples
      iex> Geometry.apoapsis(1.0, 0.5)
      1.5
  """
  @spec apoapsis(number, number) :: float
  def apoapsis(a, e), do: a * (1 + e)

  @doc """
  Orbital eccentricity from periapsis and apoapsis distances: e = (r_a - r_p)/(r_a + r_p).

  ## Examples
      iex> Geometry.eccentricity_from_apsides(0.5, 1.5) |> Float.round(4)
      0.5
  """
  @spec eccentricity_from_apsides(number, number) :: float
  def eccentricity_from_apsides(r_p, r_a), do: (r_a - r_p) / (r_a + r_p)

  @doc """
  Eccentric anomaly from mean anomaly (solves Kepler's equation via Newton-Raphson).

  Solves M = E - e sin(E)

  ## Parameters
    - m:   Mean anomaly (radians)
    - e:   Eccentricity
    - tol: Convergence tolerance (default 1.0e-10)

  ## Returns
    Eccentric anomaly in radians

  ## Examples
      iex> Geometry.eccentric_anomaly(0.0, 0.5) |> Float.round(4)
      0.0
  """
  @spec eccentric_anomaly(number, number, number) :: float
  def eccentric_anomaly(m, e, tol \\ 1.0e-10) do
    do_kepler(m, e, m, tol)
  end

  defp do_kepler(m, e, big_e, tol) do
    delta = (m - big_e + e * :math.sin(big_e)) / (1 - e * :math.cos(big_e))
    if abs(delta) < tol, do: big_e + delta, else: do_kepler(m, e, big_e + delta, tol)
  end

  @doc """
  True anomaly from eccentric anomaly: tan(ν/2) = √((1+e)/(1-e)) tan(E/2).

  ## Parameters
    - big_e: Eccentric anomaly (radians)
    - e:     Eccentricity

  ## Returns
    True anomaly in radians [0, 2π)

  ## Examples
      iex> Geometry.true_anomaly(0.0, 0.5) |> Float.round(4)
      0.0
  """
  @spec true_anomaly(number, number) :: float
  def true_anomaly(big_e, e) do
    half_nu =
      :math.atan2(
        :math.sqrt(1 + e) * :math.sin(big_e / 2),
        :math.sqrt(1 - e) * :math.cos(big_e / 2)
      )

    :math.fmod(2 * half_nu + 2 * @pi, 2 * @pi)
  end

  # ---------------------------------------------------------------------------
  # Gravitational Lensing
  # ---------------------------------------------------------------------------

  @doc """
  Einstein radius for gravitational lensing: θ_E = √(4GM D_LS / (c² D_L D_S))

  ## Parameters
    - mass:  Lens mass (kg)
    - d_l:   Distance from observer to lens (m)
    - d_s:   Distance from observer to source (m)
    - d_ls:  Distance from lens to source (m)

  ## Returns
    Einstein radius in radians

  ## Examples
      iex> Geometry.einstein_radius(1.989e30, 1.0e22, 2.0e22, 1.0e22) > 0
      true
  """
  @spec einstein_radius(number, number, number, number) :: float
  def einstein_radius(mass, d_l, d_s, d_ls) do
    g = 6.674e-11
    c = 2.99792458e8
    :math.sqrt(4 * g * mass * d_ls / (:math.pow(c, 2) * d_l * d_s))
  end

  @doc """
  Microlensing magnification from impact parameter: A(u) = (u² + 2) / (u √(u² + 4))

  ## Parameters
    - u: Impact parameter in units of the Einstein radius

  ## Returns
    Magnification factor A

  ## Examples
      iex> Geometry.microlensing_magnification(1.0) |> Float.round(3)
      1.342
  """
  @spec microlensing_magnification(number) :: float
  def microlensing_magnification(u) do
    u2 = u * u
    (u2 + 2) / (u * :math.sqrt(u2 + 4))
  end

  @doc """
  Parallax to distance conversion: d = 1 / p  (where p is in arcseconds, d in parsecs).

  ## Parameters
    - parallax_arcsec: Parallax angle (arcseconds)

  ## Returns
    Distance in parsecs

  ## Examples
      iex> Geometry.parallax_to_parsecs(0.1) |> Float.round(1)
      10.0
  """
  @spec parallax_to_parsecs(number) :: float
  def parallax_to_parsecs(parallax_arcsec), do: 1.0 / parallax_arcsec

  @doc """
  Proper motion magnitude from RA and Dec components.

  μ = √(μ_α*² + μ_δ²)  where μ_α* = μ_α cos δ

  ## Parameters
    - mu_alpha_star: Proper motion in RA × cos(Dec) (mas/yr or arcsec/yr)
    - mu_delta:      Proper motion in Dec (same units)

  ## Returns
    Total proper motion magnitude

  ## Examples
      iex> Geometry.proper_motion_magnitude(3.0, 4.0) |> Float.round(1)
      5.0
  """
  @spec proper_motion_magnitude(number, number) :: float
  def proper_motion_magnitude(mu_alpha_star, mu_delta) do
    :math.sqrt(mu_alpha_star * mu_alpha_star + mu_delta * mu_delta)
  end

  # ---------------------------------------------------------------------------
  # Coordinate Conversions
  # ---------------------------------------------------------------------------

  @doc """
  Converts equatorial (RA, Dec) to a unit Cartesian vector.

  ## Parameters
    - ra:  Right ascension (radians)
    - dec: Declination (radians)

  ## Returns
    {x, y, z} unit vector on the celestial sphere

  ## Examples
      iex> elem(Geometry.equatorial_to_cartesian(0.0, 0.0), 0) |> Float.round(4)
      1.0
  """
  @spec equatorial_to_cartesian(number, number) :: {float, float, float}
  def equatorial_to_cartesian(ra, dec) do
    {
      :math.cos(dec) * :math.cos(ra),
      :math.cos(dec) * :math.sin(ra),
      :math.sin(dec)
    }
  end

  @doc """
  Converts a unit Cartesian vector to equatorial (RA, Dec) in radians.

  ## Parameters
    - x, y, z: Cartesian coordinates

  ## Returns
    {ra, dec} in radians with ra in [0, 2π)

  ## Examples
      iex> elem(Geometry.cartesian_to_equatorial(1.0, 0.0, 0.0), 0) |> Float.round(4)
      0.0
  """
  @spec cartesian_to_equatorial(number, number, number) :: {float, float}
  def cartesian_to_equatorial(x, y, z) do
    dec = :math.asin(z)
    ra = :math.fmod(:math.atan2(y, x) + 2 * @pi, 2 * @pi)
    {ra, dec}
  end
end
