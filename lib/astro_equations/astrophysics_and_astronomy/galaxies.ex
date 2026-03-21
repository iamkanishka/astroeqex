defmodule AstroEquations.AstrophysicsAndAstronomy.Galaxies do
  @moduledoc """
  A collection of functions for calculating galaxy properties including:
  - Hubble elliptical galaxy classification
  - Sérsic light profiles
  - Exponential disk density distributions
  - Tully-Fisher and Faber-Jackson relations
  - Galaxy rotation curves (flat and Keplerian)
  - Dark matter halo profiles (NFW, enclosed mass)
  - Star formation rates (UV, Hα, IR calibrations)
  - Luminosity functions (Schechter)
  - Virial and dynamical mass estimates
  - Concentration parameter and halo mass
  - Galaxy merger timescale
  - Mass-to-light ratio

  All calculations use consistent astronomical units unless stated otherwise.
  """

  # ---------------------------------------------------------------------------
  # Classification
  # ---------------------------------------------------------------------------

  @doc """
  Classifies an elliptical galaxy based on its apparent flattening (Hubble scheme).

  Type = round(10 × (a − b) / a), clamped to E0–E7.

  ## Parameters
    - a: Semi-major axis length
    - b: Semi-minor axis length

  ## Returns
    Hubble classification string (E0–E7)

  ## Examples
      iex> Galaxies.hubble_classify(10, 10)
      "E0"
      iex> Galaxies.hubble_classify(10, 5)
      "E5"
  """
  @spec hubble_classify(number, number) :: String.t()
  def hubble_classify(a, b) do
    type = round((a - b) / a * 10)
    "E#{min(max(type, 0), 7)}"
  end

  # ---------------------------------------------------------------------------
  # Surface Brightness Profiles
  # ---------------------------------------------------------------------------

  @doc """
  Calculates surface brightness using the Sérsic profile.

  I(r) = I₀ × exp(−b_n × ((r/r_e)^(1/n) − 1))

  where b_n ≈ 2n − 1/3 + 0.009876/n (Capaccioli 1989 approximation, valid for n > 0.36).
  Special cases: n = 1 (exponential disk), n = 4 (de Vaucouleurs).

  ## Parameters
    - i0:  Central surface brightness
    - r:   Radius from center
    - r_e: Effective (half-light) radius
    - n:   Sérsic index

  ## Returns
    Surface brightness at radius r

  ## Examples
      iex> Galaxies.sersic_profile(100, 1, 1, 1) |> Float.round(4)
      100.0
      iex> Galaxies.sersic_profile(100, 2, 1, 4) |> Float.round(4)
      23.4412
  """
  @spec sersic_profile(number, number, number, number) :: float
  def sersic_profile(i0, r, r_e, n) do
    bn = 2 * n - 1 / 3 + 0.009876 / n
    i0 * :math.exp(-bn * (:math.pow(r / r_e, 1 / n) - 1))
  end

  @doc """
  Calculates surface brightness using the de Vaucouleurs (r^1/4) profile.

  Equivalent to the Sérsic profile with n = 4 and b₄ = 7.669.

  I(r) = I_e × exp(−7.669 × ((r/r_e)^(1/4) − 1))

  ## Parameters
    - i_e: Surface brightness at the effective radius
    - r:   Radius
    - r_e: Effective radius

  ## Returns
    Surface brightness at radius r

  ## Examples
      iex> Galaxies.de_vaucouleurs_profile(100, 1, 1) |> Float.round(4)
      100.0
  """
  @spec de_vaucouleurs_profile(number, number, number) :: float
  def de_vaucouleurs_profile(i_e, r, r_e) do
    i_e * :math.exp(-7.669 * (:math.pow(r / r_e, 0.25) - 1))
  end

  @doc """
  Mean surface brightness within the effective radius for a Sérsic profile.

  ⟨I⟩_e = I_e × e^b_n × Γ(2n+1) / (2 n b_n^(2n))  (approximation)

  For a quick estimate this uses the ratio I₀/2 (midpoint approximation).

  ## Parameters
    - i0:  Central surface brightness
    - n:   Sérsic index

  ## Returns
    Approximate mean surface brightness within r_e

  ## Examples
      iex> Galaxies.mean_surface_brightness_within_re(100.0, 1.0) > 0
      true
  """
  @spec mean_surface_brightness_within_re(number, number) :: float
  def mean_surface_brightness_within_re(i0, n) do
    bn = 2 * n - 1 / 3 + 0.009876 / n
    i0 * :math.exp(-bn) * :math.exp(bn) / (2 * n)
  end

  # ---------------------------------------------------------------------------
  # Density Profiles
  # ---------------------------------------------------------------------------

  @doc """
  Stellar density in a disk galaxy (double-exponential profile).

  ρ(R, z) = ρ₀ × exp(−|z|/z₀) × exp(−R/h)

  ## Parameters
    - p0: Central density
    - r:  Cylindrical radius from galactic center
    - z:  Height above the disk plane
    - h:  Radial scale length
    - z0: Vertical scale height

  ## Returns
    Stellar density at (R, z)

  ## Examples
      iex> Galaxies.disk_density(1, 0, 0, 1, 1)
      1.0
      iex> Galaxies.disk_density(1, 1, 1, 1, 1) |> Float.round(4)
      0.1353
  """
  @spec disk_density(number, number, number, number, number) :: float
  def disk_density(p0, r, z, h, z0) do
    p0 * :math.exp(-abs(z) / z0) * :math.exp(-r / h)
  end

  @doc """
  NFW (Navarro-Frenk-White) dark matter halo density profile.

  ρ(r) = ρ_s / ((r/r_s) × (1 + r/r_s)²)

  ## Parameters
    - rho_s: Characteristic density
    - r:     Radius
    - r_s:   Scale radius

  ## Returns
    Dark matter density at radius r

  ## Examples
      iex> Galaxies.nfw_profile(1.0e7, 3, 3) > 0
      true
  """
  @spec nfw_profile(number, number, number) :: float
  def nfw_profile(rho_s, r, r_s) do
    x = r / r_s
    rho_s / (x * :math.pow(1 + x, 2))
  end

  @doc """
  Mass enclosed within radius r for an NFW profile.

  M(r) = 4π ρ_s r_s³ × [ln(1 + r/r_s) − r/r_s / (1 + r/r_s)]

  ## Parameters
    - rho_s: Characteristic density
    - r:     Radius
    - r_s:   Scale radius

  ## Returns
    Enclosed mass in the same units as rho_s × r_s³

  ## Examples
      iex> Galaxies.nfw_enclosed_mass(1.0e7, 10, 3) > 0
      true
  """
  @spec nfw_enclosed_mass(number, number, number) :: float
  def nfw_enclosed_mass(rho_s, r, r_s) do
    x = r / r_s
    4 * :math.pi() * rho_s * r_s ** 3 * (:math.log(1 + x) - x / (1 + x))
  end

  @doc """
  NFW concentration parameter: c = r₂₀₀ / r_s

  where r₂₀₀ is the radius enclosing 200× the critical density.

  ## Parameters
    - r200: Virial radius r₂₀₀
    - r_s:  NFW scale radius

  ## Returns
    Concentration parameter c

  ## Examples
      iex> Galaxies.concentration_parameter(200.0, 20.0) |> Float.round(4)
      10.0
  """
  @spec concentration_parameter(number, number) :: float
  def concentration_parameter(r200, r_s), do: r200 / r_s

  # ---------------------------------------------------------------------------
  # Rotation Curves
  # ---------------------------------------------------------------------------

  @doc """
  Keplerian (point-mass) rotation velocity at radius r: v = √(G M / r)

  ## Parameters
    - mass:   Total mass interior to radius r (kg)
    - radius: Radius (m)

  ## Returns
    Orbital velocity in m/s

  ## Examples
      iex> Galaxies.keplerian_rotation(1.989e30, 1.0e11) > 0
      true
  """
  @spec keplerian_rotation(number, number) :: float
  def keplerian_rotation(mass, radius) do
    :math.sqrt(6.674e-11 * mass / radius)
  end

  @doc """
  Flat rotation curve — returns the constant velocity for the flat part.

  Observed in spiral galaxies due to the dark matter halo; independent of radius.

  ## Parameters
    - v_flat: The flat (asymptotic) rotation velocity

  ## Returns
    v_flat unchanged (identity; useful for modelling)

  ## Examples
      iex> Galaxies.flat_rotation_curve(220.0)
      220.0
  """
  @spec flat_rotation_curve(number) :: float
  def flat_rotation_curve(v_flat), do: v_flat

  @doc """
  Galaxy dynamical mass from flat rotation speed and radius: M = v² r / G

  ## Parameters
    - v_rot:  Rotation velocity (m/s)
    - radius: Galactocentric radius (m)

  ## Returns
    Enclosed dynamical mass in kg

  ## Examples
      iex> Galaxies.dynamical_mass_from_rotation(2.2e5, 1.0e21) > 0
      true
  """
  @spec dynamical_mass_from_rotation(number, number) :: float
  def dynamical_mass_from_rotation(v_rot, radius) do
    v_rot * v_rot * radius / 6.674e-11
  end

  # ---------------------------------------------------------------------------
  # Tully-Fisher & Faber-Jackson
  # ---------------------------------------------------------------------------

  @doc """
  Estimates luminosity from rotation velocity using the Tully-Fisher relation.

  L = A × v_rot^α  (baryonic Tully-Fisher: L ∝ v^4)

  ## Parameters
    - v_rot: Maximum rotation velocity in km/s
    - a:     Normalisation (default: 50.0 for L in L☉, v in km/s)
    - alpha: Power-law exponent (default: 4.0)

  ## Returns
    Estimated luminosity in solar luminosities

  ## Examples
      iex> Galaxies.tully_fisher(100) > 0
      true
  """
  @spec tully_fisher(number, number, number) :: float
  def tully_fisher(v_rot, a \\ 50.0, alpha \\ 4.0), do: a * :math.pow(v_rot, alpha)

  @doc """
  Estimates luminosity from velocity dispersion using the Faber-Jackson relation.

  L ∝ σ^4  (elliptical galaxies)

  ## Parameters
    - sigma: Stellar velocity dispersion in km/s
    - a:     Normalisation (default: 1.0)
    - alpha: Exponent (default: 4.0)

  ## Returns
    Estimated luminosity in solar luminosities

  ## Examples
      iex> Galaxies.faber_jackson(200) > 0
      true
  """
  @spec faber_jackson(number, number, number) :: float
  def faber_jackson(sigma, a \\ 1.0, alpha \\ 4.0), do: a * :math.pow(sigma, alpha)

  # ---------------------------------------------------------------------------
  # Star Formation
  # ---------------------------------------------------------------------------

  @doc """
  Specific star formation rate: sSFR = SFR / M_stellar

  ## Parameters
    - sfr:    Star formation rate in M☉/yr
    - m_star: Stellar mass in M☉

  ## Returns
    sSFR in yr⁻¹

  ## Examples
      iex> Galaxies.specific_sfr(10.0, 1.0e11)
      1.0e-10
  """
  @spec specific_sfr(number, number) :: float
  def specific_sfr(sfr, m_star), do: sfr / m_star

  @doc """
  Estimates SFR from UV luminosity (Kennicutt 1998):
  SFR [M☉/yr] = L_UV [erg/s] / 8.0×10²⁷

  ## Examples
      iex> Galaxies.sfr_from_uv(8.0e28) > 0
      true
  """
  @spec sfr_from_uv(number) :: float
  def sfr_from_uv(l_uv), do: l_uv / 8.0e27

  @doc """
  Estimates SFR from Hα luminosity (Kennicutt 1998):
  SFR [M☉/yr] = L_Hα [erg/s] / 1.26×10⁴¹

  ## Examples
      iex> Galaxies.sfr_from_halpha(1.26e41) |> Float.round(2)
      1.0
  """
  @spec sfr_from_halpha(number) :: float
  def sfr_from_halpha(l_ha), do: l_ha / 1.26e41

  @doc """
  Estimates SFR from total infrared luminosity (Kennicutt 1998):
  SFR [M☉/yr] = L_IR [erg/s] / 5.8×10⁴³

  Infrared traces star formation reprocessed by dust and is less sensitive
  to extinction than UV or optical tracers.

  ## Examples
      iex> Galaxies.sfr_from_ir(5.8e43) |> Float.round(2)
      1.0
  """
  @spec sfr_from_ir(number) :: float
  def sfr_from_ir(l_ir), do: l_ir / 5.8e43

  @doc """
  Kennicutt-Schmidt law: SFR surface density ∝ gas surface density^N.

  Σ_SFR = A × Σ_gas^n  (n ≈ 1.4 for the Kennicutt 1998 fit)

  ## Parameters
    - sigma_gas: Gas surface density in M☉/pc²
    - a:         Normalisation (default: 2.5e-4 M☉/yr/pc²)
    - n:         Power-law index (default: 1.4)

  ## Returns
    Star formation rate surface density in M☉/yr/pc²

  ## Examples
      iex> Galaxies.kennicutt_schmidt(10.0) > 0
      true
  """
  @spec kennicutt_schmidt(number, number, number) :: float
  def kennicutt_schmidt(sigma_gas, a \\ 2.5e-4, n \\ 1.4) do
    a * :math.pow(sigma_gas, n)
  end

  # ---------------------------------------------------------------------------
  # Luminosity Function
  # ---------------------------------------------------------------------------

  @doc """
  Evaluates the Schechter luminosity function φ(L).

  φ(L) = (φ*/L*) × (L/L*)^α × exp(−L/L*)

  ## Parameters
    - phi_star: Normalisation density (Mpc⁻³)
    - l:        Luminosity
    - l_star:   Characteristic luminosity
    - alpha:    Faint-end slope

  ## Returns
    Number density per unit luminosity

  ## Examples
      iex> Galaxies.schechter_function(0.016, 1.0e10, 1.0e10, -1.07) > 0
      true
  """
  @spec schechter_function(number, number, number, number) :: float
  def schechter_function(phi_star, l, l_star, alpha) do
    x = l / l_star
    phi_star / l_star * :math.pow(x, alpha) * :math.exp(-x)
  end

  # ---------------------------------------------------------------------------
  # Virial Theorem & Dynamics
  # ---------------------------------------------------------------------------

  @doc """
  Estimates dynamical mass from the virial theorem.

  M_vir = 5 σ² R / G

  ## Parameters
    - sigma:  Line-of-sight velocity dispersion in m/s
    - radius: Half-light radius in meters

  ## Returns
    Virial mass in kg

  ## Examples
      iex> Galaxies.virial_mass(2.0e5, 1.0e21) > 0
      true
  """
  @spec virial_mass(number, number) :: float
  def virial_mass(sigma, radius) do
    5 * :math.pow(sigma, 2) * radius / 6.674e-11
  end

  @doc """
  Approximate galaxy merger timescale using dynamical friction (Lacey & Cole 1993).

  τ_merge ≈ f_dyn × (M_host / M_satellite) × t_dyn

  where t_dyn = r / σ is the dynamical time and f_dyn ≈ 0.3–1.0.

  ## Parameters
    - m_host:     Host halo mass (kg or M☉ — consistent units)
    - m_satellite: Satellite mass (same units)
    - r_orbit:    Orbital radius (m)
    - sigma_host: Velocity dispersion of host (m/s)
    - f_dyn:      Dynamical friction factor (default: 0.5)

  ## Returns
    Merger timescale in seconds

  ## Examples
      iex> Galaxies.galaxy_merger_timescale(1.0e12, 1.0e10, 1.0e22, 2.0e5) > 0
      true
  """
  @spec galaxy_merger_timescale(number, number, number, number, number) :: float
  def galaxy_merger_timescale(m_host, m_satellite, r_orbit, sigma_host, f_dyn \\ 0.5) do
    t_dyn = r_orbit / sigma_host
    f_dyn * (m_host / m_satellite) * t_dyn
  end

  @doc """
  Mass-to-light ratio: Υ = M / L

  ## Parameters
    - mass:       Total mass (M☉ or kg)
    - luminosity: Total luminosity (L☉ or W — same system)

  ## Returns
    Mass-to-light ratio Υ in M☉/L☉ (or corresponding units)

  ## Examples
      iex> Galaxies.mass_to_light_ratio(1.0e11, 2.0e10) |> Float.round(4)
      5.0
  """
  @spec mass_to_light_ratio(number, number) :: float
  def mass_to_light_ratio(mass, luminosity), do: mass / luminosity
end
