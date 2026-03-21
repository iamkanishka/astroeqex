defmodule AstroEquations.Mathematics.Notation do
  @moduledoc """
  Astronomical notation, unit conversions, and physical constants used
  throughout the AstroEquations package.

  Provides:
  - Physical and astronomical constants (SI)
  - Unit conversion functions (length, mass, time, angle, energy, temperature)
  - Sexagesimal (hh:mm:ss / dd:mm:ss) ↔ decimal conversions
  - Julian Date conversions
  - Scientific notation helpers
  """

  # ---------------------------------------------------------------------------
  # Physical Constants (SI)
  # ---------------------------------------------------------------------------

  @doc "Speed of light in vacuum: c = 2.99792458×10⁸ m/s (exact)."
  @spec speed_of_light() :: float
  def speed_of_light, do: 2.99792458e8

  @doc "Newtonian gravitational constant: G = 6.6743×10⁻¹¹ m³ kg⁻¹ s⁻²."
  @spec gravitational_constant() :: float
  def gravitational_constant, do: 6.67430e-11

  @doc "Planck constant: h = 6.62607015×10⁻³⁴ J·s (exact)."
  @spec planck_constant() :: float
  def planck_constant, do: 6.62607015e-34

  @doc "Reduced Planck constant: ħ = h/(2π) ≈ 1.054571817×10⁻³⁴ J·s."
  @spec hbar() :: float
  def hbar, do: 1.054571817e-34

  @doc "Boltzmann constant: k_B = 1.380649×10⁻²³ J/K (exact)."
  @spec boltzmann_constant() :: float
  def boltzmann_constant, do: 1.380649e-23

  @doc "Stefan-Boltzmann constant: σ = 5.670374×10⁻⁸ W m⁻² K⁻⁴."
  @spec stefan_boltzmann() :: float
  def stefan_boltzmann, do: 5.670374419e-8

  @doc "Radiation constant: a = 4σ/c = 7.566×10⁻¹⁶ J m⁻³ K⁻⁴."
  @spec radiation_constant() :: float
  def radiation_constant, do: 7.565723e-16

  @doc "Proton mass: m_p = 1.67262×10⁻²⁷ kg."
  @spec proton_mass() :: float
  def proton_mass, do: 1.67262192369e-27

  @doc "Electron mass: m_e = 9.10938×10⁻³¹ kg."
  @spec electron_mass() :: float
  def electron_mass, do: 9.1093837015e-31

  @doc "Elementary charge: e = 1.602176634×10⁻¹⁹ C (exact)."
  @spec elementary_charge() :: float
  def elementary_charge, do: 1.602176634e-19

  @doc "Thomson cross-section: σ_T = 6.6524×10⁻²⁹ m²."
  @spec thomson_cross_section() :: float
  def thomson_cross_section, do: 6.6524587158e-29

  @doc "Wien displacement constant: b = 2.897771955×10⁻³ m·K."
  @spec wien_constant() :: float
  def wien_constant, do: 2.897771955e-3

  @doc "Avogadro constant: N_A = 6.02214076×10²³ mol⁻¹ (exact)."
  @spec avogadro() :: float
  def avogadro, do: 6.02214076e23

  # ---------------------------------------------------------------------------
  # Astronomical Constants
  # ---------------------------------------------------------------------------

  @doc "Solar mass: M☉ = 1.98892×10³⁰ kg."
  @spec solar_mass() :: float
  def solar_mass, do: 1.98892e30

  @doc "Solar radius: R☉ = 6.957×10⁸ m."
  @spec solar_radius() :: float
  def solar_radius, do: 6.957e8

  @doc "Solar luminosity: L☉ = 3.828×10²⁶ W."
  @spec solar_luminosity() :: float
  def solar_luminosity, do: 3.828e26

  @doc "Solar effective surface temperature: T☉ = 5778 K."
  @spec solar_temperature() :: float
  def solar_temperature, do: 5778.0

  @doc "Astronomical unit: AU = 1.495978707×10¹¹ m (exact)."
  @spec astronomical_unit() :: float
  def astronomical_unit, do: 1.495978707e11

  @doc "Parsec: 1 pc = 3.085677581×10¹⁶ m."
  @spec parsec() :: float
  def parsec, do: 3.085677581e16

  @doc "Light-year: 1 ly = 9.460730473×10¹⁵ m."
  @spec light_year() :: float
  def light_year, do: 9.4607304725808e15

  @doc "Earth mass: M⊕ = 5.9722×10²⁴ kg."
  @spec earth_mass() :: float
  def earth_mass, do: 5.9722e24

  @doc "Earth equatorial radius: R⊕ = 6.371×10⁶ m."
  @spec earth_radius() :: float
  def earth_radius, do: 6.371e6

  @doc "Jupiter mass: M_J = 1.8982×10²⁷ kg."
  @spec jupiter_mass() :: float
  def jupiter_mass, do: 1.8982e27

  @doc "Hubble constant: H₀ ≈ 67.36 km/s/Mpc (Planck 2018)."
  @spec hubble_constant() :: float
  def hubble_constant, do: 67.36

  # ---------------------------------------------------------------------------
  # Length Conversions
  # ---------------------------------------------------------------------------

  @doc "Converts parsecs to metres."
  @spec pc_to_m(number) :: float
  def pc_to_m(pc), do: pc * parsec()

  @doc "Converts metres to parsecs."
  @spec m_to_pc(number) :: float
  def m_to_pc(m), do: m / parsec()

  @doc "Converts kiloparsecs to parsecs: 1 kpc = 1000 pc."
  @spec kpc_to_pc(number) :: float
  def kpc_to_pc(kpc), do: kpc * 1_000.0

  @doc "Converts megaparsecs to parsecs: 1 Mpc = 10⁶ pc."
  @spec mpc_to_pc(number) :: float
  def mpc_to_pc(mpc), do: mpc * 1_000_000.0

  @doc "Converts astronomical units to metres."
  @spec au_to_m(number) :: float
  def au_to_m(au), do: au * astronomical_unit()

  @doc "Converts light-years to metres."
  @spec ly_to_m(number) :: float
  def ly_to_m(ly), do: ly * light_year()

  @doc "Converts light-years to parsecs: 1 ly ≈ 0.3066 pc."
  @spec ly_to_pc(number) :: float
  def ly_to_pc(ly), do: ly_to_m(ly) / parsec()

  @doc "Converts parsecs to light-years."
  @spec pc_to_ly(number) :: float
  def pc_to_ly(pc), do: pc_to_m(pc) / light_year()

  # ---------------------------------------------------------------------------
  # Angular Conversions
  # ---------------------------------------------------------------------------

  @doc """
  Converts degrees to radians.

  ## Examples
      iex> Notation.deg_to_rad(180) |> Float.round(4)
      3.1416
  """
  @spec deg_to_rad(number) :: float
  def deg_to_rad(deg), do: deg * :math.pi() / 180.0

  @doc """
  Converts radians to degrees.

  ## Examples
      iex> Notation.rad_to_deg(:math.pi()) |> Float.round(4)
      180.0
  """
  @spec rad_to_deg(number) :: float
  def rad_to_deg(rad), do: rad * 180.0 / :math.pi()

  @doc """
  Converts arcseconds to radians.

  ## Examples
      iex> Notation.arcsec_to_rad(206_265.0) |> Float.round(4)
      1.0
  """
  @spec arcsec_to_rad(number) :: float
  def arcsec_to_rad(arcsec), do: arcsec / 206_265.0

  @doc """
  Converts radians to arcseconds.

  ## Examples
      iex> Notation.rad_to_arcsec(1.0) |> Float.round(0)
      206265.0
  """
  @spec rad_to_arcsec(number) :: float
  def rad_to_arcsec(rad), do: rad * 206_265.0

  @doc "Converts arcminutes to degrees: deg = arcmin / 60."
  @spec arcmin_to_deg(number) :: float
  def arcmin_to_deg(arcmin), do: arcmin / 60.0

  @doc "Converts arcseconds to degrees: deg = arcsec / 3600."
  @spec arcsec_to_deg(number) :: float
  def arcsec_to_deg(arcsec), do: arcsec / 3600.0

  # ---------------------------------------------------------------------------
  # Energy Conversions
  # ---------------------------------------------------------------------------

  @doc "Converts electron-volts to joules: J = eV × 1.602176634×10⁻¹⁹."
  @spec ev_to_j(number) :: float
  def ev_to_j(ev), do: ev * elementary_charge()

  @doc "Converts joules to electron-volts."
  @spec j_to_ev(number) :: float
  def j_to_ev(j), do: j / elementary_charge()

  @doc "Converts ergs to joules: J = erg × 10⁻⁷."
  @spec erg_to_j(number) :: float
  def erg_to_j(erg), do: erg * 1.0e-7

  @doc "Converts joules to ergs: erg = J × 10⁷."
  @spec j_to_erg(number) :: float
  def j_to_erg(j), do: j * 1.0e7

  # ---------------------------------------------------------------------------
  # Mass Conversions
  # ---------------------------------------------------------------------------

  @doc "Converts kilograms to solar masses."
  @spec kg_to_solar_masses(number) :: float
  def kg_to_solar_masses(kg), do: kg / solar_mass()

  @doc "Converts solar masses to kilograms."
  @spec solar_masses_to_kg(number) :: float
  def solar_masses_to_kg(msun), do: msun * solar_mass()

  # ---------------------------------------------------------------------------
  # Time Conversions
  # ---------------------------------------------------------------------------

  @doc "Returns the number of seconds in a Julian year (365.25 days)."
  @spec seconds_per_year() :: float
  def seconds_per_year, do: 3.15576e7

  @doc "Converts Julian years to seconds."
  @spec yr_to_s(number) :: float
  def yr_to_s(yr), do: yr * seconds_per_year()

  @doc "Converts gigayears (10⁹ yr) to seconds."
  @spec gyr_to_s(number) :: float
  def gyr_to_s(gyr), do: gyr * 1.0e9 * seconds_per_year()

  # ---------------------------------------------------------------------------
  # Temperature Conversions
  # ---------------------------------------------------------------------------

  @doc "Converts Kelvin to Celsius: °C = K − 273.15."
  @spec k_to_celsius(number) :: float
  def k_to_celsius(k), do: k - 273.15

  @doc "Converts Celsius to Kelvin: K = °C + 273.15."
  @spec celsius_to_k(number) :: float
  def celsius_to_k(c), do: c + 273.15

  # ---------------------------------------------------------------------------
  # Sexagesimal ↔ Decimal
  # ---------------------------------------------------------------------------

  @doc """
  Converts sexagesimal hours (hh:mm:ss) to decimal degrees.

  ## Parameters
    - h: Hours
    - m: Minutes
    - s: Seconds

  ## Returns
    Decimal degrees

  ## Examples
      iex> Notation.hms_to_deg(12, 30, 0)
      187.5
  """
  @spec hms_to_deg(number, number, number) :: float
  def hms_to_deg(h, m, s) do
    (h + m / 60.0 + s / 3600.0) * 15.0
  end

  @doc """
  Converts sexagesimal degrees (dd:mm:ss) to decimal degrees.

  ## Parameters
    - d: Degrees (negative for south/west)
    - m: Arcminutes
    - s: Arcseconds

  ## Returns
    Decimal degrees

  ## Examples
      iex> Notation.dms_to_deg(-30, 15, 30)
      -30.258333333333333
  """
  @spec dms_to_deg(number, number, number) :: float
  def dms_to_deg(d, m, s) do
    sign = if d < 0, do: -1, else: 1
    sign * (abs(d) + m / 60.0 + s / 3600.0)
  end

  @doc """
  Converts decimal degrees to sexagesimal degrees (dd, mm, ss).

  ## Returns
    {degrees, arcminutes, arcseconds}

  ## Examples
      iex> elem(Notation.deg_to_dms(-30.258333), 0)
      -30
  """
  @spec deg_to_dms(number) :: {integer, integer, float}
  def deg_to_dms(deg) do
    sign = if deg < 0, do: -1, else: 1
    abs_d = abs(deg)
    d = trunc(abs_d)
    rem_m = (abs_d - d) * 60
    m = trunc(rem_m)
    s = (rem_m - m) * 60
    {sign * d, m, s}
  end

  @doc """
  Converts decimal degrees to sexagesimal hours (hh, mm, ss).

  ## Returns
    {hours, minutes, seconds}

  ## Examples
      iex> Notation.deg_to_hms(187.5)
      {12, 30, 0.0}
  """
  @spec deg_to_hms(number) :: {integer, integer, float}
  def deg_to_hms(deg) do
    h_dec = deg / 15.0
    h = trunc(h_dec)
    rem_m = (h_dec - h) * 60
    m = trunc(rem_m)
    s = (rem_m - m) * 60
    {h, m, s}
  end

  # ---------------------------------------------------------------------------
  # Julian Date
  # ---------------------------------------------------------------------------

  @doc """
  Converts a Gregorian calendar date to Julian Date.

  Valid for dates from 1582-10-15 onwards.

  ## Parameters
    - year, month, day: Gregorian calendar date
    - hour, minute, second: UT time (default: noon = 12:00:00)

  ## Returns
    Julian Date

  ## Examples
      iex> Notation.calendar_to_jd(2000, 1, 1.5)
      2451545.0
  """
  @spec calendar_to_jd(number, number, number, number, number, number) :: float
  def calendar_to_jd(year, month, day, hour \\ 12, minute \\ 0, second \\ 0) do
    y = if month <= 2, do: year - 1, else: year
    m = if month <= 2, do: month + 12, else: month

    a = trunc(y / 100)
    b = 2 - a + trunc(a / 4)

    day_frac = day + (hour + minute / 60.0 + second / 3600.0) / 24.0
    trunc(365.25 * (y + 4716)) + trunc(30.6001 * (m + 1)) + day_frac + b - 1524.5
  end

  @doc """
  Converts Julian Date to Modified Julian Date: MJD = JD - 2400000.5

  ## Examples
      iex> Notation.jd_to_mjd(2_451_545.0) |> Float.round(1)
      51544.5
  """
  @spec jd_to_mjd(number) :: float
  def jd_to_mjd(jd), do: jd - 2_400_000.5

  @doc """
  Converts Modified Julian Date to Julian Date.

  ## Examples
      iex> Notation.mjd_to_jd(51_544.5) |> Float.round(1)
      2451545.0
  """
  @spec mjd_to_jd(number) :: float
  def mjd_to_jd(mjd), do: mjd + 2_400_000.5

  @doc "Returns the Julian Date of the J2000.0 epoch: 2451545.0 (2000 Jan 1.5 TT)."
  @spec j2000() :: float
  def j2000, do: 2_451_545.0

  @doc "Julian centuries elapsed since J2000.0 for a given Julian Date."
  @spec julian_centuries(number) :: float
  def julian_centuries(jd), do: (jd - j2000()) / 36_525.0
end
