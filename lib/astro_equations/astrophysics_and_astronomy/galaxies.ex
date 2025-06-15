defmodule AstroEquations.AstrophysicsAndAstronomy.Galaxies do
  @moduledoc """
  A collection of functions for calculating galaxy properties including:
  - Hubble elliptical galaxy classification
  - Sérsic light profiles
  - Stellar density distributions in disk galaxies

  All calculations use consistent astronomical units.
  """

  @doc """
  Classifies an elliptical galaxy based on its apparent flattening.

  ## Parameters
    - a: Semi-major axis length
    - b: Semi-minor axis length

  ## Returns
    Hubble classification type as a string (E0 to E7)

  ## Examples
      iex> GalaxyProperties.hubble_classify(10, 10)
      "E0"
      iex> GalaxyProperties.hubble_classify(10, 5)
      "E5"
  """
  @spec hubble_classify(number, number) :: String.t()
  def hubble_classify(a, b) do
    ratio = (a - b) / a
    type = round(ratio * 10)
    "E#{min(max(type, 0), 7)}"  # Clamped between E0 and E7
  end

  @doc """
  Calculates the surface brightness using the Sérsic profile.

  ## Parameters
    - i0: Central surface brightness
    - r: Radius from center
    - r_e: Effective radius (half-light radius)
    - n: Sérsic index (shape parameter)

  ## Returns
    Surface brightness at radius r

  ## Examples
      iex> GalaxyProperties.sersic_profile(100, 1, 1, 1) |> Float.round(4)
      100.0
      iex> GalaxyProperties.sersic_profile(100, 2, 1, 4) |> Float.round(4)
      13.5335
  """
  @spec sersic_profile(number, number, number, number) :: float
  def sersic_profile(i0, r, r_e, n) do
    # Approximation of bn for n > 0.36
    bn = 2 * n - 1/3 + 0.009876 / n
    i0 * :math.exp(-bn * (:math.pow(r / r_e, 1/n) - 1))
  end

  @doc """
  Calculates the stellar density in a disk galaxy using exponential profiles.

  ## Parameters
    - p0: Central density
    - r: Cylindrical radius from center
    - z: Height above the plane
    - h: Radial scale length
    - z0: Vertical scale height

  ## Returns
    Stellar density at position (R,z)

  ## Examples
      iex> GalaxyProperties.disk_density(1, 0, 0, 1, 1) |> Float.round(4)
      1.0
      iex> GalaxyProperties.disk_density(1, 1, 1, 1, 1) |> Float.round(4)
      0.1353
  """
  @spec disk_density(number, number, number, number, number) :: float
  def disk_density(p0, r, z, h, z0) do
    p0 * :math.exp(-abs(z)/z0) * :math.exp(-r/h)
  end

end
