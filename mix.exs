defmodule AstroEquations.MixProject do
  use Mix.Project

  @version "0.2.0"
  @source_url "https://github.com/your-org/astroequations"

  def project do
    [
      app: :astroequations,
      version: @version,
      elixir: "~> 1.18",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      description: description(),
      package: package(),
      docs: docs(),
      dialyzer: [
        plt_add_apps: [:mix],
        flags: [:error_handling, :underspecs]
      ]
    ]
  end

  def application do
    [extra_applications: [:logger]]
  end

  defp deps do
    [
      {:ex_doc, "~> 0.34", only: :dev, runtime: false},
      {:dialyxir, "~> 1.4", only: [:dev, :test], runtime: false},
      {:credo, "~> 1.7", only: [:dev, :test], runtime: false}
    ]
  end

  defp description do
    """
    A comprehensive Elixir library of scientific and astronomical equations spanning
    physics, astrophysics, mathematics, and statistics. 627+ functions in SI units.
    """
  end

  defp package do
    [
      name: "astroequations",
      licenses: ["MIT"],
      links: %{"GitHub" => @source_url},
      maintainers: ["Your Name"],
      files: ~w(lib .formatter.exs mix.exs README.md CHANGELOG.md LICENSE)
    ]
  end

  defp docs do
    [
      main: "readme",
      name: "AstroEquations",
      source_url: @source_url,
      extras: ["README.md", "CHANGELOG.md"]
    ]
  end
end
