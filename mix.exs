defmodule AstroEquations.MixProject do
  use Mix.Project

  @version "0.2.0"
  @source_url "https://github.com/iamkanishka/astroeqex"

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
      aliases: aliases(),
      test_coverage: [tool: ExCoveralls],
      preferred_cli_env: [
        coveralls: :test,
        "coveralls.detail": :test,
        "coveralls.post": :test,
        "coveralls.html": :test,
        "coveralls.github": :test
      ],
      dialyzer: dialyzer(),
      elixirc_options: [warnings_as_errors: Mix.env() == :prod]
    ]
  end

  def application do
    [extra_applications: [:logger]]
  end

  # ---------------------------------------------------------------------------
  # Dependencies
  # ---------------------------------------------------------------------------

  defp deps do
    [
      # Documentation
      {:ex_doc, "~> 0.34", only: :dev, runtime: false},

      # Static analysis
      {:dialyxir, "~> 1.4", only: [:dev, :test], runtime: false},
      {:credo, "~> 1.7", only: [:dev, :test], runtime: false},

      # Test coverage
      {:excoveralls, "~> 0.18", only: :test, runtime: false},

      # Property-based testing
      {:stream_data, "~> 1.1", only: [:dev, :test]}
    ]
  end

  # ---------------------------------------------------------------------------
  # Hex metadata
  # ---------------------------------------------------------------------------

  defp description do
    """
    A comprehensive Elixir library of scientific and astronomical equations spanning
    astrophysics, physics, mathematics, and statistics. 696+ functions covering 23 modules —
    all in SI units with full @spec, @type, @doc, and input guard coverage.
    """
  end

  defp package do
    [
      name: "astroequations",
      licenses: ["MIT"],
      links: %{
        "GitHub" => @source_url,
        "Changelog" => "#{@source_url}/blob/master/CHANGELOG.md"
      },
      maintainers: ["Kanishka Naik"],
      files: ~w(lib .formatter.exs mix.exs README.md CHANGELOG.md LICENSE)
    ]
  end

  # ---------------------------------------------------------------------------
  # ExDoc
  # ---------------------------------------------------------------------------

  defp docs do
    [
      main: "readme",
      name: "AstroEquations",
      source_url: @source_url,
      source_ref: "v#{@version}",
      api_reference: true,
      extras: [
        "README.md": [title: "Overview"],
        "CHANGELOG.md": [title: "Changelog"]
      ],
      groups_for_modules: [
        "Astrophysics & Astronomy": [
          AstroEquations.AstrophysicsAndAstronomy.Astrometry,
          AstroEquations.AstrophysicsAndAstronomy.BlackHole,
          AstroEquations.AstrophysicsAndAstronomy.Galaxies,
          AstroEquations.AstrophysicsAndAstronomy.Instrumentation,
          AstroEquations.AstrophysicsAndAstronomy.Stars
        ],
        Physics: [
          AstroEquations.Physics.Electromagnetism,
          AstroEquations.Physics.Energy,
          AstroEquations.Physics.Forces,
          AstroEquations.Physics.GeneralRelativity,
          AstroEquations.Physics.Materials,
          AstroEquations.Physics.Motion,
          AstroEquations.Physics.NewtonGravity,
          AstroEquations.Physics.Oscillations,
          AstroEquations.Physics.QuantumMechanics,
          AstroEquations.Physics.SpecialRelativity,
          AstroEquations.Physics.Thermodynamics,
          AstroEquations.Physics.Waves
        ],
        Mathematics: [
          AstroEquations.Mathematics.Calculus,
          AstroEquations.Mathematics.Geometry,
          AstroEquations.Mathematics.Notation,
          AstroEquations.Mathematics.Trigonometry
        ],
        Statistics: [
          AstroEquations.Statistics.StandardDeviation,
          AstroEquations.Statistics.Variance
        ],
        Deprecated: [
          AstroEquations.AstrophysicsAndAstronomy.Instrumentaion,
          AstroEquations.Mathematics.Trignometry
        ]
      ]
    ]
  end

  # ---------------------------------------------------------------------------
  # Dialyzer
  # ---------------------------------------------------------------------------

  defp dialyzer do
    [
      plt_add_apps: [:mix, :ex_unit],
      plt_local_path: "priv/plts/project.plt",
      plt_core_path: "priv/plts/core.plt",
      flags: [
        :error_handling,
        :missing_return,
        :unmatched_returns,
        :unknown
      ],
      ignore_warnings: ".dialyzer_ignore.exs",
      list_unused_filters: true
    ]
  end

  # ---------------------------------------------------------------------------
  # Aliases
  # ---------------------------------------------------------------------------

  defp aliases do
    [
      # Setup
      setup: ["deps.get", "compile"],

      # Full quality check — run before every commit
      check: [
        "format --check-formatted",
        "credo --strict",
        "dialyzer"
        # "test --cover"
      ],

      # Auto-fix formatting then strict lint
      lint: ["format", "credo --strict"],

      # Rebuild Dialyzer PLTs from scratch after dep upgrades
      "dialyzer.clean": ["cmd rm -rf priv/plts", "dialyzer --plt"],

      # Generate and open docs locally
      "docs.open": ["docs", "cmd open doc/index.html"],

      # CI pipeline — no interactive prompts
      ci: [
        "format --check-formatted",
        "deps.unlock --check-unused",
        "credo --strict",
        "dialyzer --quiet",
        "test --cover"
      ]
    ]
  end
end
