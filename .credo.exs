%{
  configs: [
    %{
      name: "default",
      files: %{
        included: ["lib/", "src/", "test/"],
        excluded: []
      },
      checks: [
        {Credo.Check.Readability.ModuleDoc, false}
      ]
    }
  ]
}
