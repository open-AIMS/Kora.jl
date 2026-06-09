using Documenter, DocumenterVitepress, Kora

makedocs(;
    sitename="Kora.jl",
    modules=[Kora],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="https://github.com/ConnectedSystems/Kora.jl",
        devbranch="main",
        devurl="dev"
    ),
    warnonly=Symbol[],
    checkdocs=:exports,
    pages=[
        "Start Here" => "what-can-kora-tell-me.md",
        "Getting Started" => "getting-started.md",
        "Model Overview" => "model-overview.md",
        "Concepts" => [
            "Decision Support Under Uncertainty" => "concepts/decision-support-under-uncertainty.md"
        ],
        "Tutorials" => [
            "Running Simulations" => "tutorials/running-simulations.md",
            "Visualization" => "tutorials/visualization.md",
            "Ensemble Analysis" => "tutorials/ensemble-analysis.md",
            "Fitting from EcoRRAP" => "tutorials/fitting-from-ecorrap.md",
            "Restoration Scenarios" => "tutorials/restoration-scenarios.md"
        ],
        "Calibration" => [
            "Model Calibration" => "calibration/model-calibration.md",
            "Ensemble Assessment" => "calibration/ensemble-assessment.md"
        ],
        "API Reference" => [
            "Reef State" => "reference/api-reef-state.md",
            "Simulation" => "reference/api-simulation.md",
            "Coral Models" => "reference/api-models.md",
            "Model I/O" => "reference/api-interface.md",
            "Coral Dynamics" => "reference/api-coral-dynamics.md",
            "Metrics" => "reference/api-metrics.md"
        ],
        "Background" => [
            "Coral Biology" => "concepts/coral-biology.md"
        ],
        "Glossary" => "glossary.md",
        "Contributing" => "contributing.md"
    ]
)

deploydocs(;
    repo="github.com/ConnectedSystems/Kora.jl",
    target="build",
    branch="gh-pages",
    devbranch="main",
    push_preview=true
)
