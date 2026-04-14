using Documenter
using KWave

DocMeta.setdocmeta!(KWave, :DocTestSetup, :(using KWave); recursive=true)

makedocs(
    sitename = "KWave.jl",
    authors  = "Joshua Clingo and contributors",
    repo     = "https://github.com/JClingo/k-wave-julia/blob/{commit}{path}#{line}",
    format   = Documenter.HTML(
        prettyurls    = get(ENV, "CI", nothing) == "true",
        canonical     = "https://jclingo.github.io/KWave.jl",
        collapselevel = 2,
        warn_outdated = true,
        size_threshold_ignore = ["api/solvers.md", "api/geometry.md"],
    ),
    modules   = [KWave],
    warnonly  = [:missing_docs, :cross_references],
    checkdocs = :exports,
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Theory"          => "theory.md",
        "Manual" => [
            "Grid"                => "manual/grid.md",
            "Medium"              => "manual/medium.md",
            "Source & Sensor"     => "manual/source_sensor.md",
            "Running Simulations" => "manual/simulation.md",
            "PML Configuration"   => "manual/pml.md",
            "Visualization"       => "manual/visualization.md",
            "GPU Acceleration"    => "manual/gpu.md",
            "Unitful Integration" => "manual/unitful.md",
            "Python Interop"      => "manual/python_interop.md",
            "HDF5 File Format"    => "manual/hdf5_format.md",
        ],
        "Examples" => [
            "Gallery"                        => "examples/index.md",
        ],
        "API Reference" => [
            "Overview"             => "api/index.md",
            "Solvers"              => "api/solvers.md",
            "Grid"                 => "api/grid.md",
            "Medium"               => "api/medium.md",
            "Source & Sensor"      => "api/source_sensor.md",
            "Signal Processing"    => "api/signal.md",
            "Geometry"             => "api/geometry.md",
            "Arrays"               => "api/array.md",
            "Transducer"           => "api/transducer.md",
            "Reconstruction"       => "api/reconstruction.md",
            "Visualization"        => "api/visualization.md",
            "Utilities"            => "api/utils.md",
            "I/O"                  => "api/io.md",
            "Reference Solutions"  => "api/reference.md",
            "Python Interop"       => "api/interop.md",
        ],
    ],
)

deploydocs(
    repo      = "github.com/JClingo/k-wave-julia",
    target    = "build",
    branch    = "gh-pages",
    devbranch = "main",
)
