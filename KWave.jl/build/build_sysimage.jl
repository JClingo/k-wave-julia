"""
KWave.jl — PackageCompiler.jl sysimage build script
=====================================================
Compiles KWave.jl into a Julia system image that eliminates JIT compilation
overhead on first use.  Benchmark data shows 1487× overhead for small 1D
scenarios and 40× for small 2D scenarios on a cold start — this script
removes that latency.

Unlike `build_standalone.jl` (which produces a self-contained app binary),
this produces a `.dylib` / `.so` sysimage that is loaded by any Julia session:

    julia --sysimage build/kwave-sysimage.so -t auto

Usage:
    julia --project=.. build/build_sysimage.jl [--output PATH]

Prerequisites (one-time install):
    ] add PackageCompiler

Expected build time: 3–10 minutes depending on hardware.

Tips:
  - Use `-t auto` when building to parallelise compilation.
  - The resulting sysimage bakes in FFTW plan caches for the precompile
    workload shapes — other grid sizes will still plan dynamically on first
    use but remain JIT-free for all KWave Julia code.
"""

using Pkg

try
    @eval using PackageCompiler
catch
    @error "PackageCompiler.jl is required. Install it with: ] add PackageCompiler"
    exit(1)
end

# ─────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────

output_path = joinpath(@__DIR__, "kwave-sysimage.so")
for i in eachindex(ARGS)
    if ARGS[i] == "--output" && i < length(ARGS)
        output_path = ARGS[i + 1]
    end
end

project_dir = dirname(@__DIR__)

# ─────────────────────────────────────────────────────────────
# Precompile execution script
# ─────────────────────────────────────────────────────────────
# This script is executed during image compilation.  The more code paths it
# exercises, the more JIT work is baked into the sysimage.

precompile_script = joinpath(@__DIR__, "precompile_sysimage.jl")

open(precompile_script, "w") do io
    write(io, """
    using KWave
    using FFTW
    FFTW.set_num_threads(Threads.nthreads())

    # ── 1D ──────────────────────────────────────────────────────
    let
        kgrid = KWaveGrid(256, 1e-4)
        make_time!(kgrid, 1500.0; t_end=2e-6)
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = zeros(256); p0[128] = 1.0
        source = KWaveSource(p0=p0)
        mask = falses(256); mask[64] = true
        sensor = KWaveSensor(mask=mask, record=[:p])
        kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        println("precompile: 1D F64 done")
    end

    # ── 1D absorbing ────────────────────────────────────────────
    let
        kgrid = KWaveGrid(256, 1e-4)
        make_time!(kgrid, 1500.0; t_end=2e-6)
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0,
                             alpha_coeff=0.75, alpha_power=1.5)
        p0 = zeros(256); p0[128] = 1.0
        source = KWaveSource(p0=p0)
        mask = falses(256); mask[64] = true
        sensor = KWaveSensor(mask=mask, record=[:p])
        kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        println("precompile: 1D absorbing done")
    end

    # ── 2D ──────────────────────────────────────────────────────
    let
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        make_time!(kgrid, 1500.0; t_end=1e-6)
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_disc(64, 64, 32, 32, 4))
        source = KWaveSource(p0=p0)
        mask = falses(64, 64); mask[1, :] .= true
        sensor = KWaveSensor(mask=mask, record=[:p])
        kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        println("precompile: 2D F64 done")
    end

    # ── 2D F32 ──────────────────────────────────────────────────
    let
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        make_time!(kgrid, 1500.0; t_end=1e-6)
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_disc(64, 64, 32, 32, 4))
        source = KWaveSource(p0=p0)
        mask = falses(64, 64); mask[1, :] .= true
        sensor = KWaveSensor(mask=mask, record=[:p])
        kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false, data_cast=Float32)
        println("precompile: 2D F32 done")
    end

    # ── 3D ──────────────────────────────────────────────────────
    let
        kgrid = KWaveGrid(32, 1e-4, 32, 1e-4, 32, 1e-4)
        make_time!(kgrid, 1500.0; t_end=5e-7)
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_ball(32, 32, 32, 16, 16, 16, 3))
        source = KWaveSource(p0=p0)
        mask = falses(32, 32, 32); mask[1, 16, 16] = true
        sensor = KWaveSensor(mask=mask, record=[:p])
        kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        println("precompile: 3D F64 done")
    end

    # ── 2D heterogeneous ────────────────────────────────────────
    let
        using Random
        rng = MersenneTwister(42)
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        make_time!(kgrid, 1500.0; t_end=1e-6)
        c_map = 1500.0 .* (1.0 .+ 0.1 .* randn(rng, 64, 64))
        medium = KWaveMedium(sound_speed=c_map, density=1000.0)
        p0 = Float64.(make_disc(64, 64, 32, 32, 4))
        source = KWaveSource(p0=p0)
        mask = falses(64, 64); mask[1, :] .= true
        sensor = KWaveSensor(mask=mask, record=[:p])
        kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        println("precompile: 2D heterogeneous done")
    end

    # ── utility functions ────────────────────────────────────────
    make_disc(32, 32, 16, 16, 5)
    make_ball(16, 16, 16, 8, 8, 8, 3)
    tone_burst(1 / 2e-8, 1e6, 3)
    println("precompile: utilities done")
    """)
end

# ─────────────────────────────────────────────────────────────
# Build
# ─────────────────────────────────────────────────────────────

println("Building KWave sysimage...")
println("  Project:  $project_dir")
println("  Output:   $output_path")
println("  This may take several minutes...")
println()

create_sysimage(
    [:KWave, :FFTW];
    sysimage_path=output_path,
    precompile_execution_file=precompile_script,
    project=project_dir,
)

println()
println("Build complete!")
println()
println("Use the sysimage to eliminate JIT overhead:")
println("  julia --sysimage $output_path -t auto")
println()
println("Benchmark with sysimage:")
println("  julia --sysimage $output_path -t auto --project=benchmarks/julia benchmarks/julia/run_benchmarks.jl")
