"""
KWave.jl GPU Performance Benchmark Suite
==========================================
Benchmarks Metal (Apple Silicon) and CUDA GPU backends against the CPU Float32
baseline using matched scenarios.

Usage:
    # Metal (Mac)
    julia --project=benchmarks/julia benchmarks/julia/run_benchmarks_gpu.jl metal

    # CUDA (NVIDIA)
    julia --project=benchmarks/julia benchmarks/julia/run_benchmarks_gpu.jl cuda

    # Auto-detect
    julia --project=benchmarks/julia benchmarks/julia/run_benchmarks_gpu.jl

Results written to:
    benchmarks/results/julia/gpu_simulations.csv
"""

using Pkg
Pkg.activate(@__DIR__)

# ─────────────────────────────────────────────────────────────
# GPU backend detection / loading
# ─────────────────────────────────────────────────────────────

const GPU_ARG = get(ARGS, 1, "auto")

global gpu_type   = :none
global gpu_to_f32 = identity   # array → GPU Float32
global gpu_to_cpu = identity   # GPU array → CPU

if GPU_ARG in ("metal", "auto")
    try
        using Metal
        if Metal.functional()
            global gpu_type   = :metal
            global gpu_to_f32 = x -> x isa AbstractArray ? MtlArray{Float32}(x) :
                                     x isa Real           ? Float32(x)           : x
            global gpu_to_cpu = x -> x isa MtlArray ? Array(x) : x
            println("Backend: Metal (Apple GPU)")
        else
            println("Metal loaded but not functional.")
        end
    catch e
        println("Metal unavailable: $e")
        GPU_ARG == "metal" && exit(1)
    end
end

if gpu_type == :none && GPU_ARG in ("cuda", "auto")
    try
        using CUDA
        if CUDA.functional()
            global gpu_type   = :cuda
            global gpu_to_f32 = x -> x isa AbstractArray ? CuArray{Float32}(x) :
                                     x isa Real           ? Float32(x)          : x
            global gpu_to_cpu = x -> x isa CuArray ? Array(x) : x
            println("Backend: CUDA ($(CUDA.name(CUDA.device())))")
        else
            println("CUDA loaded but not functional.")
        end
    catch e
        println("CUDA unavailable: $e")
        GPU_ARG == "cuda" && exit(1)
    end
end

if gpu_type == :none
    println("No functional GPU backend found.  Exiting.")
    exit(0)
end

using KWave
using Statistics
using Random
using Printf
using Dates
using CSV
using DataFrames

const N_WARMUP  = 2   # GPU warmup matters more (shader/kernel compilation)
const N_SAMPLES = parse(Int, get(ENV, "BENCH_SAMPLES", "5"))
const RNG_SEED  = 42
const DX        = 1e-4
const C0        = 1500.0
const RHO0      = 1000.0
const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "julia")
mkpath(RESULTS_DIR)

println("=" ^ 65)
println("KWave.jl GPU Benchmark  ($(gpu_type))")
println("Julia: $(VERSION)  Date: $(now())")
println("N_WARMUP=$N_WARMUP  N_SAMPLES=$N_SAMPLES")
println("=" ^ 65)

# ─────────────────────────────────────────────────────────────
# GPU medium / source conversion
# ─────────────────────────────────────────────────────────────

function to_gpu_medium(medium::KWaveMedium)
    _g = gpu_to_f32
    return KWaveMedium{Float32}(
        _g(medium.sound_speed),
        _g(medium.density),
        medium.alpha_coeff  === nothing ? nothing : _g(medium.alpha_coeff),
        medium.alpha_power  === nothing ? nothing : Float32(medium.alpha_power),
        medium.alpha_mode,
        medium.BonA === nothing ? nothing : _g(medium.BonA),
    )
end

function to_gpu_source(source::KWaveSource)
    _g  = gpu_to_f32
    _gb = x -> x === nothing ? nothing :
               (x isa AbstractArray{Bool} ? (gpu_type == :cuda ? CuArray(x) : MtlArray(x)) :
                _g(x))
    return KWaveSource{Float32}(
        source.p0     === nothing ? nothing : _g(source.p0),
        source.p_mask === nothing ? nothing : _gb(source.p_mask),
        source.p      === nothing ? nothing : _g(source.p),
        source.p_mode,
        source.u_mask === nothing ? nothing : _gb(source.u_mask),
        source.ux     === nothing ? nothing : _g(source.ux),
        source.uy     === nothing ? nothing : _g(source.uy),
        source.uz     === nothing ? nothing : _g(source.uz),
        source.u_mode,
    )
end

# ─────────────────────────────────────────────────────────────
# GPU scenarios (subset of CPU suite — focus on compute-heavy cases)
# ─────────────────────────────────────────────────────────────
# Each entry: (label, dim, Nx, Ny, Nz, t_end, pml_size)

GPU_SCENARIOS = [
    # 2D — meaningful sizes for GPU
    ("2d_small",   2,  64,  64,   1,  4e-6, 20),
    ("2d_medium",  2, 256, 256,   1,  4e-6, 20),
    ("2d_large",   2, 512, 512,   1,  2e-6, 20),
    # 3D — GPU shines here
    ("3d_small",   3,  32,  32,  32,  2e-6,  8),
    ("3d_medium",  3,  64,  64,  64,  2e-6, 12),
    ("3d_large",   3, 128, 128, 128,  1e-6, 20),
]

# ─────────────────────────────────────────────────────────────
# Build simulation components (CPU arrays, converted later)
# ─────────────────────────────────────────────────────────────

function build_sim_cpu(label, dim, Nx, Ny, Nz, t_end, pml_size)
    medium_cpu = KWaveMedium(sound_speed=C0, density=RHO0)

    if dim == 2
        kgrid  = KWaveGrid(Nx, DX, Ny, DX)
        make_time!(kgrid, C0; t_end=t_end)
        p0     = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5))
        source_cpu = KWaveSource(p0=p0)
        smask  = falses(Nx, Ny)
        smask[Nx÷4,   Ny÷2] = true
        smask[3*Nx÷4, Ny÷2] = true
    else  # 3D
        kgrid  = KWaveGrid(Nx, DX, Ny, DX, Nz, DX)
        make_time!(kgrid, C0; t_end=t_end)
        p0     = Float64.(make_ball(Nx, Ny, Nz, Nx÷2, Ny÷2, Nz÷2, 3))
        source_cpu = KWaveSource(p0=p0)
        smask  = falses(Nx, Ny, Nz)
        smask[Nx÷4,   Ny÷2, Nz÷2] = true
        smask[3*Nx÷4, Ny÷2, Nz÷2] = true
    end

    sensor = KWaveSensor(mask=smask, record=[:p])
    return kgrid, medium_cpu, source_cpu, sensor
end

# ─────────────────────────────────────────────────────────────
# Timing harness
# ─────────────────────────────────────────────────────────────

function timed_runs(f; n_warmup=N_WARMUP, n_samples=N_SAMPLES)
    # First call: JIT + GPU kernel compilation
    first_stats = @timed f()
    t_first     = first_stats.time

    for _ in 2:n_warmup; f(); end

    times  = Vector{Float64}(undef, n_samples)
    allocs = Vector{Int}(undef,    n_samples)
    for i in 1:n_samples
        s = @timed f()
        times[i]  = s.time
        allocs[i] = s.bytes
    end

    return (
        first_run_s = t_first,
        min_s       = minimum(times),
        median_s    = median(times),
        mean_s      = mean(times),
        max_s       = maximum(times),
        std_s       = std(times),
        median_alloc = Int(round(median(allocs))),
    )
end

# ─────────────────────────────────────────────────────────────
# CPU Float32 baseline (same scenarios, for direct comparison)
# ─────────────────────────────────────────────────────────────

println("\n─── CPU Float32 baseline ───")
cpu_rows = []

for (label, dim, Nx, Ny, Nz, t_end, pml_size) in GPU_SCENARIOS
    kgrid, medium_cpu, source_cpu, sensor = build_sim_cpu(label, dim, Nx, Ny, Nz, t_end, pml_size)
    Nt = kgrid.Nt[]

    f = () -> kspace_first_order(kgrid, medium_cpu, source_cpu, sensor;
                                  smooth_p0=false, plot_sim=false,
                                  pml_size=pml_size, data_cast=Float32)

    @printf("  %-14s  Nx=%-4d Nt=%-4d  ", label, Nx, Nt)
    flush(stdout)
    try
        r = timed_runs(f)
        @printf("first=%6.3fs  median=%6.3fs\n", r.first_run_s, r.median_s)
        push!(cpu_rows, (
            timestamp="$(now())", language="julia", backend="cpu_f32",
            scenario=label, dim=dim, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt,
            medium_type="homogeneous", r...
        ))
    catch e; println(" ERROR: $e"); end
end

# ─────────────────────────────────────────────────────────────
# GPU Float32
# ─────────────────────────────────────────────────────────────

println("\n─── GPU Float32  ($(gpu_type)) ───")
gpu_rows = []

for (label, dim, Nx, Ny, Nz, t_end, pml_size) in GPU_SCENARIOS
    kgrid, medium_cpu, source_cpu, sensor = build_sim_cpu(label, dim, Nx, Ny, Nz, t_end, pml_size)
    Nt = kgrid.Nt[]

    medium_gpu = to_gpu_medium(medium_cpu)
    source_gpu = to_gpu_source(source_cpu)

    f = () -> kspace_first_order(kgrid, medium_gpu, source_gpu, sensor;
                                  smooth_p0=false, plot_sim=false,
                                  pml_size=pml_size, data_cast=Float32)

    @printf("  %-14s  Nx=%-4d Nt=%-4d  ", label, Nx, Nt)
    flush(stdout)
    try
        r = timed_runs(f)
        @printf("first=%6.3fs  median=%6.3fs\n", r.first_run_s, r.median_s)
        push!(gpu_rows, (
            timestamp="$(now())", language="julia", backend="$(gpu_type)_f32",
            scenario=label, dim=dim, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt,
            medium_type="homogeneous", r...
        ))
    catch e
        println(" ERROR: $e")
        Base.showerror(stdout, e, catch_backtrace()); println()
    end
end

# ─────────────────────────────────────────────────────────────
# Print speedup summary
# ─────────────────────────────────────────────────────────────

println("\n─── GPU Speedup Summary  (cpu_f32 → $(gpu_type)_f32) ───")
@printf("  %-14s  %10s  %10s  %8s\n", "scenario", "cpu_f32", "$(gpu_type)_f32", "speedup")
println("  " * "─"^50)

cpu_dict = Dict(r.scenario => r.median_s for r in cpu_rows)
for r in gpu_rows
    if haskey(cpu_dict, r.scenario)
        speedup = cpu_dict[r.scenario] / r.median_s
        @printf("  %-14s  %8.3f s  %8.3f s  %6.2fx\n",
                r.scenario, cpu_dict[r.scenario], r.median_s, speedup)
    end
end

# ─────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────

df = DataFrame(vcat(cpu_rows, gpu_rows))
path = joinpath(RESULTS_DIR, "gpu_simulations.csv")
CSV.write(path, df)
println("\nSaved: $path")
println("Done.  $(now())")
