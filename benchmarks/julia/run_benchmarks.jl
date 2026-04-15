"""
KWave.jl CPU Performance Benchmark Suite
=========================================
Benchmarks full simulations (1D / 2D / 3D) and FFT component throughput
for both Float64 and Float32 on CPU.

Usage:
    # Run with all CPU threads (STRONGLY RECOMMENDED — gives 4-8x FFT speedup)
    julia -t auto --project=benchmarks/julia benchmarks/julia/run_benchmarks.jl

    # Explicit thread count (replace N with desired number)
    julia -t N --project=benchmarks/julia benchmarks/julia/run_benchmarks.jl

    # Override sample count:
    BENCH_SAMPLES=10 julia -t auto --project=benchmarks/julia benchmarks/julia/run_benchmarks.jl

Results are written to:
    benchmarks/results/julia/full_simulations.csv
    benchmarks/results/julia/fft_components.csv

Note: FFTW multi-threading is set to Threads.nthreads() inside create_fft_plans.
Running with a single Julia thread (no -t flag) pins FFT to 1 thread and will
show severely degraded 2D/3D performance vs MATLAB.
"""

using Pkg
Pkg.activate(@__DIR__)

# Load Apple Accelerate before KWave so the KWaveAppleAccelerateExt extension
# is triggered.  On Apple Silicon this replaces FFTW (which ships without ARM
# SIMD flags) with vDSP for 1D FFTs — critical for 1D scenario performance.
try
    using AppleAccelerate
catch e
    @info "AppleAccelerate not available: $(e). Using FFTW for all FFT operations."
end

using KWave
using FFTW
using LinearAlgebra: mul!
using BenchmarkTools
using Statistics
using Random
using Printf
using Dates
using CSV
using DataFrames

# ─────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────

const N_WARMUP  = 1      # runs before timing begins (warms JIT + FFTW planner cache)
const N_SAMPLES = parse(Int, get(ENV, "BENCH_SAMPLES", "5"))
const RNG_SEED  = 42
const DX        = 1e-4   # 0.1 mm — common to all scenarios
const C0        = 1500.0 # m/s — water
const RHO0      = 1000.0 # kg/m³

const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "julia")
mkpath(RESULTS_DIR)

# Apply FFTW threading at startup so the planning stage also uses multiple threads.
# create_fft_plans() also calls this, but setting it here ensures the benchmark
# header prints the correct thread count before any plan is created.
FFTW.set_num_threads(Threads.nthreads())

println("=" ^ 65)
println("KWave.jl CPU Benchmark Suite")
println("Julia:        $(VERSION)")
println("Julia threads: $(Threads.nthreads())")
println("FFTW threads:  $(FFTW.get_num_threads())")
println("Date:         $(now())")
println("N_WARMUP=$N_WARMUP  N_SAMPLES=$N_SAMPLES")
println("=" ^ 65)

if Threads.nthreads() == 1
    @warn "Running with 1 Julia thread. Re-run with `julia -t auto` for multi-threaded FFT performance."
end

# ─────────────────────────────────────────────────────────────
# Scenario descriptor
# ─────────────────────────────────────────────────────────────

struct Scenario
    label::String
    dim::Int
    Nx::Int
    Ny::Int          # 1 for 1D
    Nz::Int          # 1 for 1D/2D
    t_end::Float64
    pml_size::Int
    medium_type::Symbol  # :homogeneous | :absorbing | :heterogeneous
end

# Convenience constructors
Scenario1D(label, Nx, t_end, mtype=:homogeneous) =
    Scenario(label, 1, Nx, 1, 1, t_end, 20, mtype)

Scenario2D(label, Nx, Ny, t_end, mtype=:homogeneous; pml=20) =
    Scenario(label, 2, Nx, Ny, 1, t_end, pml, mtype)

Scenario3D(label, Nx, Ny, Nz, t_end, mtype=:homogeneous; pml=12) =
    Scenario(label, 3, Nx, Ny, Nz, t_end, pml, mtype)

# Estimated total grid points (for display)
total_pts(s::Scenario) = s.Nx * max(1, s.Ny) * max(1, s.Nz)

# ─────────────────────────────────────────────────────────────
# Scenario list
# ─────────────────────────────────────────────────────────────
#
# t_end chosen so that:
#   dt ≈ CFL * dx / c0  =  0.3 * 1e-4 / 1500  ≈ 2e-8 s
#   Nt = t_end / dt
#
#   t_end=4e-6 → Nt≈200   t_end=2e-6 → Nt≈100   t_end=1e-6 → Nt≈50

SCENARIOS = [
    # ── 1D ──────────────────────────────────────────────────
    Scenario1D("1d_small",        256,    4e-6),
    Scenario1D("1d_medium",      1024,    4e-6),
    Scenario1D("1d_large",       4096,    4e-6),
    Scenario1D("1d_absorbing",   1024,    4e-6, :absorbing),
    Scenario1D("1d_hetero",      1024,    4e-6, :heterogeneous),

    # ── 2D ──────────────────────────────────────────────────
    Scenario2D("2d_small",         64,   64,  4e-6),
    Scenario2D("2d_medium",       256,  256,  4e-6),
    Scenario2D("2d_large",        512,  512,  2e-6),
    Scenario2D("2d_absorbing",    256,  256,  4e-6, :absorbing),
    Scenario2D("2d_hetero",       256,  256,  4e-6, :heterogeneous),

    # ── 3D ──────────────────────────────────────────────────
    Scenario3D("3d_small",         32,  32,  32,  2e-6; pml=8),
    Scenario3D("3d_medium",        64,  64,  64,  2e-6; pml=12),
    Scenario3D("3d_large",        128, 128, 128,  1e-6; pml=20),
    Scenario3D("3d_absorbing",     64,  64,  64,  2e-6, :absorbing; pml=12),
    Scenario3D("3d_hetero",        64,  64,  64,  2e-6, :heterogeneous; pml=12),
]

# ─────────────────────────────────────────────────────────────
# Medium factory
# ─────────────────────────────────────────────────────────────

function make_medium(s::Scenario, rng::AbstractRNG)
    if s.medium_type == :homogeneous
        return KWaveMedium(sound_speed=C0, density=RHO0)

    elseif s.medium_type == :absorbing
        # 0.75 dB/(MHz^y cm) — typical soft tissue
        return KWaveMedium(sound_speed=C0, density=RHO0,
                           alpha_coeff=0.75, alpha_power=1.5)

    elseif s.medium_type == :heterogeneous
        dims = s.dim == 1 ? (s.Nx,) :
               s.dim == 2 ? (s.Nx, s.Ny) :
                            (s.Nx, s.Ny, s.Nz)
        # ±10% speed of sound perturbation, uniform density
        c_map   = C0   .* (1.0 .+ 0.10 .* randn(rng, dims...))
        rho_map = RHO0 .* ones(dims...)
        return KWaveMedium(sound_speed=c_map, density=rho_map)
    end
end

# ─────────────────────────────────────────────────────────────
# Simulation setup helpers (return (kgrid, medium, source, sensor))
# ─────────────────────────────────────────────────────────────

function build_sim(s::Scenario, rng::AbstractRNG)
    medium = make_medium(s, rng)

    if s.dim == 1
        kgrid = KWaveGrid(s.Nx, DX)
        make_time!(kgrid, C0; t_end=s.t_end)

        p0        = zeros(Float64, s.Nx)
        p0[s.Nx÷2-5:s.Nx÷2+5] .= 1.0
        source    = KWaveSource(p0=p0)

        smask     = falses(s.Nx)
        smask[s.Nx÷4]     = true
        smask[3*s.Nx÷4]   = true
        sensor    = KWaveSensor(mask=smask, record=[:p])

    elseif s.dim == 2
        kgrid = KWaveGrid(s.Nx, DX, s.Ny, DX)
        make_time!(kgrid, C0; t_end=s.t_end)

        p0     = Float64.(make_disc(s.Nx, s.Ny, s.Nx÷2, s.Ny÷2, 5))
        source = KWaveSource(p0=p0)

        smask  = falses(s.Nx, s.Ny)
        smask[s.Nx÷4,   s.Ny÷2] = true
        smask[3*s.Nx÷4, s.Ny÷2] = true
        sensor = KWaveSensor(mask=smask, record=[:p])

    else  # 3D
        kgrid = KWaveGrid(s.Nx, DX, s.Ny, DX, s.Nz, DX)
        make_time!(kgrid, C0; t_end=s.t_end)

        p0     = Float64.(make_ball(s.Nx, s.Ny, s.Nz,
                                    s.Nx÷2, s.Ny÷2, s.Nz÷2, 3))
        source = KWaveSource(p0=p0)

        smask  = falses(s.Nx, s.Ny, s.Nz)
        smask[s.Nx÷4,   s.Ny÷2, s.Nz÷2] = true
        smask[3*s.Nx÷4, s.Ny÷2, s.Nz÷2] = true
        sensor = KWaveSensor(mask=smask, record=[:p])
    end

    return kgrid, medium, source, sensor
end

# ─────────────────────────────────────────────────────────────
# Timing harness
# ─────────────────────────────────────────────────────────────

"""
Run `f()` once for JIT warmup, then `N_WARMUP-1` additional warmup passes,
then `N_SAMPLES` timed passes. Returns a NamedTuple of timing stats.
"""
function timed_runs(f; n_warmup=N_WARMUP, n_samples=N_SAMPLES)
    # ── first run (captures JIT compilation time) ──
    first_stats = @timed f()
    t_first      = first_stats.time
    alloc_first  = first_stats.bytes

    # ── additional warmup ──
    for _ in 2:n_warmup
        f()
    end

    # ── timed samples ──
    times     = Vector{Float64}(undef, n_samples)
    allocs    = Vector{Int}(undef, n_samples)
    gc_times  = Vector{Float64}(undef, n_samples)

    for i in 1:n_samples
        s = @timed f()
        times[i]    = s.time
        allocs[i]   = s.bytes
        gc_times[i] = s.gctime
    end

    return (
        first_run_s      = t_first,
        first_run_alloc  = alloc_first,
        min_s            = minimum(times),
        median_s         = median(times),
        mean_s           = mean(times),
        max_s            = maximum(times),
        std_s            = std(times),
        median_alloc     = Int(round(median(allocs))),
        median_gc_s      = median(gc_times),
    )
end

# ─────────────────────────────────────────────────────────────
# Run full-simulation benchmarks
# ─────────────────────────────────────────────────────────────

function run_sim_benchmarks(; data_cast::Type=Float64, backend::String="cpu_f64")
    println("\n─── Full-simulation benchmarks  backend=$(backend) ───")

    rows = []

    for s in SCENARIOS
        rng = MersenneTwister(RNG_SEED)
        kgrid, medium, source, sensor = build_sim(s, rng)
        Nt = kgrid.Nt[]

        f = () -> kspace_first_order(kgrid, medium, source, sensor;
                                      smooth_p0=false,
                                      plot_sim=false,
                                      show_progress=false,
                                      pml_size=s.pml_size,
                                      data_cast=data_cast)

        @printf("  %-22s  Nx=%-5d Nt=%-5d  ", s.label, s.Nx, Nt)
        flush(stdout)

        try
            r = timed_runs(f)
            @printf("first=%6.3fs  median=%6.3fs  alloc=%5.1f MB\n",
                    r.first_run_s, r.median_s, r.median_alloc / 1e6)

            push!(rows, (
                timestamp   = string(now()),
                language    = "julia",
                backend     = backend,
                scenario    = s.label,
                dim         = s.dim,
                Nx          = s.Nx,
                Ny          = s.Ny,
                Nz          = s.Nz,
                Nt          = Nt,
                medium_type = string(s.medium_type),
                r...
            ))
        catch e
            println(" ERROR: $e")
        end
    end

    return rows
end

# ─────────────────────────────────────────────────────────────
# FFT component micro-benchmarks
# ─────────────────────────────────────────────────────────────

function run_fft_benchmarks()
    println("\n─── FFT component micro-benchmarks ───")
    println("    (using mul! with pre-allocated output — allocation-free timing)")

    rows = []

    # 1D rfft (matches solver: real → complex half-spectrum)
    # Single-threaded FFTW — small 1D sizes don't benefit from parallelism,
    # and 1 thread is already faster than MATLAB on these sizes.
    for N in [128, 256, 512, 1024, 2048, 4096, 8192, 16384]
        x     = rand(Float64, N)
        out   = zeros(ComplexF64, N ÷ 2 + 1)
        flops = 5 * N * log2(N)

        # ── FFTW plan (single-threaded for small N) ──
        FFTW.set_num_threads(1)
        plan_fftw = FFTW.plan_rfft(x; flags=FFTW.MEASURE)
        b_fftw    = @benchmark mul!($out, $plan_fftw, $x)  samples=1000  evals=5
        med_fftw  = median(b_fftw).time
        @printf("  rfft 1D FFTW  N=%-6d  median=%7.2f µs  %.2f GFLOP/s\n",
                N, med_fftw/1e3, flops/med_fftw)
        push!(rows, (language="julia", backend="cpu_f64", label="fft_1d",
                     dim=1, N=N, Nx=N, Ny=1, Nz=1,
                     median_ns=med_fftw, min_ns=minimum(b_fftw).time,
                     throughput_gflops=flops/med_fftw))

    end
    # Restore thread count for subsequent (larger) benchmarks
    FFTW.set_num_threads(Threads.nthreads())

    # 2D rfft
    for (Nx, Ny) in [(64,64),(128,128),(256,256),(512,512),(1024,1024)]
        x    = rand(Float64, Nx, Ny)
        plan = plan_rfft(x; flags=FFTW.MEASURE)
        out  = zeros(ComplexF64, Nx ÷ 2 + 1, Ny)
        b    = @benchmark mul!($out, $plan, $x)  samples=200  evals=3
        N    = Nx * Ny
        flops = 5 * N * log2(N)
        med_ns = median(b).time
        @printf("  rfft 2D  %4dx%-4d  median=%7.2f ms  %.2f GFLOP/s\n",
                Nx, Ny, med_ns/1e6, flops/med_ns)
        push!(rows, (language="julia", backend="cpu_f64", label="fft_2d",
                     dim=2, N=N, Nx=Nx, Ny=Ny, Nz=1,
                     median_ns=med_ns, min_ns=minimum(b).time,
                     throughput_gflops=flops/med_ns))
    end

    # 3D rfft
    for (Nx, Ny, Nz) in [(32,32,32),(64,64,64),(128,128,128),(256,256,256)]
        x    = rand(Float64, Nx, Ny, Nz)
        plan = plan_rfft(x; flags=FFTW.MEASURE)
        out  = zeros(ComplexF64, Nx ÷ 2 + 1, Ny, Nz)
        b    = @benchmark mul!($out, $plan, $x)  samples=50  evals=3
        N    = Nx * Ny * Nz
        flops = 5 * N * log2(N)
        med_ns = median(b).time
        @printf("  rfft 3D  %dx%dx%-3d  median=%7.2f ms  %.2f GFLOP/s\n",
                Nx, Ny, Nz, med_ns/1e6, flops/med_ns)
        push!(rows, (language="julia", backend="cpu_f64", label="fft_3d",
                     dim=3, N=N, Nx=Nx, Ny=Ny, Nz=Nz,
                     median_ns=med_ns, min_ns=minimum(b).time,
                     throughput_gflops=flops/med_ns))
    end

    return rows
end

# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

# ── CPU Float64 (default) ──────────────────────────────────
sim_rows_f64 = run_sim_benchmarks(data_cast=Float64, backend="cpu_f64")

# ── CPU Float32 (homogeneous only — matches MATLAB 'single' cast) ──
homog_scenarios_idx = findall(s -> s.medium_type == :homogeneous, SCENARIOS)
f32_rows = let rows = []
    println("\n─── Full-simulation benchmarks  backend=cpu_f32 ───")
    for idx in homog_scenarios_idx
        s = SCENARIOS[idx]
        rng = MersenneTwister(RNG_SEED)
        kgrid, medium, source, sensor = build_sim(s, rng)
        Nt = kgrid.Nt[]

        f = () -> kspace_first_order(kgrid, medium, source, sensor;
                                      smooth_p0=false,
                                      plot_sim=false,
                                      show_progress=false,
                                      pml_size=s.pml_size,
                                      data_cast=Float32)

        @printf("  %-22s  Nx=%-5d Nt=%-5d  ", s.label * "_f32", s.Nx, Nt)
        flush(stdout)
        try
            r = timed_runs(f)
            @printf("first=%6.3fs  median=%6.3fs  alloc=%5.1f MB\n",
                    r.first_run_s, r.median_s, r.median_alloc / 1e6)
            push!(rows, (
                timestamp   = string(now()),
                language    = "julia",
                backend     = "cpu_f32",
                scenario    = s.label,
                dim         = s.dim,
                Nx          = s.Nx,
                Ny          = s.Ny,
                Nz          = s.Nz,
                Nt          = Nt,
                medium_type = string(s.medium_type),
                r...
            ))
        catch e
            println(" ERROR: $e")
        end
    end
    rows
end

# ── FFT micro-benchmarks ──────────────────────────────────
fft_rows = run_fft_benchmarks()

# ─────────────────────────────────────────────────────────────
# Save results
# ─────────────────────────────────────────────────────────────

sim_df = DataFrame(vcat(sim_rows_f64, f32_rows))
fft_df = DataFrame(fft_rows)

sim_path = joinpath(RESULTS_DIR, "full_simulations.csv")
fft_path = joinpath(RESULTS_DIR, "fft_components.csv")
CSV.write(sim_path, sim_df)
CSV.write(fft_path, fft_df)

println("\n─── Saved ───")
println("  $sim_path")
println("  $fft_path")
println("\nDone.  $(now())")
