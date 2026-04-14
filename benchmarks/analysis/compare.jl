"""
KWave Benchmark Comparison — Julia vs MATLAB
=============================================
Loads CSV results from both runtimes, prints a formatted speedup table,
and generates plots.

Usage:
    julia --project=benchmarks/julia benchmarks/analysis/compare.jl

Expected inputs  (all optional — only present CSVs are loaded):
    benchmarks/results/julia/full_simulations.csv
    benchmarks/results/matlab/full_simulations.csv
    benchmarks/results/julia/fft_components.csv
    benchmarks/results/matlab/fft_components.csv
    benchmarks/results/julia/gpu_simulations.csv
    benchmarks/results/matlab/gpu_simulations.csv  (if available)

Outputs:
    benchmarks/results/plots/speedup_full.png
    benchmarks/results/plots/scaling_cpu.png
    benchmarks/results/plots/fft_throughput.png
    benchmarks/results/plots/gpu_speedup.png
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "julia"))

using CSV
using DataFrames
using Statistics
using Printf
using CairoMakie

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
const PLOTS_DIR   = joinpath(RESULTS_DIR, "plots")
mkpath(PLOTS_DIR)

# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

function load_csv(path)
    isfile(path) ? CSV.read(path, DataFrame) : nothing
end

function section(title)
    println()
    println("=" ^ 65)
    println("  $title")
    println("=" ^ 65)
end

# ─────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────

julia_sim  = load_csv(joinpath(RESULTS_DIR, "julia",  "full_simulations.csv"))
matlab_sim = load_csv(joinpath(RESULTS_DIR, "matlab", "full_simulations.csv"))
julia_fft  = load_csv(joinpath(RESULTS_DIR, "julia",  "fft_components.csv"))
matlab_fft = load_csv(joinpath(RESULTS_DIR, "matlab", "fft_components.csv"))
julia_gpu  = load_csv(joinpath(RESULTS_DIR, "julia",  "gpu_simulations.csv"))
matlab_gpu = load_csv(joinpath(RESULTS_DIR, "matlab", "gpu_simulations.csv"))

# ─────────────────────────────────────────────────────────────
# Full-simulation speedup table  (cpu_f64 comparison)
# ─────────────────────────────────────────────────────────────

if !isnothing(julia_sim) && !isnothing(matlab_sim)
    section("Full Simulation Speedup  (cpu_f64, median time)")

    jdf = filter(r -> r.backend == "cpu_f64", julia_sim)
    mdf = filter(r -> r.backend == "cpu_f64", matlab_sim)

    # Build lookup: scenario → median time
    j_times = Dict(r.scenario => r.median_s for r in eachrow(jdf))
    m_times = Dict(r.scenario => r.median_s for r in eachrow(mdf))

    scenarios = sort(collect(intersect(keys(j_times), keys(m_times))))

    @printf("  %-22s  %9s  %9s  %8s  %6s\n",
            "Scenario", "Julia(s)", "MATLAB(s)", "Speedup", "Dim")
    println("  " * "─"^60)

    speedups = Float64[]
    labels   = String[]
    dims     = Int[]

    for sc in scenarios
        jt = j_times[sc]
        mt = m_times[sc]
        sp = mt / jt   # >1 means Julia is faster
        dim_row = filter(r -> r.scenario == sc, jdf)
        d = isempty(dim_row) ? 0 : first(dim_row).dim

        @printf("  %-22s  %9.3f  %9.3f  %7.2fx  %5dD\n",
                sc, jt, mt, sp, d)
        push!(speedups, sp)
        push!(labels, sc)
        push!(dims, d)
    end

    println("  " * "─"^60)
    @printf("  %-22s  %9s  %9s  %7.2fx\n",
            "GEOMETRIC MEAN", "", "", exp(mean(log.(speedups))))

    # ── Plot: speedup by scenario ────────────────────────────────────────
    fig = Figure(size=(900, 500))
    ax  = Axis(fig[1,1],
               title    = "Julia vs MATLAB Speedup — Full Simulations (cpu_f64)",
               xlabel   = "Scenario",
               ylabel   = "Speedup (MATLAB time / Julia time)",
               xticklabelrotation = π/3,
               xticklabelsize = 11)

    colors = map(d -> d==1 ? :steelblue : d==2 ? :darkorange : :forestgreen, dims)
    xs     = 1:length(labels)
    barplot!(ax, xs, speedups; color=colors, bar_labels=nothing)
    hlines!(ax, [1.0]; color=:red, linestyle=:dash, linewidth=1.5,
            label="parity")
    ax.xticks = (xs, labels)

    # Legend for dimension colours
    elem1 = PolyElement(color=:steelblue)
    elem2 = PolyElement(color=:darkorange)
    elem3 = PolyElement(color=:forestgreen)
    Legend(fig[1,2], [elem1, elem2, elem3], ["1D", "2D", "3D"],
           "Dimension"; framevisible=false)

    save(joinpath(PLOTS_DIR, "speedup_full.png"), fig; px_per_unit=2)
    println("\n  Saved: plots/speedup_full.png")

    # ── Plot: Float32 speedup comparison ───────────────────────────────
    jf32 = filter(r -> r.backend == "cpu_f32", julia_sim)
    mf32 = filter(r -> r.backend == "cpu_f32", matlab_sim)
    if !isempty(jf32) && !isempty(mf32)
        j32_times = Dict(r.scenario => r.median_s for r in eachrow(jf32))
        m32_times = Dict(r.scenario => r.median_s for r in eachrow(mf32))
        sc32 = sort(collect(intersect(keys(j32_times), keys(m32_times))))

        section("Full Simulation Speedup  (cpu_f32, median time)")
        @printf("  %-22s  %9s  %9s  %8s\n",
                "Scenario", "Julia_f32", "MATLAB_f32", "Speedup")
        println("  " * "─"^55)
        sp32 = Float64[]
        for sc in sc32
            sp = m32_times[sc] / j32_times[sc]
            @printf("  %-22s  %9.3f  %9.3f  %7.2fx\n",
                    sc, j32_times[sc], m32_times[sc], sp)
            push!(sp32, sp)
        end
        println("  " * "─"^55)
        @printf("  %-22s  %9s  %9s  %7.2fx\n",
                "GEOMETRIC MEAN", "", "", exp(mean(log.(sp32))))
    end
end

# ─────────────────────────────────────────────────────────────
# Scaling plot  (time vs N grid points, log–log)
# ─────────────────────────────────────────────────────────────

if !isnothing(julia_sim)
    section("Scaling Analysis  (Julia cpu_f64, time vs total grid points)")

    jdf = filter(r -> r.backend == "cpu_f64", julia_sim)
    jdf.total_pts = jdf.Nx .* max.(1, jdf.Ny) .* max.(1, jdf.Nz)
    # normalise by Nt so we get time-per-step
    jdf.time_per_step = jdf.median_s ./ jdf.Nt

    fig = Figure(size=(800, 500))
    ax  = Axis(fig[1,1],
               title  = "Julia Scaling: time/step vs grid points (cpu_f64, log–log)",
               xlabel = "Total grid points",
               ylabel = "Wall time per step (s)",
               xscale = log10,
               yscale = log10)

    dim_colors = Dict(1 => :steelblue, 2 => :darkorange, 3 => :forestgreen)
    dim_markers = Dict(1 => :circle, 2 => :rect, 3 => :diamond)

    for d in [1, 2, 3]
        sub = filter(r -> r.dim == d && r.medium_type == "homogeneous", jdf)
        isempty(sub) && continue
        sort!(sub, :total_pts)
        scatter!(ax, sub.total_pts, sub.time_per_step;
                 color=dim_colors[d], marker=dim_markers[d],
                 markersize=12, label="$(d)D")
        lines!(ax, sub.total_pts, sub.time_per_step;
               color=dim_colors[d], linewidth=1.5)
    end

    # Reference slope: O(N log N) guideline
    pts = exp10.(range(2, 7, length=100))
    t_ref = pts .* log.(pts) .* 1e-10   # arbitrary scaling
    lines!(ax, pts, t_ref; color=:grey, linestyle=:dash,
           linewidth=1, label="O(N log N)")
    axislegend(ax; position=:lt)

    save(joinpath(PLOTS_DIR, "scaling_cpu.png"), fig; px_per_unit=2)
    println("  Saved: plots/scaling_cpu.png")
end

# ─────────────────────────────────────────────────────────────
# FFT throughput comparison
# ─────────────────────────────────────────────────────────────

if !isnothing(julia_fft) && !isnothing(matlab_fft)
    section("FFT Throughput Comparison (GFLOP/s)")

    all_fft = vcat(julia_fft, matlab_fft)

    @printf("  %-20s  %8s  %8s  %8s\n", "FFT size", "Julia", "MATLAB", "Speedup")
    println("  " * "─"^50)

    for dim in [1, 2, 3]
        lbl = dim == 1 ? "fft_1d" : dim == 2 ? "fft_2d" : "fft_3d"
        jrows = filter(r -> r.language == "julia"  && r.label == lbl, all_fft)
        mrows = filter(r -> r.language == "matlab" && r.label == lbl, all_fft)
        isempty(jrows) || isempty(mrows) && continue
        sort!(jrows, :N); sort!(mrows, :N)

        j_dict = Dict(r.N => r.throughput_gflops for r in eachrow(jrows))
        m_dict = Dict(r.N => r.throughput_gflops for r in eachrow(mrows))
        Ns = sort(collect(intersect(keys(j_dict), keys(m_dict))))

        println("  ── $(dim)D ──")
        for N in Ns
            sp = j_dict[N] / m_dict[N]
            @printf("  %-20s  %8.2f  %8.2f  %7.2fx\n",
                    "N=$(N)", j_dict[N], m_dict[N], sp)
        end
    end

    # Plot: 1D FFT throughput
    j1d = filter(r -> r.language == "julia"  && r.label == "fft_1d", all_fft)
    m1d = filter(r -> r.language == "matlab" && r.label == "fft_1d", all_fft)

    if !isempty(j1d) && !isempty(m1d)
        sort!(j1d, :N); sort!(m1d, :N)

        fig = Figure(size=(700, 400))
        ax  = Axis(fig[1,1],
                   title  = "1D FFT Throughput: Julia (FFTW) vs MATLAB",
                   xlabel = "FFT size N",
                   ylabel = "GFLOP/s",
                   xscale = log2)

        lines!(ax, j1d.N, j1d.throughput_gflops; color=:steelblue,
               linewidth=2, label="Julia (FFTW)")
        scatter!(ax, j1d.N, j1d.throughput_gflops; color=:steelblue, markersize=8)

        lines!(ax, m1d.N, m1d.throughput_gflops; color=:tomato,
               linewidth=2, label="MATLAB")
        scatter!(ax, m1d.N, m1d.throughput_gflops; color=:tomato, markersize=8)

        axislegend(ax; position=:lt)
        save(joinpath(PLOTS_DIR, "fft_throughput.png"), fig; px_per_unit=2)
        println("  Saved: plots/fft_throughput.png")
    end
end

# ─────────────────────────────────────────────────────────────
# GPU speedup  (Julia: Metal/CUDA vs cpu_f32)
# ─────────────────────────────────────────────────────────────

if !isnothing(julia_gpu)
    section("GPU Speedup Summary  (Julia)")

    gpu_backends = unique(filter(b -> contains(b, "f32") && b != "cpu_f32",
                                 julia_gpu.backend))

    for gb in gpu_backends
        println("\n  ── $gb vs cpu_f32 ──")
        @printf("  %-14s  %8s  %8s  %8s  %6s\n",
                "Scenario", "cpu_f32", gb, "Speedup", "Dim")
        println("  " * "─"^55)

        cpu_rows = filter(r -> r.backend == "cpu_f32", julia_gpu)
        gpu_rows = filter(r -> r.backend == gb, julia_gpu)

        c_dict = Dict(r.scenario => r.median_s for r in eachrow(cpu_rows))
        g_dict = Dict(r.scenario => r.median_s for r in eachrow(gpu_rows))
        dim_d  = Dict(r.scenario => r.dim      for r in eachrow(cpu_rows))

        sc_gpu = sort(collect(intersect(keys(c_dict), keys(g_dict))))
        sp_vals = Float64[]
        sc_lbls = String[]
        dims_g  = Int[]

        for sc in sc_gpu
            sp = c_dict[sc] / g_dict[sc]
            push!(sp_vals, sp); push!(sc_lbls, sc); push!(dims_g, get(dim_d, sc, 0))
            @printf("  %-14s  %8.3f  %8.3f  %7.2fx  %4dD\n",
                    sc, c_dict[sc], g_dict[sc], sp, dim_d[sc])
        end

        isempty(sp_vals) && continue

        # Plot
        fig = Figure(size=(700, 400))
        ax  = Axis(fig[1,1],
                   title             = "Julia GPU Speedup: $gb vs cpu_f32",
                   ylabel            = "Speedup",
                   xticklabelrotation = π/3, xticklabelsize=11)

        dim_colors = map(d -> d==2 ? :darkorange : :forestgreen, dims_g)
        xs = 1:length(sc_lbls)
        barplot!(ax, xs, sp_vals; color=dim_colors)
        hlines!(ax, [1.0]; color=:red, linestyle=:dash, linewidth=1.5)
        ax.xticks = (xs, sc_lbls)
        fname = "gpu_speedup_$(replace(gb, '/' => '_')).png"
        save(joinpath(PLOTS_DIR, fname), fig; px_per_unit=2)
        println("  Saved: plots/$fname")
    end
end

# ─────────────────────────────────────────────────────────────
# JIT compilation overhead
# ─────────────────────────────────────────────────────────────

if !isnothing(julia_sim)
    section("Julia JIT Compilation Overhead  (first_run vs median)")

    jdf = filter(r -> r.backend == "cpu_f64", julia_sim)
    @printf("  %-22s  %10s  %10s  %8s\n",
            "Scenario", "first_run(s)", "median(s)", "overhead")
    println("  " * "─"^57)
    for r in eachrow(jdf)
        overhead = r.first_run_s / r.median_s
        @printf("  %-22s  %10.3f  %10.3f  %7.2fx\n",
                r.scenario, r.first_run_s, r.median_s, overhead)
    end
end

println("\nDone.  All plots saved to: $PLOTS_DIR")
