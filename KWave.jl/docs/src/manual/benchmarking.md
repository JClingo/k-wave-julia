# Benchmarking

The `benchmarks/` directory contains scripts for measuring KWave.jl performance
and comparing it against the reference MATLAB k-Wave toolbox.

## Directory layout

```
benchmarks/
├── julia/
│   ├── Project.toml          # isolated env (pins KWave.jl from source)
│   ├── run_benchmarks.jl     # CPU benchmark suite
│   └── run_benchmarks_gpu.jl # GPU benchmark suite
├── matlab/
│   └── run_benchmarks.m      # MATLAB reference suite
├── analysis/
│   └── compare.jl            # Julia vs MATLAB comparison + plots
└── results/
    ├── julia/                 # CSV output from Julia runs
    ├── matlab/                # CSV output from MATLAB runs
    └── plots/                 # PNG plots from compare.jl
```

## Julia CPU benchmarks

Runs full simulations (1D / 2D / 3D) in Float64 and Float32, plus FFT
micro-benchmarks. Requires Julia ≥ 1.11.

```bash
julia --project=benchmarks/julia benchmarks/julia/run_benchmarks.jl
```

Override sample count (default 5):

```bash
BENCH_SAMPLES=10 julia --project=benchmarks/julia benchmarks/julia/run_benchmarks.jl
```

**Output:**
- `benchmarks/results/julia/full_simulations.csv`
- `benchmarks/results/julia/fft_components.csv`

## Julia GPU benchmarks

Benchmarks Metal (Apple Silicon) and CUDA backends against the CPU Float32
baseline.

```bash
# Metal (Mac)
julia --project=benchmarks/julia benchmarks/julia/run_benchmarks_gpu.jl metal

# CUDA (NVIDIA)
julia --project=benchmarks/julia benchmarks/julia/run_benchmarks_gpu.jl cuda

# Auto-detect
julia --project=benchmarks/julia benchmarks/julia/run_benchmarks_gpu.jl
```

**Output:** `benchmarks/results/julia/gpu_simulations.csv`

## MATLAB reference benchmarks

Requires MATLAB R2021a+ with the k-Wave toolbox on the path.

```matlab
% CPU only
run_benchmarks

% CPU + GPU (requires Parallel Computing Toolbox)
run_benchmarks('gpu')
```

**Output:**
- `benchmarks/results/matlab/full_simulations.csv`
- `benchmarks/results/matlab/fft_components.csv`
- `benchmarks/results/matlab/gpu_simulations.csv` (if GPU requested)

## Comparison and plots

After collecting results from both runtimes, run:

```bash
julia --project=benchmarks/julia benchmarks/analysis/compare.jl
```

Prints speedup tables to stdout and writes plots to `benchmarks/results/plots/`:

| Plot file | Contents |
|-----------|----------|
| `speedup_full.png` | Julia vs MATLAB speedup per scenario (cpu\_f64) |
| `scaling_cpu.png` | Time/step vs grid points, log–log (1D/2D/3D) |
| `fft_throughput.png` | 1D FFT GFLOP/s: Julia (FFTW) vs MATLAB |
| `gpu_speedup_*.png` | GPU backend vs cpu\_f32 speedup |

Comparison is graceful — only CSVs that exist are loaded. You can compare
Julia-only results without MATLAB data.

## Scenarios

All suites run the same 15 scenarios at matched grid sizes and time spans:

| Name | Dim | Grid | Medium |
|------|-----|------|--------|
| `1d_small` | 1D | 256 | homogeneous |
| `1d_medium` | 1D | 1024 | homogeneous |
| `1d_large` | 1D | 4096 | homogeneous |
| `1d_absorbing` | 1D | 1024 | absorbing |
| `1d_hetero` | 1D | 1024 | heterogeneous ±10% |
| `2d_small` | 2D | 64×64 | homogeneous |
| `2d_medium` | 2D | 256×256 | homogeneous |
| `2d_large` | 2D | 512×512 | homogeneous |
| `2d_absorbing` | 2D | 256×256 | absorbing |
| `2d_hetero` | 2D | 256×256 | heterogeneous ±10% |
| `3d_small` | 3D | 32³ | homogeneous |
| `3d_medium` | 3D | 64³ | homogeneous |
| `3d_large` | 3D | 128³ | homogeneous |
| `3d_absorbing` | 3D | 64³ | absorbing |
| `3d_hetero` | 3D | 64³ | heterogeneous ±10% |

Float32 runs cover homogeneous scenarios only (matches MATLAB `DataCast 'single'`).

## Timing methodology

Each scenario runs:
1. One first-run pass (captures JIT/planner overhead, recorded separately).
2. `N_WARMUP - 1` additional warmup passes (default: 0 extra).
3. `N_SAMPLES` timed passes (default: 5). Median is the primary metric.

Reported fields: `first_run_s`, `min_s`, `median_s`, `mean_s`, `max_s`,
`std_s`, `median_alloc` (bytes, Julia only), `median_gc_s` (Julia only).

## First-run vs steady-state

Julia's first-run time includes LLVM compilation and FFTW plan construction.
The comparison script reports the ratio `first_run / median` per scenario so
you can see JIT overhead separately from steady-state throughput.
