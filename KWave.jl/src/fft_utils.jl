# ============================================================================
# KWave.jl — FFT utilities and spectral gradient computation
# ============================================================================

"""
    FFTPlans{F,I}

Pre-computed real-to-complex FFT plans for efficient repeated transforms.

Uses `plan_rfft` / `plan_irfft` (real ↔ complex, exploiting conjugate symmetry):
- Forward plan maps a real array of shape `real_dims` to a complex array of
  shape `rfft_dims`, where `rfft_dims[1] = real_dims[1] ÷ 2 + 1`.
- Inverse plan maps the complex array back to the real array.

This halves the memory and compute cost vs a complex-to-complex FFT on real
data. For an (Nx, Ny) grid the complex scratch is only (Nx÷2+1, Ny).
"""
struct FFTPlans{F, I}
    forward::F           # plan_rfft:  real (Nx,…) → complex (Nx÷2+1,…)
    inverse::I           # plan_irfft: complex (Nx÷2+1,…) → real (Nx,…)
    real_dims::Dims      # original real dimensions, e.g. (Nx, Ny, Nz)
    rfft_dims::Dims      # rfft output dimensions, e.g. (Nx÷2+1, Ny, Nz)
end

"""Return the first-dimension size of the rfft output (= Nx÷2+1)."""
rfft_first_dim(p::FFTPlans) = p.rfft_dims[1]

"""
    create_fft_plans(dims; data_cast=Float64)

Create out-of-place rfft / irfft plans for arrays of the given dimensions.

Uses `FFTW.MEASURE` for good plans while keeping setup time reasonable.
Threads are set to `Threads.nthreads()` — launch Julia with `-t auto` (or
`-t N`) to exploit multi-core FFT parallelism.

# Arguments
- `dims`: Tuple of grid dimensions, e.g., `(Nx, Ny)`
- `data_cast`: Floating-point element type (default: `Float64`)

# Returns
An `FFTPlans` struct with forward (rfft) and inverse (irfft) plans.
"""
function create_fft_plans(dims::Tuple; data_cast::Type{T}=Float64) where T <: AbstractFloat
    # Use multi-threaded FFTW only for large arrays — for small grids the thread
    # synchronisation overhead exceeds the compute benefit and can make Julia
    # 10× slower than MATLAB on 1D scenarios.  Threshold chosen empirically:
    # single-threaded is faster below ~64 k points; multi-threaded above.
    nthreads = prod(dims) >= 65_536 ? Threads.nthreads() : 1
    FFTW.set_num_threads(nthreads)
    real_tmp    = zeros(T, dims...)
    rfft_dims   = (dims[1] ÷ 2 + 1, dims[2:end]...)
    complex_tmp = zeros(Complex{T}, rfft_dims...)
    fwd = plan_rfft(real_tmp;    flags=FFTW.MEASURE)
    inv = plan_irfft(complex_tmp, dims[1]; flags=FFTW.MEASURE)
    return FFTPlans(fwd, inv, dims, rfft_dims)
end

# ── Extension hook for 1D FFT backend ──────────────────────────────────────
#
# Extensions (e.g. KWaveAppleAccelerateExt) set this Ref in their __init__ to
# replace the default FFTW planner with a faster backend (e.g. vDSP).
# The factory signature is: (Nx::Int, ::Type{T}) -> FFTPlans
# Mutated only at runtime (__init__), never during precompilation — this avoids
# the "method overwriting during precompilation" error that arises from
# specializing create_fft_plans(::KWaveGrid1D) in both the base module and an
# extension.
const _CREATE_1D_PLANS = Ref{Any}(nothing)

# Grid-dispatching convenience wrappers

function create_fft_plans(grid::KWaveGrid1D; data_cast::Type{T}=Float64) where T
    factory = _CREATE_1D_PLANS[]
    factory !== nothing && return factory(grid.Nx, T)
    create_fft_plans((grid.Nx,); data_cast=T)
end

create_fft_plans(grid::KWaveGrid2D; data_cast::Type{T}=Float64) where T =
    create_fft_plans((grid.Nx, grid.Ny); data_cast=T)

create_fft_plans(grid::KWaveGrid3D; data_cast::Type{T}=Float64) where T =
    create_fft_plans((grid.Nx, grid.Ny, grid.Nz); data_cast=T)

# ============================================================================
# Spectral gradient computation
# ============================================================================

"""
    spectral_gradient!(df, f, k_vec, shift_op, scratch, plans, dim, ndims_total)

Compute the spatial gradient of the real field `f` along dimension `dim`
using the real-to-complex FFT (rfft):

    df/dx = irfft( im * kx_r .* shift_r .* rfft(f) )

where `kx_r` and `shift_r` are the wavenumber / shift vectors truncated to
the rfft output length along the first dimension (only for `dim == 1`).

# Arguments
- `df`          : Output real array (overwritten). Same shape as `f`.
- `f`           : Input real field array.
- `k_vec`       : Wavenumber vector for this dimension (any float element type).
                  For `dim == 1` only the first `rfft_first_dim(plans)` elements
                  are used (positive frequencies); for other dims all are used.
- `shift_op`    : Phase-shift operator for staggered grid (any complex element type).
                  Same truncation rule as `k_vec`.
- `scratch`     : Pre-allocated complex work array of shape `plans.rfft_dims`.
                  Element type determines the working precision T.
                  **Destroyed** by the inverse irfft — do not read after return.
- `plans`       : `FFTPlans` struct.
- `dim`         : Spatial dimension to differentiate (1, 2, or 3).
- `ndims_total` : Total number of spatial dimensions.

# Notes
- `k_vec` and `shift_op` are auto-converted to the precision `T` of `scratch`
  when they are Float64 and T is Float32. For the default Float64 case the
  grid vectors are used directly (zero allocation).
- The irfft destroys `scratch`; the caller must not use `scratch` after the
  call without first re-populating it.
"""
function spectral_gradient!(df::AbstractArray{T},
                            f::AbstractArray{T},
                            k_vec::AbstractVector,
                            shift_op::AbstractVector,
                            scratch::AbstractArray{Complex{T}},
                            plans::FFTPlans,
                            dim::Int,
                            ndims_total::Int) where T <: AbstractFloat

    # Forward rfft (out-of-place): f (real, full dims) → scratch (complex, rfft dims)
    mul!(scratch, plans.forward, f)

    # Ensure k and shift are T-typed for the broadcast.
    # For the default Float64 case the grid vectors already match T, so no allocation.
    k = k_vec  isa AbstractVector{T}         ? k_vec   : convert(Vector{T},          k_vec)
    s = shift_op isa AbstractVector{Complex{T}} ? shift_op : convert(Vector{Complex{T}}, shift_op)

    # The rfft output has Nx÷2+1 elements along the first dimension.
    # For dim == 1, only those positive-frequency elements of k / s are relevant;
    # for other dims, the full vector is used.
    n1 = size(scratch, 1)   # = Nx÷2+1

    if ndims_total == 1
        kv = reshape(@view(k[1:n1]), n1)
        sv = reshape(@view(s[1:n1]), n1)
        @. scratch = im * kv * sv * scratch

    elseif ndims_total == 2
        if dim == 1
            kv = reshape(@view(k[1:n1]), n1, 1)
            sv = reshape(@view(s[1:n1]), n1, 1)
        else  # dim == 2
            kv = reshape(k, 1, :)
            sv = reshape(s, 1, :)
        end
        @. scratch = im * kv * sv * scratch

    else  # ndims_total == 3
        if dim == 1
            kv = reshape(@view(k[1:n1]), n1, 1, 1)
            sv = reshape(@view(s[1:n1]), n1, 1, 1)
        elseif dim == 2
            kv = reshape(k, 1, :, 1)
            sv = reshape(s, 1, :, 1)
        else  # dim == 3
            kv = reshape(k, 1, 1, :)
            sv = reshape(s, 1, 1, :)
        end
        @. scratch = im * kv * sv * scratch
    end

    # Inverse irfft (out-of-place): scratch (complex) → df (real, full dims).
    # Note: irfft destroys the input (scratch) — this is expected and acceptable.
    mul!(df, plans.inverse, scratch)

    return df
end

# ============================================================================
# Pre-cast spectral operators
# ============================================================================

"""
    SpectralOps{T}

Wavenumber vectors and staggered-grid phase-shift operators pre-cast to the
working precision `T`.  Created once before the time loop so that per-step
calls to `spectral_gradient!` never trigger a Float64→T conversion (and the
associated heap allocation) when `T = Float32`.

For 1D grids `ky`, `kz` and the `ddy_*`/`ddz_*` fields are empty vectors.
For 2D grids `kz` and the `ddz_*` fields are empty.
"""
struct SpectralOps{T<:AbstractFloat}
    kx::Vector{T}
    ky::Vector{T}
    kz::Vector{T}
    ddx_pos::Vector{Complex{T}}
    ddx_neg::Vector{Complex{T}}
    ddy_pos::Vector{Complex{T}}
    ddy_neg::Vector{Complex{T}}
    ddz_pos::Vector{Complex{T}}
    ddz_neg::Vector{Complex{T}}
end

function create_spectral_ops(grid::KWaveGrid1D; data_cast::Type{T}=Float64) where T<:AbstractFloat
    SpectralOps{T}(
        T.(grid.kx_vec), T[], T[],
        Complex{T}.(grid.ddx_k_shift_pos), Complex{T}.(grid.ddx_k_shift_neg),
        Complex{T}[], Complex{T}[], Complex{T}[], Complex{T}[],
    )
end

function create_spectral_ops(grid::KWaveGrid2D; data_cast::Type{T}=Float64) where T<:AbstractFloat
    SpectralOps{T}(
        T.(grid.kx_vec), T.(grid.ky_vec), T[],
        Complex{T}.(grid.ddx_k_shift_pos), Complex{T}.(grid.ddx_k_shift_neg),
        Complex{T}.(grid.ddy_k_shift_pos), Complex{T}.(grid.ddy_k_shift_neg),
        Complex{T}[], Complex{T}[],
    )
end

function create_spectral_ops(grid::KWaveGrid3D; data_cast::Type{T}=Float64) where T<:AbstractFloat
    SpectralOps{T}(
        T.(grid.kx_vec), T.(grid.ky_vec), T.(grid.kz_vec),
        Complex{T}.(grid.ddx_k_shift_pos), Complex{T}.(grid.ddx_k_shift_neg),
        Complex{T}.(grid.ddy_k_shift_pos), Complex{T}.(grid.ddy_k_shift_neg),
        Complex{T}.(grid.ddz_k_shift_pos), Complex{T}.(grid.ddz_k_shift_neg),
    )
end
