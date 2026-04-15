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
    FFTW.set_num_threads(Threads.nthreads())
    real_tmp    = zeros(T, dims...)
    rfft_dims   = (dims[1] ÷ 2 + 1, dims[2:end]...)
    complex_tmp = zeros(Complex{T}, rfft_dims...)
    fwd = plan_rfft(real_tmp;    flags=FFTW.MEASURE)
    inv = plan_irfft(complex_tmp, dims[1]; flags=FFTW.MEASURE)
    return FFTPlans(fwd, inv, dims, rfft_dims)
end

# Grid-dispatching convenience wrappers

create_fft_plans(grid::KWaveGrid1D; data_cast::Type{T}=Float64) where T =
    create_fft_plans((grid.Nx,); data_cast=T)

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
