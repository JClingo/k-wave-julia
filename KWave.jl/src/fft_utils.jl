# ============================================================================
# KWave.jl — FFT utilities and spectral gradient computation
# ============================================================================

"""
    FFTPlans

Pre-computed FFT plans for efficient repeated transforms.
"""
struct FFTPlans
    forward::Any   # FFTW forward plan (complex in-place)
    inverse::Any   # FFTW inverse plan (complex in-place)
end

"""
    create_fft_plans(dims; data_cast=Float64)

Create in-place FFT/IFFT plans for arrays of the given dimensions.

Uses `FFTW.MEASURE` flag for optimal performance (plans are cached internally by FFTW).

# Arguments
- `dims`: Tuple of grid dimensions, e.g., `(Nx, Ny)`
- `data_cast`: Floating point type (default: Float64)

# Returns
An `FFTPlans` struct with forward and inverse plans.
"""
function create_fft_plans(dims::Tuple; data_cast::Type{T}=Float64) where T<:AbstractFloat
    tmp = zeros(Complex{T}, dims...)
    fwd = plan_fft!(tmp; flags=FFTW.MEASURE)
    inv = plan_ifft!(tmp; flags=FFTW.MEASURE)
    return FFTPlans(fwd, inv)
end

"""
    create_fft_plans(grid::KWaveGrid2D; data_cast=Float64)

Create FFT plans sized for the given grid.
"""
create_fft_plans(grid::KWaveGrid2D; data_cast::Type{T}=Float64) where T =
    create_fft_plans((grid.Nx, grid.Ny); data_cast=data_cast)

create_fft_plans(grid::KWaveGrid1D; data_cast::Type{T}=Float64) where T =
    create_fft_plans((grid.Nx,); data_cast=data_cast)

create_fft_plans(grid::KWaveGrid3D; data_cast::Type{T}=Float64) where T =
    create_fft_plans((grid.Nx, grid.Ny, grid.Nz); data_cast=data_cast)

# ============================================================================
# Spectral gradient computation
# ============================================================================

"""
    spectral_gradient!(df, f, k_vec, shift_op, scratch, plans, dim, ndims_total)

Compute the spatial gradient of `f` along dimension `dim` using the FFT.

    df/dx = real(IFFT(im * kx .* shift_op .* FFT(f)))

This computes the gradient on a staggered grid when `shift_op` contains
the appropriate phase shift factors.

# Arguments
- `df`: Output array (overwritten with gradient values)
- `f`: Input field array (real-valued)
- `k_vec`: 1D wavenumber vector for this dimension
- `shift_op`: 1D complex shift operator for staggered grid
- `scratch`: Pre-allocated complex scratch array (same size as `f`)
- `plans`: FFTPlans struct
- `dim`: Dimension along which to differentiate (1, 2, or 3)
- `ndims_total`: Total number of dimensions (1, 2, or 3)
"""
function spectral_gradient!(df::AbstractArray, f::AbstractArray,
                            k_vec::Vector{Float64}, shift_op::Vector{ComplexF64},
                            scratch::AbstractArray{<:Complex}, plans::FFTPlans,
                            dim::Int, ndims_total::Int)
    # Copy real data into complex scratch
    scratch .= complex.(f)

    # Forward FFT (in-place)
    plans.forward * scratch

    # Multiply by im * k * shift in frequency domain
    # k_vec and shift_op are 1D — reshape for broadcasting along the correct dimension
    if ndims_total == 1
        @. scratch = im * k_vec * shift_op * scratch
    elseif ndims_total == 2
        if dim == 1
            # kx is along dim 1: shape (Nx, 1)
            @. scratch = im * k_vec * shift_op * scratch  # broadcasts kx along columns
        else
            # ky is along dim 2: shape (1, Ny)
            k_row = reshape(k_vec, 1, :)
            s_row = reshape(shift_op, 1, :)
            @. scratch = im * k_row * s_row * scratch
        end
    elseif ndims_total == 3
        if dim == 1
            @. scratch = im * k_vec * shift_op * scratch
        elseif dim == 2
            k_r = reshape(k_vec, 1, :, 1)
            s_r = reshape(shift_op, 1, :, 1)
            @. scratch = im * k_r * s_r * scratch
        else
            k_r = reshape(k_vec, 1, 1, :)
            s_r = reshape(shift_op, 1, 1, :)
            @. scratch = im * k_r * s_r * scratch
        end
    end

    # Inverse FFT (in-place)
    plans.inverse * scratch

    # Extract real part
    @. df = real(scratch)

    return df
end
