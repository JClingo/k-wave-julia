# ============================================================================
# KWave.jl — Perfectly Matched Layer (PML) implementation
# ============================================================================

"""
    get_pml(N, d, pml_size, pml_alpha, dt; staggered=false)

Generate a 1D PML absorption profile for a grid of `N` points.

The PML uses a 4th-order polynomial absorption profile that decays from the
domain edges inward. The profile is symmetric — applied to both ends of the axis.

# Arguments
- `N`: Number of grid points along this axis
- `d`: Grid spacing [m]
- `pml_size`: Number of PML grid points on each side
- `pml_alpha`: PML absorption coefficient
- `dt`: Time step [s]
- `staggered`: If true, shift profile by half a grid point for staggered grids

# Returns
A `Vector{Float64}` of length `N` with values in (0, 1].
Interior points have value 1.0, PML regions have values < 1.0.
"""
function get_pml(N::Int, d::Float64, pml_size::Int, pml_alpha::Float64, dt::Float64;
                 staggered::Bool=false)
    # Initialize PML to ones (no absorption)
    pml = ones(Float64, N)

    if pml_size == 0
        return pml
    end

    # Grid point indices within PML (1-based, measured from edge)
    # Left PML: points 1 to pml_size
    # Right PML: points (N - pml_size + 1) to N
    for i in 1:pml_size
        # Normalized distance from the inner PML boundary (0 at inner edge, 1 at outer edge)
        # i=1 is the outermost grid point on the left; i=pml_size is innermost on the left
        if staggered
            pos = (pml_size - i + 0.5) / pml_size
        else
            pos = (pml_size - i + 1) / pml_size
        end

        # 4th-order polynomial absorption profile
        absorption = pml_alpha * pos^4
        pml_val = exp(-absorption * dt)

        # Apply symmetrically to both ends
        pml[i] = pml_val
        pml[N - i + 1] = pml_val
    end

    return pml
end

"""
    get_optimal_pml_size(grid_size; pml_range=10:60, axisymmetric=false)

Find the PML size that results in the most efficient FFT grid dimensions.

The total grid size (including PML) should have small prime factors
(ideally powers of 2, 3, 5) for optimal FFT performance.

# Arguments
- `grid_size`: Tuple of grid dimensions (Nx,) or (Nx, Ny) or (Nx, Ny, Nz)
- `pml_range`: Range of PML sizes to search
- `axisymmetric`: If true, only add PML to one side of the first dimension

# Returns
A tuple of optimal PML sizes, one per dimension.
"""
function get_optimal_pml_size(grid_size::Tuple; pml_range::AbstractRange=10:60)
    return Tuple(
        _optimal_pml_for_dim(n, pml_range) for n in grid_size
    )
end

function get_optimal_pml_size(n::Int; pml_range::AbstractRange=10:60)
    return _optimal_pml_for_dim(n, pml_range)
end

"""
Score a number by its largest prime factor (smaller is better for FFT).
"""
function _max_prime_factor(n::Int)
    if n <= 1
        return 1
    end
    max_factor = 1
    d = 2
    m = n
    while d * d <= m
        while m % d == 0
            max_factor = max(max_factor, d)
            m ÷= d
        end
        d += 1
    end
    if m > 1
        max_factor = max(max_factor, m)
    end
    return max_factor
end

function _optimal_pml_for_dim(N::Int, pml_range::AbstractRange)
    best_pml = first(pml_range)
    best_score = typemax(Int)

    for pml_size in pml_range
        total = N + 2 * pml_size
        score = _max_prime_factor(total)
        if score < best_score
            best_score = score
            best_pml = pml_size
        end
    end

    return best_pml
end
