# ============================================================================
# KWave.jl — KWaveGrid types and constructors
# ============================================================================

"""
    AbstractKWaveGrid

Abstract supertype for all k-Wave computational grids.
"""
abstract type AbstractKWaveGrid end

# --------------------------------------------------------------------------
# KWaveGrid1D
# --------------------------------------------------------------------------

"""
    KWaveGrid1D <: AbstractKWaveGrid

1D computational grid for k-space simulations.
"""
struct KWaveGrid1D <: AbstractKWaveGrid
    Nx::Int
    dx::Float64
    # Spatial coordinate vector
    x_vec::Vector{Float64}
    # Wavenumber vector
    kx_vec::Vector{Float64}
    k::Vector{Float64}
    # Staggered grid shift operators
    ddx_k_shift_pos::Vector{ComplexF64}  # exp(+i*kx*dx/2)
    ddx_k_shift_neg::Vector{ComplexF64}  # exp(-i*kx*dx/2)
    # Time parameters (mutable)
    dt::Base.RefValue{Float64}
    Nt::Base.RefValue{Int}
    t_array::Vector{Float64}
end

# --------------------------------------------------------------------------
# KWaveGrid2D
# --------------------------------------------------------------------------

"""
    KWaveGrid2D <: AbstractKWaveGrid

2D computational grid for k-space simulations.
"""
struct KWaveGrid2D <: AbstractKWaveGrid
    Nx::Int
    dx::Float64
    Ny::Int
    dy::Float64
    # Spatial coordinate vectors
    x_vec::Vector{Float64}
    y_vec::Vector{Float64}
    # Wavenumber vectors (1D, broadcast during use)
    kx_vec::Vector{Float64}
    ky_vec::Vector{Float64}
    # Wavenumber magnitude matrix
    k::Matrix{Float64}
    # Staggered grid shift operators (1D vectors)
    ddx_k_shift_pos::Vector{ComplexF64}
    ddx_k_shift_neg::Vector{ComplexF64}
    ddy_k_shift_pos::Vector{ComplexF64}
    ddy_k_shift_neg::Vector{ComplexF64}
    # Time parameters (mutable)
    dt::Base.RefValue{Float64}
    Nt::Base.RefValue{Int}
    t_array::Vector{Float64}
end

# --------------------------------------------------------------------------
# KWaveGrid3D
# --------------------------------------------------------------------------

"""
    KWaveGrid3D <: AbstractKWaveGrid

3D computational grid for k-space simulations.
"""
struct KWaveGrid3D <: AbstractKWaveGrid
    Nx::Int
    dx::Float64
    Ny::Int
    dy::Float64
    Nz::Int
    dz::Float64
    # Spatial coordinate vectors
    x_vec::Vector{Float64}
    y_vec::Vector{Float64}
    z_vec::Vector{Float64}
    # Wavenumber vectors (1D, broadcast during use)
    kx_vec::Vector{Float64}
    ky_vec::Vector{Float64}
    kz_vec::Vector{Float64}
    # Wavenumber magnitude 3D array
    k::Array{Float64, 3}
    # Staggered grid shift operators (1D vectors)
    ddx_k_shift_pos::Vector{ComplexF64}
    ddx_k_shift_neg::Vector{ComplexF64}
    ddy_k_shift_pos::Vector{ComplexF64}
    ddy_k_shift_neg::Vector{ComplexF64}
    ddz_k_shift_pos::Vector{ComplexF64}
    ddz_k_shift_neg::Vector{ComplexF64}
    # Time parameters (mutable)
    dt::Base.RefValue{Float64}
    Nt::Base.RefValue{Int}
    t_array::Vector{Float64}
end

# ============================================================================
# Wavenumber vector generation
# ============================================================================

"""
    _make_wavenumber_vec(N, d)

Generate the wavenumber vector for a grid of `N` points with spacing `d`.
Uses the MATLAB/FFTW convention: [0, 1, ..., N/2-1, -N/2, ..., -1] * (2π / (N*d))
"""
function _make_wavenumber_vec(N::Int, d::Float64)
    dk = 2π / (N * d)
    if iseven(N)
        # [0, 1, ..., N/2-1, -N/2, ..., -1]
        kv = [collect(0:N÷2-1); collect(-N÷2:-1)] .* dk
    else
        # [0, 1, ..., (N-1)/2, -(N-1)/2, ..., -1]
        kv = [collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)] .* dk
    end
    return kv
end

"""
    _make_spatial_vec(N, d)

Generate a centered spatial coordinate vector: -(N-1)/2*d to (N-1)/2*d
"""
function _make_spatial_vec(N::Int, d::Float64)
    return collect(((0:N-1) .- (N - 1) / 2.0) .* d)
end

"""
    _make_shift_operators(kv, d)

Compute staggered grid shift operators for wavenumber vector `kv` and spacing `d`.
Returns (shift_pos, shift_neg) where:
  shift_pos = exp(+im * kv * d/2)
  shift_neg = exp(-im * kv * d/2)
"""
function _make_shift_operators(kv::Vector{Float64}, d::Float64)
    shift_pos = exp.(im .* kv .* d / 2)
    shift_neg = exp.(-im .* kv .* d / 2)
    return shift_pos, shift_neg
end

# ============================================================================
# Constructors
# ============================================================================

"""
    KWaveGrid(Nx, dx) -> KWaveGrid1D
    KWaveGrid(Nx, dx, Ny, dy) -> KWaveGrid2D
    KWaveGrid(Nx, dx, Ny, dy, Nz, dz) -> KWaveGrid3D

Create a k-Wave computational grid.

# Arguments
- `Nx`, `Ny`, `Nz`: Number of grid points in each dimension
- `dx`, `dy`, `dz`: Grid spacing in each dimension [m]

# Returns
Returns [`KWaveGrid1D`](@ref), [`KWaveGrid2D`](@ref), or [`KWaveGrid3D`](@ref) based on
the number of arguments. Call [`make_time!`](@ref) on the result before running a solver.

# Example
```julia
kgrid = KWaveGrid(256, 0.1e-3, 256, 0.1e-3)
make_time!(kgrid, 1500.0)
```

# See Also
[`make_time!`](@ref), [`k_max`](@ref), [`total_grid_points`](@ref), [`get_optimal_pml_size`](@ref)
"""
function KWaveGrid(Nx::Int, dx::Real)
    dx = Float64(dx)
    x_vec = _make_spatial_vec(Nx, dx)
    kx_vec = _make_wavenumber_vec(Nx, dx)
    k = abs.(kx_vec)
    ddx_pos, ddx_neg = _make_shift_operators(kx_vec, dx)
    return KWaveGrid1D(
        Nx, dx, x_vec, kx_vec, k,
        ddx_pos, ddx_neg,
        Ref(0.0), Ref(0), Float64[]
    )
end

function KWaveGrid(Nx::Int, dx::Real, Ny::Int, dy::Real)
    dx, dy = Float64(dx), Float64(dy)
    x_vec = _make_spatial_vec(Nx, dx)
    y_vec = _make_spatial_vec(Ny, dy)
    kx_vec = _make_wavenumber_vec(Nx, dx)
    ky_vec = _make_wavenumber_vec(Ny, dy)

    # Wavenumber magnitude: k[i,j] = sqrt(kx[i]^2 + ky[j]^2)
    k = sqrt.(kx_vec.^2 .+ ky_vec'.^2)

    ddx_pos, ddx_neg = _make_shift_operators(kx_vec, dx)
    ddy_pos, ddy_neg = _make_shift_operators(ky_vec, dy)

    return KWaveGrid2D(
        Nx, dx, Ny, dy,
        x_vec, y_vec,
        kx_vec, ky_vec, k,
        ddx_pos, ddx_neg, ddy_pos, ddy_neg,
        Ref(0.0), Ref(0), Float64[]
    )
end

function KWaveGrid(Nx::Int, dx::Real, Ny::Int, dy::Real, Nz::Int, dz::Real)
    dx, dy, dz = Float64(dx), Float64(dy), Float64(dz)
    x_vec = _make_spatial_vec(Nx, dx)
    y_vec = _make_spatial_vec(Ny, dy)
    z_vec = _make_spatial_vec(Nz, dz)
    kx_vec = _make_wavenumber_vec(Nx, dx)
    ky_vec = _make_wavenumber_vec(Ny, dy)
    kz_vec = _make_wavenumber_vec(Nz, dz)

    # Wavenumber magnitude: k[i,j,k] = sqrt(kx[i]^2 + ky[j]^2 + kz[k]^2)
    k = sqrt.(
        reshape(kx_vec, :, 1, 1).^2 .+
        reshape(ky_vec, 1, :, 1).^2 .+
        reshape(kz_vec, 1, 1, :).^2
    )

    ddx_pos, ddx_neg = _make_shift_operators(kx_vec, dx)
    ddy_pos, ddy_neg = _make_shift_operators(ky_vec, dy)
    ddz_pos, ddz_neg = _make_shift_operators(kz_vec, dz)

    return KWaveGrid3D(
        Nx, dx, Ny, dy, Nz, dz,
        x_vec, y_vec, z_vec,
        kx_vec, ky_vec, kz_vec, k,
        ddx_pos, ddx_neg, ddy_pos, ddy_neg, ddz_pos, ddz_neg,
        Ref(0.0), Ref(0), Float64[]
    )
end

# ============================================================================
# Convenience methods
# ============================================================================

Base.ndims(g::KWaveGrid1D) = 1
Base.ndims(g::KWaveGrid2D) = 2
Base.ndims(g::KWaveGrid3D) = 3

"""
    total_grid_points(grid)

Return the total number of grid points.
"""
total_grid_points(g::KWaveGrid1D) = g.Nx
total_grid_points(g::KWaveGrid2D) = g.Nx * g.Ny
total_grid_points(g::KWaveGrid3D) = g.Nx * g.Ny * g.Nz

"""
    grid_size(grid)

Return a tuple of the grid dimensions.
"""
grid_size(g::KWaveGrid1D) = (g.Nx,)
grid_size(g::KWaveGrid2D) = (g.Nx, g.Ny)
grid_size(g::KWaveGrid3D) = (g.Nx, g.Ny, g.Nz)

"""
    grid_spacing(grid)

Return a tuple of the grid spacings.
"""
grid_spacing(g::KWaveGrid1D) = (g.dx,)
grid_spacing(g::KWaveGrid2D) = (g.dx, g.dy)
grid_spacing(g::KWaveGrid3D) = (g.dx, g.dy, g.dz)

"""
    k_max(grid)

Return the maximum supported wavenumber (Nyquist limit) [rad/m].

`k_max = π / dx` (for the finest-spaced dimension).
Used internally for stability checks; also indicates the highest spatial frequency
resolvable on the grid.

# See Also
[`KWaveGrid`](@ref), [`make_time!`](@ref), [`total_grid_points`](@ref)
"""
k_max(g::KWaveGrid1D) = π / g.dx
k_max(g::KWaveGrid2D) = min(π / g.dx, π / g.dy)
k_max(g::KWaveGrid3D) = min(π / g.dx, π / g.dy, π / g.dz)

# ============================================================================
# make_time!
# ============================================================================

"""
    make_time!(grid, sound_speed; cfl=0.3, t_end=nothing)

Compute the time array for the simulation based on the CFL stability condition.

# Arguments
- `grid`: KWaveGrid (any dimension)
- `sound_speed`: Maximum sound speed in the medium [m/s] (scalar or array)
- `cfl`: CFL number (default: 0.3, must be ≤ 1 for stability)
- `t_end`: End time [s] (if nothing, computed from grid diagonal / sound speed)

# Notes
Modifies `grid.dt[]`, `grid.Nt[]`, and `grid.t_array` in-place.
Pass `maximum(medium.sound_speed)` as `sound_speed` for heterogeneous media.

# See Also
[`KWaveGrid`](@ref), [`kspace_first_order`](@ref), [`k_max`](@ref)
"""
function make_time!(grid::AbstractKWaveGrid, sound_speed; cfl::Float64=0.3, t_end::Union{Nothing, Float64}=nothing)
    # Get maximum sound speed
    c_max = maximum(sound_speed)
    c_min = minimum(sound_speed)

    # Get minimum grid spacing
    ds = grid_spacing(grid)
    d_min = minimum(ds)

    # Compute time step from CFL condition
    dt = cfl * d_min / c_max
    grid.dt[] = dt

    # Compute end time if not specified
    # Default: time for wave to traverse the grid diagonal twice
    if t_end === nothing
        gs = grid_size(grid)
        # Grid diagonal length
        diag = sqrt(sum((n * d)^2 for (n, d) in zip(gs, ds)))
        # Time for slowest wave to cross diagonal
        t_end = 2.0 * diag / c_min
    end

    # Number of time steps
    Nt = ceil(Int, t_end / dt) + 1
    grid.Nt[] = Nt

    # Update time array
    resize!(grid.t_array, Nt)
    grid.t_array .= (0:Nt-1) .* dt

    return dt, Nt
end
