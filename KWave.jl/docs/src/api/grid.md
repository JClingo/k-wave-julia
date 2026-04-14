# Grid

The grid defines the spatial discretisation and pre-computes all wavenumber operators needed by the solvers.

## Grid Types

| Type | Dimensions | Fields |
|---|---|---|
| [`KWaveGrid1D`](@ref) | 1D | `Nx`, `dx`, `x_vec`, `kx_vec`, `k` |
| [`KWaveGrid2D`](@ref) | 2D | `Nx/Ny`, `dx/dy`, `x/y_vec`, `kx/ky_vec`, `k` |
| [`KWaveGrid3D`](@ref) | 3D | `Nx/Ny/Nz`, `dx/dy/dz`, `x/y/z_vec`, `kx/ky/kz_vec`, `k` |

All grid types share mutable time fields `dt`, `Nt`, `t_array` (set by `make_time!`).

## Public Fields (all grid types)

| Field | Type | Description |
|---|---|---|
| `Nx`, `Ny`, `Nz` | `Int` | Grid point counts per dimension |
| `dx`, `dy`, `dz` | `Float64` | Grid spacing [m] |
| `x_vec`, `y_vec`, `z_vec` | `Vector{Float64}` | Centred spatial coordinate vectors [m] |
| `kx_vec`, `ky_vec`, `kz_vec` | `Vector{Float64}` | Wavenumber vectors [rad/m] |
| `k` | `Array{Float64}` | Wavenumber magnitude (shape matches grid) |
| `dt` | `Ref{Float64}` | Time step [s] (set by `make_time!`) |
| `Nt` | `Ref{Int}` | Number of time steps (set by `make_time!`) |
| `t_array` | `Vector{Float64}` | Time vector [s] (set by `make_time!`) |

---

```@docs
AbstractKWaveGrid
KWaveGrid1D
KWaveGrid2D
KWaveGrid3D
KWaveGrid
make_time!
total_grid_points
grid_size
grid_spacing
k_max
```
