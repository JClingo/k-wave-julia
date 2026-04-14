# Grid

## Why Uniform Grids?

The k-space pseudospectral method evaluates spatial derivatives exactly in the Fourier domain.
This requires a **uniform Cartesian grid** — non-uniform spacing would break the DFT basis.
All three dimensions may have different spacings (`dx ≠ dy ≠ dz`), but spacing within each dimension is constant.

## Creating a Grid

Use the `KWaveGrid` constructor; it dispatches on argument count:

```julia
# 1D
kgrid = KWaveGrid(Nx, dx)

# 2D
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# 3D
kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz)
```

Grid spacings must satisfy the **points-per-wavelength** requirement:

```
dx ≤ c_min / (2 * f_max)
```

A minimum of **3–4 points per wavelength** at the highest frequency of interest is recommended.

## Setting the Time Step

```julia
make_time!(kgrid, c_max)
```

`make_time!` computes the time step and number of steps using the **CFL condition**:

```
dt = CFL * dx / c_max     (CFL = 0.3 by default)
Nt = ceil(t_end / dt)
```

The grid fields `kgrid.dt[]`, `kgrid.Nt[]`, and `kgrid.t_array` are populated in-place.

!!! note "Call make_time! before running a solver"
    All solver functions require `kgrid.dt[]` to be set. A call to `kspace_first_order` before
    `make_time!` will throw `UndefRefError`.

## Wavenumber Layout

Wavenumber vectors follow the FFTW / MATLAB convention:

```
kx_vec = [0, 1, ..., Nx/2-1, -Nx/2, ..., -1] * (2π / (Nx * dx))   (even Nx)
```

The wavenumber magnitude `kgrid.k` is:

| Dimension | Shape | Definition |
|---|---|---|
| 1D | `(Nx,)` | `|kx|` |
| 2D | `(Nx, Ny)` | `sqrt(kx² + ky²)` |
| 3D | `(Nx, Ny, Nz)` | `sqrt(kx² + ky² + kz²)` |

`k_max(kgrid)` returns the maximum resolvable wavenumber: `π / dx` (for the finest-spaced dimension).
This is used internally for stability checks.

## Example — 3D Grid

```julia
Nx, dx = 128, 0.5e-3   # 6.4 cm at 0.5 mm
Ny, dy = 128, 0.5e-3
Nz, dz = 64, 1.0e-3    # coarser in z

kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz)
make_time!(kgrid, 1500.0; t_end=50e-6)

println("dt  = ", kgrid.dt[], " s")
println("Nt  = ", kgrid.Nt[])
println("k_max = ", k_max(kgrid), " rad/m")
```

## See Also

[`KWaveGrid`](@ref), [`make_time!`](@ref), [`k_max`](@ref), [`get_optimal_pml_size`](@ref)
