# Solvers

KWave.jl provides solvers for acoustic, elastic, and thermal wave propagation.

## Solver Selection Guide

| Use case | Solver |
|---|---|
| 1D/2D/3D acoustic time-domain | [`kspace_first_order`](@ref) |
| Axisymmetric (cylindrical) acoustic | [`kspace_first_order_as`](@ref) |
| Continuous-wave pressure field (Green's function) | [`acoustic_field_propagator`](@ref) |
| CW field via angular spectrum | [`angular_spectrum_cw`](@ref) |
| Elastic (P + S) waves 2D | [`pstd_elastic_2d`](@ref) |
| Elastic (P + S) waves 3D | [`pstd_elastic_3d`](@ref) |
| Thermal diffusion / bioheat | [`kwave_diffusion`](@ref) |
| Analytical bioheat solution | [`bioheat_exact`](@ref) |

---

## Acoustic Time-Domain Solvers

### kspace_first_order

**Purpose** — Run a k-space pseudospectral time-domain simulation in 1D, 2D, or 3D.
Dispatches on `kgrid` type; the same call signature works for all dimensionalities.

**Syntax**
```julia
output = kspace_first_order(kgrid, medium, source, sensor)
output = kspace_first_order(kgrid, medium, source, sensor; pml_size=20, data_cast=Float32)
```

```@docs
kspace_first_order
```

**Notes**
- Call [`make_time!`](@ref) on `kgrid` before this function.
- `pml_inside=true` (default) keeps the PML within the declared grid; `false` expands the grid, increasing memory use.
- Set `data_cast=Float32` when using GPU backends (CUDA/Metal/AMDGPU).

---

### kspace_first_order_as

**Purpose** — Axisymmetric (2D-AS) acoustic simulation using a cylindrical coordinate system.
The x-axis is the axial direction; the y-axis is the radial direction (y ≥ 0).

```@docs
kspace_first_order_as
```

---

## CW / Frequency-Domain Solvers

### acoustic_field_propagator

**Purpose** — Compute the steady-state CW pressure field using a Green's function approach.
Returns amplitude and phase maps at a single frequency. Faster than time-domain for CW-only problems.

```@docs
acoustic_field_propagator
```

---

### angular_spectrum_cw

**Purpose** — Propagate a CW pressure field using the angular spectrum method.
Valid in the paraxial (forward-propagating) regime.

```@docs
angular_spectrum_cw
```

---

## Elastic Wave Solvers

### pstd_elastic_2d

**Purpose** — 2D elastic wave simulation (P and S waves) using the pseudospectral time-domain method.
Requires [`ElasticMedium`](@ref) and [`ElasticSource`](@ref).

```@docs
pstd_elastic_2d
```

---

### pstd_elastic_3d

**Purpose** — 3D elastic wave simulation.

```@docs
pstd_elastic_3d
```

---

## Thermal / Diffusion Solvers

### kwave_diffusion

**Purpose** — Simulate thermal diffusion using the Pennes bioheat equation.
Requires [`ThermalMedium`](@ref) and [`ThermalSource`](@ref).

```@docs
kwave_diffusion
```

---

### bioheat_exact

**Purpose** — Analytical solution to the Pennes bioheat equation for homogeneous tissue.
Useful for validation against [`kwave_diffusion`](@ref).

```@docs
bioheat_exact
```
