# Running Simulations

## Basic Call

```julia
output = kspace_first_order(kgrid, medium, source, sensor)
```

Grid dimensionality is dispatched automatically; the same call works for 1D, 2D, and 3D grids.

## Keyword Arguments

| Keyword | Type | Default | Description |
|---|---|---|---|
| `pml_inside` | `Bool` | `true` | Place PML inside grid (true) or expand grid (false) |
| `pml_size` | `Int` | `20` | PML thickness in grid points |
| `pml_alpha` | `Float64` | `2.0` | PML absorption coefficient |
| `smooth_p0` | `Bool` | `true` | Smooth initial pressure |
| `smooth_c0` | `Bool` | `false` | Smooth sound speed |
| `smooth_rho0` | `Bool` | `false` | Smooth density |
| `data_cast` | `Type` | `Float64` | Float type for arrays (use `Float32` for GPU) |
| `save_to_disk` | `String\|Nothing` | `nothing` | Write HDF5 instead of running |
| `plot_sim` | `Bool` | `false` | Display real-time simulation progress |
| `plot_layout` | `Symbol` | `:default` | Display layout |
| `plot_scale` | `Symbol` | `:auto` | Display colour scale |
| `record_movie` | `String\|Nothing` | `nothing` | Path to save movie file |
| `progress_callback` | `Function\|Nothing` | `nothing` | Called each time step |

## Progress Callback

The `progress_callback` is called at every time step with signature:

```julia
callback(t_index::Int, Nt::Int, p_field::AbstractArray)
```

Example — print progress every 100 steps:

```julia
cb = (t, Nt, p) -> (t % 100 == 0) && @info "Step $t / $Nt"
output = kspace_first_order(kgrid, medium, source, sensor; progress_callback=cb)
```

## Saving to HDF5 (C++ Binary Mode)

To offload computation to the k-Wave C++ binary solver:

```julia
# Step 1: write input file
kspace_first_order(kgrid, medium, source, sensor; save_to_disk="input.h5")

# Step 2: run C++ binary externally
# $ kspaceFirstOrder3D-CUDA -i input.h5 -o output.h5

# Step 3: read results back
output = read_output("output.h5", sensor)
```

## GPU Acceleration

Load a GPU backend extension and set `data_cast=Float32`:

```julia
using KWave, CUDA

output = kspace_first_order(kgrid, medium, source, sensor;
                             data_cast=Float32)
```

Arrays are automatically moved to the GPU. Results are returned as CPU arrays.

## Axisymmetric Solver

Use [`kspace_first_order_as`](@ref) for cylindrically symmetric problems.
The x-axis is axial; y-axis is radial (y ≥ 0 only):

```julia
output = kspace_first_order_as(kgrid, medium, source, sensor)
```

## See Also

[`kspace_first_order`](@ref), [`kspace_first_order_as`](@ref),
[`make_time!`](@ref), [`KWaveMedium`](@ref), [`KWaveSource`](@ref), [`KWaveSensor`](@ref),
[`write_grid`](@ref), [`read_output`](@ref)
