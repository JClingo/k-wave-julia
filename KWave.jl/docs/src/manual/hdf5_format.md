# HDF5 File Format

KWave.jl uses an HDF5 file format compatible with the [k-Wave C++ binary solver](http://www.k-wave.org/documentation/k-wave_cpp_solver_user_manual.pdf).
This allows using Julia for simulation setup/post-processing while offloading computation to the
high-performance C++ binary.

## Writing an Input File

Pass `save_to_disk` to `kspace_first_order`:

```julia
kspace_first_order(kgrid, medium, source, sensor; save_to_disk="input.h5")
```

This calls the following write functions in sequence:

| Function | Datasets written |
|---|---|
| [`write_matrix`](@ref) | Source arrays (`p_source_input`, `ux_source_input`, …) |
| [`write_grid`](@ref) | `Nx`, `Ny`, `Nz`, `dx`, `dy`, `dz`, `dt`, `Nt` |
| [`write_flags`](@ref) | `p_source_flag`, `u_source_flag`, `p0_source_flag`, `sensor_mask_type`, … |
| [`write_attributes`](@ref) | `created_by`, `file_description`, `major_version`, `minor_version` |

## Running the C++ Binary

```bash
kspaceFirstOrder3D-OMP -i input.h5 -o output.h5 --p_rms --p_max
```

See the k-Wave C++ documentation for available flags.

## Reading the Output

```julia
output = read_output("output.h5", sensor)
p_rms = output[:p_rms]
```

`read_output` maps C++ output dataset names to the same `Symbol` keys used by `SimulationOutput`.

## Dataset Reference

### Grid datasets (`write_grid`)

| Dataset | Type | Units |
|---|---|---|
| `Nx`, `Ny`, `Nz` | `uint64` | grid points |
| `dx`, `dy`, `dz` | `float32` | m |
| `dt` | `float32` | s |
| `Nt` | `uint64` | steps |

### Flag datasets (`write_flags`)

| Dataset | Values |
|---|---|
| `p_source_flag` | 0 or 1 |
| `u_source_flag` | 0 or 1 |
| `p0_source_flag` | 0 or 1 |
| `sensor_mask_type` | 1 = index, 2 = corners |
| `nonlinear_flag` | 0 or 1 |
| `absorbing_flag` | 0, 1, or 2 |

## See Also

[`write_matrix`](@ref), [`write_grid`](@ref), [`write_flags`](@ref),
[`write_attributes`](@ref), [`read_output`](@ref)
