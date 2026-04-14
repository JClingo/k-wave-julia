# I/O

KWave.jl uses HDF5 files compatible with the k-Wave C++ binary solver format.
Pass `save_to_disk="/path/to/input.h5"` to [`kspace_first_order`](@ref) to write an input file
instead of running the simulation in Julia. The C++ binary can then be run separately and the
output read back with [`read_output`](@ref).

## Writing

```@docs
write_matrix
write_grid
write_flags
write_attributes
```

## Reading

```@docs
read_matrix
read_output
```

## File Format

The HDF5 file layout mirrors the k-Wave C++ binary specification:

| Dataset | Written by | Content |
|---|---|---|
| `ux_source_input`, `p_source_input`, … | `write_matrix` | Source time series arrays |
| `Nx`, `Ny`, `Nz`, `dx`, `dy`, `dz`, `dt`, `Nt` | `write_grid` | Grid parameters |
| `p_source_flag`, `u_source_flag`, … | `write_flags` | Simulation flags |
| `created_by`, `file_description`, … | `write_attributes` | Metadata strings |
