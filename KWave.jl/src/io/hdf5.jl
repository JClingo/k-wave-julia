# ============================================================================
# KWave.jl — HDF5 I/O (compatible with k-Wave C++ binary format)
# ============================================================================

"""
    write_matrix(filename, dataset_name, data; compression=true)

Write an array to an HDF5 file as a named dataset.

Julia arrays are already column-major (Fortran order), matching k-Wave's
C++ binary expectations, so no transposition is needed.

# Arguments
- `filename`: Path to HDF5 file
- `dataset_name`: Name of the dataset within the file
- `data`: Array to write
- `compression`: Whether to use gzip compression (default: true)
"""
function write_matrix(filename::AbstractString, dataset_name::AbstractString,
                      data::AbstractArray; compression::Bool=true)
    h5open(filename, isfile(filename) ? "r+" : "w") do fid
        if haskey(fid, dataset_name)
            delete_object(fid, dataset_name)
        end
        if compression && length(data) > 1
            fid[dataset_name, chunk=size(data), deflate=3] = data
        else
            fid[dataset_name] = data
        end
        # Write data type attribute
        attrs(fid[dataset_name])["data_type"] = eltype(data) <: Integer ? "long" : "float"
    end
end

# Scalar overload
function write_matrix(filename::AbstractString, dataset_name::AbstractString,
                      data::Real; compression::Bool=true)
    write_matrix(filename, dataset_name, [data]; compression=compression)
end

"""
    write_grid(filename, kgrid; pml_size, pml_alpha)

Write grid dimensions, spacing, time parameters, and PML configuration
to an HDF5 file in the k-Wave C++ binary format.

# Arguments
- `filename`: Path to HDF5 file
- `kgrid`: AbstractKWaveGrid
- `pml_size`: PML size (Int or Tuple)
- `pml_alpha`: PML absorption coefficient (Float64 or Tuple)
"""
function write_grid(filename::AbstractString, kgrid::AbstractKWaveGrid;
                    pml_size, pml_alpha)
    nd = ndims(kgrid)
    gs = grid_size(kgrid)
    ds = grid_spacing(kgrid)

    # Normalize pml_size and pml_alpha to tuples
    pml_sizes = pml_size isa Tuple ? pml_size : ntuple(_ -> Int(pml_size), nd)
    pml_alphas = pml_alpha isa Tuple ? pml_alpha : ntuple(_ -> Float64(pml_alpha), nd)

    h5open(filename, isfile(filename) ? "r+" : "w") do fid
        # Grid dimensions
        fid["Nx"] = Int64[gs[1]]
        if nd >= 2; fid["Ny"] = Int64[gs[2]]; end
        if nd >= 3; fid["Nz"] = Int64[gs[3]]; end

        # Grid spacing
        fid["dx"] = Float64[ds[1]]
        if nd >= 2; fid["dy"] = Float64[ds[2]]; end
        if nd >= 3; fid["dz"] = Float64[ds[3]]; end

        # Time parameters
        fid["Nt"] = Int64[kgrid.Nt[]]
        fid["dt"] = Float64[kgrid.dt[]]

        # PML parameters
        dim_names = ["x", "y", "z"]
        for d in 1:nd
            fid["pml_$(dim_names[d])_size"] = Int64[pml_sizes[d]]
            fid["pml_$(dim_names[d])_alpha"] = Float64[pml_alphas[d]]
        end
    end
end

"""
    write_flags(filename, flags)

Write simulation configuration flags to an HDF5 file.

# Arguments
- `filename`: Path to HDF5 file
- `flags`: Dict{String, Any} of flag names and values
"""
function write_flags(filename::AbstractString, flags::Dict{String, <:Any})
    h5open(filename, isfile(filename) ? "r+" : "w") do fid
        for (key, value) in flags
            if haskey(fid, key)
                delete_object(fid, key)
            end
            if value isa AbstractString
                fid[key] = value
            elseif value isa Integer
                fid[key] = Int64[value]
            elseif value isa AbstractFloat
                fid[key] = Float64[value]
            else
                fid[key] = value
            end
        end
    end
end

"""
    write_attributes(filename)

Write file-level metadata attributes to an HDF5 file.

Writes:
- `created_by`: "KWave.jl"
- `creation_date`: Current date/time string
- `file_type`: "input"
- `major_version`: "1"
- `minor_version`: "2"
"""
function write_attributes(filename::AbstractString)
    h5open(filename, isfile(filename) ? "r+" : "w") do fid
        attrs(fid)["created_by"] = "KWave.jl"
        attrs(fid)["creation_date"] = string(Dates.now())
        attrs(fid)["file_type"] = "input"
        attrs(fid)["major_version"] = "1"
        attrs(fid)["minor_version"] = "2"
    end
end

# We need Dates for the timestamp
using Dates

"""
    read_matrix(filename, dataset_name)

Read an array from an HDF5 file.

# Arguments
- `filename`: Path to HDF5 file
- `dataset_name`: Name of the dataset within the file

# Returns
Array read from the file.
"""
function read_matrix(filename::AbstractString, dataset_name::AbstractString)
    h5open(filename, "r") do fid
        return read(fid[dataset_name])
    end
end

"""
    read_output(filename)

Read a k-Wave simulation output HDF5 file and return a SimulationOutput.

# Arguments
- `filename`: Path to HDF5 output file

# Returns
`SimulationOutput` containing all recorded fields.
"""
function read_output(filename::AbstractString)
    data = Dict{Symbol, AbstractArray}()
    h5open(filename, "r") do fid
        for name in keys(fid)
            data[Symbol(name)] = read(fid[name])
        end
    end
    return SimulationOutput(data)
end
