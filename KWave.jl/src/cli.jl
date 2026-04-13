# ============================================================================
# KWave.jl — CLI tool
# ============================================================================

module CLI

using ..KWave

"""
    main(args=ARGS)

Command-line interface entry point for KWave.jl.

# Subcommands
- `run`: Run a simulation from an HDF5 input file
- `validate`: Validate an HDF5 input file
- `info`: Display information about an HDF5 file

# Usage
```
julia --project -e 'using KWave; KWave.CLI.main()' -- run --input sim.h5 --output results.h5
julia --project -e 'using KWave; KWave.CLI.main()' -- validate --input sim.h5
julia --project -e 'using KWave; KWave.CLI.main()' -- info --input sim.h5
```
"""
function main(args::Vector{String}=ARGS)
    if isempty(args)
        _print_usage()
        return 1
    end

    cmd = lowercase(args[1])
    rest = args[2:end]

    if cmd == "run"
        return _cmd_run(rest)
    elseif cmd == "validate"
        return _cmd_validate(rest)
    elseif cmd == "info"
        return _cmd_info(rest)
    elseif cmd in ["help", "-h", "--help"]
        _print_usage()
        return 0
    else
        println(stderr, "Unknown command: $cmd")
        _print_usage()
        return 1
    end
end

function _print_usage()
    println("""
    KWave.jl — k-space pseudospectral acoustic simulation

    Usage:
      kwave run --input <file.h5> --output <results.h5> [options]
      kwave validate --input <file.h5>
      kwave info --input <file.h5>
      kwave help

    Run options:
      --input, -i     Input HDF5 file path (required)
      --output, -o    Output HDF5 file path (required)
      --pml-size      PML size in grid points (default: 20)
      --pml-alpha     PML absorption coefficient (default: 2.0)
      --verbose, -v   Verbose output
    """)
end

function _parse_args(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--") || startswith(arg, "-")
            key = lstrip(arg, '-')
            # Map short flags
            if key == "i"; key = "input"; end
            if key == "o"; key = "output"; end
            if key == "v"; key = "verbose"; end

            if key in ["verbose"]
                opts[key] = "true"
            elseif i < length(args)
                i += 1
                opts[key] = args[i]
            end
        end
        i += 1
    end
    return opts
end

function _cmd_run(args::Vector{String})
    opts = _parse_args(args)

    input_file = get(opts, "input", nothing)
    output_file = get(opts, "output", nothing)

    if input_file === nothing
        println(stderr, "Error: --input is required")
        return 1
    end
    if output_file === nothing
        println(stderr, "Error: --output is required")
        return 1
    end

    if !isfile(input_file)
        println(stderr, "Error: Input file not found: $input_file")
        return 1
    end

    pml_size = parse(Int, get(opts, "pml-size", "20"))
    pml_alpha = parse(Float64, get(opts, "pml-alpha", "2.0"))
    verbose = get(opts, "verbose", "false") == "true"

    verbose && println("Loading input file: $input_file")

    # Read simulation parameters from HDF5
    try
        h5open(input_file, "r") do fid
            # Read grid dimensions
            Nx = read(fid["Nx"])[1]
            Ny = haskey(fid, "Ny") ? read(fid["Ny"])[1] : nothing
            Nz = haskey(fid, "Nz") ? read(fid["Nz"])[1] : nothing
            dx = read(fid["dx"])[1]
            dy = haskey(fid, "dy") ? read(fid["dy"])[1] : nothing
            dz = haskey(fid, "dz") ? read(fid["dz"])[1] : nothing
            Nt_val = read(fid["Nt"])[1]
            dt_val = read(fid["dt"])[1]

            verbose && println("Grid: Nx=$Nx" *
                (Ny !== nothing ? ", Ny=$Ny" : "") *
                (Nz !== nothing ? ", Nz=$Nz" : ""))

            # Create grid
            kgrid = if Nz !== nothing
                KWaveGrid(Int(Nx), dx, Int(Ny), dy, Int(Nz), dz)
            elseif Ny !== nothing
                KWaveGrid(Int(Nx), dx, Int(Ny), dy)
            else
                KWaveGrid(Int(Nx), dx)
            end

            # Set time parameters
            kgrid.dt[] = dt_val
            kgrid.Nt[] = Int(Nt_val)
            resize!(kgrid.t_array, Int(Nt_val))
            kgrid.t_array .= (0:Int(Nt_val)-1) .* dt_val

            # Read medium
            c0 = haskey(fid, "c0") ? read(fid["c0"]) : read(fid["c_ref"])[1]
            rho0 = haskey(fid, "rho0") ? read(fid["rho0"]) : 1.0

            medium = KWaveMedium(sound_speed=c0, density=rho0)

            # Read source
            p0 = haskey(fid, "p0") ? read(fid["p0"]) : nothing
            source = KWaveSource(p0=p0)

            # Sensor: record at all points
            gs = grid_size(kgrid)
            sensor_mask = trues(gs...)
            sensor = KWaveSensor(mask=sensor_mask, record=[:p_max, :p_rms])

            verbose && println("Running simulation (Nt=$Nt_val)...")

            result = kspace_first_order(kgrid, medium, source, sensor;
                                        pml_size=pml_size, pml_alpha=pml_alpha)

            # Write output
            verbose && println("Writing output to: $output_file")
            for key in keys(result)
                write_matrix(output_file, string(key), result[key])
            end
            write_attributes(output_file)

            verbose && println("Done.")
        end
    catch e
        println(stderr, "Error: $e")
        return 1
    end

    return 0
end

function _cmd_validate(args::Vector{String})
    opts = _parse_args(args)
    input_file = get(opts, "input", nothing)

    if input_file === nothing
        println(stderr, "Error: --input is required")
        return 1
    end

    if !isfile(input_file)
        println(stderr, "Error: File not found: $input_file")
        return 1
    end

    try
        h5open(input_file, "r") do fid
            required = ["Nx", "dx", "Nt", "dt"]
            missing_fields = String[]
            for field in required
                if !haskey(fid, field)
                    push!(missing_fields, field)
                end
            end

            if isempty(missing_fields)
                println("✓ Valid k-Wave input file")
                println("  Nx = ", read(fid["Nx"])[1])
                if haskey(fid, "Ny"); println("  Ny = ", read(fid["Ny"])[1]); end
                if haskey(fid, "Nz"); println("  Nz = ", read(fid["Nz"])[1]); end
                println("  Nt = ", read(fid["Nt"])[1])
                println("  dt = ", read(fid["dt"])[1])
                return 0
            else
                println(stderr, "✗ Invalid: missing fields: ", join(missing_fields, ", "))
                return 1
            end
        end
    catch e
        println(stderr, "Error reading file: $e")
        return 1
    end
end

function _cmd_info(args::Vector{String})
    opts = _parse_args(args)
    input_file = get(opts, "input", nothing)

    if input_file === nothing
        println(stderr, "Error: --input is required")
        return 1
    end

    if !isfile(input_file)
        println(stderr, "Error: File not found: $input_file")
        return 1
    end

    try
        h5open(input_file, "r") do fid
            println("File: $input_file")

            # File attributes
            file_attrs = attrs(fid)
            if length(file_attrs) > 0
                println("\nAttributes:")
                for key in keys(file_attrs)
                    println("  $key = ", read(file_attrs[key]))
                end
            end

            # Datasets
            println("\nDatasets:")
            for name in keys(fid)
                ds = fid[name]
                if ds isa HDF5.Dataset
                    sz = size(ds)
                    tp = eltype(ds)
                    println("  $name: $tp $sz")
                end
            end
        end
    catch e
        println(stderr, "Error: $e")
        return 1
    end

    return 0
end

# Allow using HDF5 functions
using HDF5

end # module CLI
