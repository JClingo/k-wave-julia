# ============================================================================
# KWave.jl — PackageCompiler.jl standalone binary build script
# ============================================================================
# Run this script to compile KWave.jl into a standalone executable:
#
#   julia --project=.. build/build_standalone.jl [--output DIR]
#
# Prerequisites:
#   ] add PackageCompiler
#
# The compiled binary provides the KWave CLI interface:
#   ./kwave-cli run --input sim.h5 --output results.h5
#   ./kwave-cli validate --input sim.h5
#   ./kwave-cli info --input sim.h5
# ============================================================================

using Pkg

# Ensure PackageCompiler is available
try
    @eval using PackageCompiler
catch
    @error "PackageCompiler.jl is required. Install it with: ] add PackageCompiler"
    exit(1)
end

# Parse arguments
output_dir = "build/kwave-cli"
for i in eachindex(ARGS)
    if ARGS[i] == "--output" && i < length(ARGS)
        output_dir = ARGS[i + 1]
    end
end

# Precompile statements to speed up the compiled binary
precompile_script = joinpath(@__DIR__, "precompile_execution.jl")

# Create precompile execution script if it doesn't exist
if !isfile(precompile_script)
    open(precompile_script, "w") do io
        write(io, """
        using KWave

        # Pre-warm grid construction
        grid1d = KWaveGrid(64, 1e-4)
        grid2d = KWaveGrid(32, 1e-4, 32, 1e-4)
        grid3d = KWaveGrid(16, 1e-4, 16, 1e-4, 16, 1e-4)

        # Pre-warm medium, source, sensor
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        make_time!(grid2d, 1500.0)

        source = KWaveSource(p0=zeros(32, 32))
        source.p0[16, 16] = 1.0

        sensor_mask = falses(32, 32)
        sensor_mask[1, :] .= true
        sensor = KWaveSensor(mask=sensor_mask)

        # Pre-warm a small simulation
        result = kspace_first_order(grid2d, medium, source, sensor)

        # Pre-warm utilities
        make_disc(32, 32, 16, 16, 5)
        tone_burst(1/grid2d.dt[], 1e6, 3)

        println("Precompilation complete.")
        """)
    end
end

println("Building KWave standalone executable...")
println("  Output directory: $output_dir")
println("  This may take several minutes...")

# Get the KWave.jl project directory
project_dir = dirname(@__DIR__)

create_app(
    project_dir,
    output_dir;
    precompile_execution_file=[precompile_script],
    executables=["kwave-cli" => "KWave.CLI.main"],
    force=true,
    include_lazy_artifacts=true,
)

println("\nBuild complete! Binary at: $(joinpath(output_dir, "bin", "kwave-cli"))")
println("Usage:")
println("  $output_dir/bin/kwave-cli run --input sim.h5 --output results.h5")
println("  $output_dir/bin/kwave-cli validate --input sim.h5")
println("  $output_dir/bin/kwave-cli info --input sim.h5")
