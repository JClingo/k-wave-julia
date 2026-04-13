# ============================================================================
# Example: Sensor Directivity in 2D
# ============================================================================
# Demonstrates the effect of finite-sized sensor elements on recorded
# pressure data using sensor directivity.
# Matches MATLAB k-Wave example: example_sd_directivity_modelling_2D.m

using KWave

# Grid
Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# Medium
c0 = 1500.0
medium = KWaveMedium(sound_speed=c0, density=1000.0)

# Time
make_time!(kgrid, c0)

# Source — point source
p0 = zeros(Nx, Ny)
p0[Nx÷2, Ny÷2] = 1.0
source = KWaveSource(p0=p0)

# Sensor with directivity — line sensor at x=10
sensor_mask = falses(Nx, Ny)
sensor_mask[10, :] .= true

# Directivity: sensor elements are sensitive along y-axis
# directivity_angle = angle of maximum sensitivity
directivity_angle = zeros(Nx, Ny)
directivity_angle[10, :] .= π/2  # sensitive in y direction

sensor = KWaveSensor(
    mask=sensor_mask,
    record=[:p],
    directivity_angle=directivity_angle,
    directivity_size=5 * dy,
)

# Run
result = kspace_first_order(kgrid, medium, source, sensor)

println("Directivity simulation complete!")
println("  Sensor data size: ", size(result[:p]))
