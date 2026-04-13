# ============================================================================
# Example: Initial Value Problem in a Heterogeneous 2D Medium
# ============================================================================
# Simulates propagation through a medium with spatially varying
# sound speed and density.
# Matches MATLAB k-Wave example: example_ivp_heterogeneous_medium_2D.m

using KWave

# Grid
Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# Heterogeneous medium — higher speed lens in the center
sound_speed = fill(1500.0, Nx, Ny)
density = fill(1000.0, Nx, Ny)

# Create a circular lens region with higher sound speed
lens = make_disc(Nx, Ny, Nx÷2 + 20, Ny÷2, 15)
sound_speed[lens] .= 1800.0
density[lens] .= 1200.0

medium = KWaveMedium(sound_speed=sound_speed, density=density)

# Time array
make_time!(kgrid, sound_speed)

# Source — initial pressure disc
p0 = zeros(Nx, Ny)
disc = make_disc(Nx, Ny, Nx÷4, Ny÷2, 5)
p0[disc] .= 1.0
source = KWaveSource(p0=p0)

# Sensor — full boundary
sensor_mask = falses(Nx, Ny)
sensor_mask[1, :] .= true
sensor_mask[end, :] .= true
sensor_mask[:, 1] .= true
sensor_mask[:, end] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Run
result = kspace_first_order(kgrid, medium, source, sensor)

println("Heterogeneous simulation complete!")
println("  Sensor data size: ", size(result[:p]))
