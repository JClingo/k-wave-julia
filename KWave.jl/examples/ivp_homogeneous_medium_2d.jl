# ============================================================================
# Example: Initial Value Problem in a Homogeneous 2D Medium
# ============================================================================
# Simulates the propagation of an initial pressure distribution (disc)
# in a 2D homogeneous medium. This is the simplest k-Wave example.
# Matches MATLAB k-Wave example: example_ivp_homogeneous_medium_2D.m

using KWave

# Grid parameters
Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# Medium properties (water)
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

# Time array
make_time!(kgrid, 1500.0)

# Initial pressure distribution — a small disc
p0 = zeros(Nx, Ny)
disc = make_disc(Nx, Ny, Nx÷2, Ny÷2, 5)
p0[disc] .= 1.0

source = KWaveSource(p0=p0)

# Sensor — record along a line
sensor_mask = falses(Nx, Ny)
sensor_mask[1, :] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max])

# Run simulation
result = kspace_first_order(kgrid, medium, source, sensor)

println("Simulation complete!")
println("  Max pressure recorded: ", maximum(result[:p_max]))
println("  Sensor data size: ", size(result[:p]))
