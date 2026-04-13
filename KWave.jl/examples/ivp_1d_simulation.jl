# ============================================================================
# Example: 1D Initial Value Problem
# ============================================================================
# Simplest possible k-Wave simulation — 1D wave propagation.
# Matches MATLAB k-Wave example: example_ivp_1D_simulation.m

using KWave

# Grid
Nx = 256; dx = 10e-6
kgrid = KWaveGrid(Nx, dx)

# Medium
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

# Time
make_time!(kgrid, 1500.0)

# Source — Gaussian initial pressure
p0 = zeros(Nx)
p0[Nx÷2] = 1.0
p0 = smooth(p0; restore_max=true)
source = KWaveSource(p0=p0)

# Sensor — two points
sensor_mask = falses(Nx)
sensor_mask[Nx÷4] = true
sensor_mask[3*Nx÷4] = true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Run
result = kspace_first_order(kgrid, medium, source, sensor)

println("1D simulation complete!")
println("  Sensor data size: ", size(result[:p]))
