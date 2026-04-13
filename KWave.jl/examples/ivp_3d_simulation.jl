# ============================================================================
# Example: 3D Initial Value Problem
# ============================================================================
# 3D wave propagation from a ball-shaped initial pressure distribution.
# Matches MATLAB k-Wave example: example_ivp_3D_simulation.m

using KWave

# Grid — smaller for speed
Nx = 64; dx = 0.1e-3
Ny = 64; dy = 0.1e-3
Nz = 64; dz = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz)

# Medium
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

# Time
make_time!(kgrid, 1500.0)

# Source — ball initial pressure
p0 = zeros(Nx, Ny, Nz)
ball = make_ball(Nx, Ny, Nz, Nx÷2, Ny÷2, Nz÷2, 5)
p0[ball] .= 1.0
source = KWaveSource(p0=p0)

# Sensor — single plane
sensor_mask = falses(Nx, Ny, Nz)
sensor_mask[:, :, Nz÷2] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p_max])

# Run
println("Running 3D simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)

println("3D simulation complete!")
println("  Peak pressure: ", maximum(result[:p_max]))
