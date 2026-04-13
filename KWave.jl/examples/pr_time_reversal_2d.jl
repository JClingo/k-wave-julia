# ============================================================================
# Example: Photoacoustic Reconstruction via Time Reversal (2D)
# ============================================================================
# Demonstrates photoacoustic image reconstruction using time reversal.
# First, forward-propagates an initial pressure distribution and records
# sensor data, then uses time reversal to reconstruct the original image.
# Matches MATLAB k-Wave example: example_pr_2D_TR_line_sensor.m

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

# Source — two small discs (photoacoustic sources)
p0 = zeros(Nx, Ny)
disc1 = make_disc(Nx, Ny, 40, 50, 4)
disc2 = make_disc(Nx, Ny, 80, 70, 6)
p0[disc1] .= 1.0
p0[disc2] .= 0.5
source = KWaveSource(p0=p0)

# Sensor — line sensor at x = 1
sensor_mask = falses(Nx, Ny)
sensor_mask[1, :] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Forward simulation
println("Running forward simulation...")
result_fwd = kspace_first_order(kgrid, medium, source, sensor)

# Time reversal reconstruction
println("Running time reversal reconstruction...")
sensor_tr = KWaveSensor(
    mask=sensor_mask,
    time_reversal_boundary_data=result_fwd[:p],
)
source_tr = KWaveSource()  # no source for reconstruction

result_tr = kspace_first_order(kgrid, medium, source_tr, sensor_tr)

println("Reconstruction complete!")
println("  Original max p0: ", maximum(p0))
println("  Reconstructed max: ", maximum(result_tr[:p_final]))
