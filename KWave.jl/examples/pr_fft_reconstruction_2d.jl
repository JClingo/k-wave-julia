# ============================================================================
# Example: FFT-based Photoacoustic Reconstruction (2D)
# ============================================================================
# Uses the k-space FFT reconstruction method to reconstruct an initial
# pressure distribution from line sensor data.
# Matches MATLAB k-Wave example: example_pr_2D_FFT_line_sensor.m

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
dt = kgrid.dt[]
Nt = kgrid.Nt[]

# Source — phantom with two objects
p0 = zeros(Nx, Ny)
disc1 = make_disc(Nx, Ny, 50, 50, 5)
disc2 = make_disc(Nx, Ny, 80, 70, 3)
p0[disc1] .= 1.0
p0[disc2] .= 0.7
source = KWaveSource(p0=p0)

# Sensor — line at x = 1
sensor_mask = falses(Nx, Ny)
sensor_mask[1, :] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Forward simulation
println("Forward simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)

# FFT reconstruction
println("FFT reconstruction...")
p_recon = kspace_line_recon(result[:p], dy, dt; c=c0)

println("Reconstruction complete!")
println("  Reconstructed image size: ", size(p_recon))
println("  Max value: ", maximum(p_recon))
