# ============================================================================
# Example: Delay-and-Sum Beamforming (2D)
# ============================================================================
# Demonstrates delay-and-sum beamforming reconstruction from
# simulated ultrasound sensor data.

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

# Source — point scatterer
p0 = zeros(Nx, Ny)
p0[Nx÷2, Ny÷2] = 1.0
source = KWaveSource(p0=p0)

# Sensor — linear array at x = 5
sensor_mask = falses(Nx, Ny)
sensor_mask[5, 20:108] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Forward simulation
println("Running forward simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)

# Beamforming reconstruction
println("Running delay-and-sum beamforming...")
sensor_positions = collect(20:108) .* dy  # sensor y-positions in meters

grid_x = range(10*dx, (Nx-10)*dx, length=100)
grid_z = range(10*dy, (Ny-10)*dy, length=100)

image = beamform_delay_and_sum(
    result[:p], sensor_positions, c0, dt,
    collect(grid_x), collect(grid_z);
    apodization=:hann,
)

println("Beamforming complete!")
println("  Image size: ", size(image))
println("  Peak value: ", maximum(abs, image))
