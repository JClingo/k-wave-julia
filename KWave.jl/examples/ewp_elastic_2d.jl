# ============================================================================
# Example: Elastic Wave Propagation in 2D
# ============================================================================
# Simulates coupled compressional and shear wave propagation in
# a 2D elastic medium.
# Matches MATLAB k-Wave example: example_ewp_elastic_wave_2D.m

using KWave

# Grid
Nx = 128; dx = 0.5e-3
Ny = 128; dy = 0.5e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# Elastic medium — bone-like material
medium = ElasticMedium(
    sound_speed_compression=3000.0,  # P-wave speed [m/s]
    sound_speed_shear=1500.0,        # S-wave speed [m/s]
    density=1800.0,                  # [kg/m³]
)

# Time
make_time!(kgrid, 3000.0)

# Source — initial isotropic pressure (explosion source)
p0 = zeros(Nx, Ny)
disc = make_disc(Nx, Ny, Nx÷2, Ny÷2, 3)
p0[disc] .= 1.0
source = ElasticSource(p0=p0)

# Sensor — ring of points
sensor_mask = falses(Nx, Ny)
circle = make_circle(Nx, Ny, Nx÷2, Ny÷2, 40)
sensor_mask[circle] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Run
println("Running 2D elastic simulation...")
result = pstd_elastic_2d(kgrid, medium, source, sensor)

println("Elastic simulation complete!")
println("  Sensor data size: ", size(result[:p]))
