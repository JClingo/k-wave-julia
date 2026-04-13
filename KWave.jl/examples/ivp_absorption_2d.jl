# ============================================================================
# Example: Acoustic Absorption in 2D
# ============================================================================
# Simulates wave propagation with power-law acoustic absorption
# and dispersion.
# Matches MATLAB k-Wave example: example_ivp_photoacoustic_waveforms.m

using KWave

# Grid
Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# Medium with absorption
medium = KWaveMedium(
    sound_speed=1500.0,
    density=1000.0,
    alpha_coeff=0.75,       # dB/(MHz^y cm)
    alpha_power=1.5,        # power law exponent
    alpha_mode=:stokes,     # include dispersion
)

# Time
make_time!(kgrid, 1500.0)

# Source
p0 = zeros(Nx, Ny)
disc = make_disc(Nx, Ny, Nx÷2, Ny÷2, 8)
p0[disc] .= 1.0
source = KWaveSource(p0=p0)

# Sensor
sensor_mask = falses(Nx, Ny)
sensor_mask[1, :] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max])

# Run
result = kspace_first_order(kgrid, medium, source, sensor)

println("Absorption simulation complete!")
println("  Peak pressure (with absorption): ", maximum(result[:p_max]))
