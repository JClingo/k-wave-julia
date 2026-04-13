# ============================================================================
# Example: Time-Varying Source Problem in a Homogeneous 2D Medium
# ============================================================================
# Simulates a pulsed ultrasound source (tone burst) transmitting into
# a homogeneous medium.
# Matches MATLAB k-Wave example: example_tvsp_homogeneous_medium_2D.m

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

# Source — line source with tone burst signal
freq = 1e6  # 1 MHz
signal = tone_burst(1/dt, freq, 3)

# Zero-pad signal to match Nt
if length(signal) < Nt
    signal = vcat(signal, zeros(Nt - length(signal)))
else
    signal = signal[1:Nt]
end

p_mask = falses(Nx, Ny)
p_mask[Nx÷4, Ny÷4:3*Ny÷4] .= true
n_sources = count(p_mask)

# Each source point gets the same signal
p_signal = repeat(signal', n_sources, 1)

source = KWaveSource(p_mask=p_mask, p=p_signal)

# Sensor
sensor_mask = falses(Nx, Ny)
sensor_mask[3*Nx÷4, :] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max])

# Run
result = kspace_first_order(kgrid, medium, source, sensor)

println("Time-varying source simulation complete!")
println("  Max pressure: ", maximum(result[:p_max]))
