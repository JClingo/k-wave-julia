# ============================================================================
# Example: Ultrasound Phased Array Simulation (3D)
# ============================================================================
# Simulates a focused ultrasound beam from a linear phased array
# transducer in 3D.
# Matches MATLAB k-Wave example: example_us_defining_transducer.m

using KWave

# Grid — smaller for demonstration
Nx = 64; dx = 0.5e-3
Ny = 64; dy = 0.5e-3
Nz = 64; dz = 0.5e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz)

# Medium
c0 = 1500.0
medium = KWaveMedium(sound_speed=c0, density=1000.0)

# Time
make_time!(kgrid, c0)

# Transducer definition
transducer = KWaveTransducer(
    number_elements=32,
    element_width=1,
    element_length=12,
    element_spacing=0,
    position=(1, Ny÷2 - 16, Nz÷2 - 6),
    radius=Inf,
    focus_distance=30.0 * dx,
    steering_angle=0.0,
    transmit_apodization=:hann,
    receive_apodization=:rectangular,
)

# Get source from transducer
dt = kgrid.dt[]
Nt = kgrid.Nt[]
freq = 1e6
signal = tone_burst(1/dt, freq, 3)
if length(signal) < Nt
    signal = vcat(signal, zeros(Nt - length(signal)))
end

source_mask = get_transducer_binary_mask(transducer, kgrid)

# Create source
p_mask_bool = source_mask .> 0
n_sources = count(p_mask_bool)
p_signal = repeat(signal', n_sources, 1)

source = KWaveSource(p_mask=p_mask_bool, p=Float64.(p_signal[:, 1:Nt]))

# Sensor — single plane
sensor_mask = falses(Nx, Ny, Nz)
sensor_mask[Nx÷2, :, :] .= true
sensor = KWaveSensor(mask=sensor_mask, record=[:p_max])

# Run
println("Running 3D phased array simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)

println("3D simulation complete!")
println("  Peak pressure: ", maximum(result[:p_max]))
