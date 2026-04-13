# ============================================================================
# Example: Nonlinear Propagation in 1D
# ============================================================================
# Simulates finite-amplitude (nonlinear) acoustic wave propagation in 1D.
# Demonstrates harmonic generation due to the B/A nonlinearity parameter.
# Matches MATLAB k-Wave example: example_tvsp_nonlinear_propagation_1D.m

using KWave

# Grid
Nx = 512; dx = 10e-6
kgrid = KWaveGrid(Nx, dx)

# Medium with nonlinearity and absorption
medium = KWaveMedium(
    sound_speed=1500.0,
    density=1000.0,
    BonA=5.0,                     # Water-like nonlinearity
    alpha_coeff=0.75,             # [dB/(MHz^y cm)]
    alpha_power=1.5,
    alpha_mode=:stokes,
)

# Time
make_time!(kgrid, 1500.0)
dt = kgrid.dt[]
Nt = kgrid.Nt[]

# Source — CW sinusoid
freq = 1e6
signal = sin.(2π * freq .* (0:Nt-1) .* dt)
signal .*= 500_000.0  # 500 kPa amplitude

p_mask = falses(Nx)
p_mask[10] = true
source = KWaveSource(p_mask=p_mask, p=reshape(signal, 1, :))

# Sensor — multiple positions along propagation axis
sensor_mask = falses(Nx)
sensor_mask[Nx÷4] = true
sensor_mask[Nx÷2] = true
sensor_mask[3*Nx÷4] = true
sensor = KWaveSensor(mask=sensor_mask, record=[:p])

# Run
result = kspace_first_order(kgrid, medium, source, sensor)

println("Nonlinear propagation complete!")
println("  Recorded at 3 positions along the axis")
