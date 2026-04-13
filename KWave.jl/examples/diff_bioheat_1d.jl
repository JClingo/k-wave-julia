# ============================================================================
# Example: Bioheat Transfer in 1D
# ============================================================================
# Simulates temperature evolution in tissue due to an external heat
# source (e.g., from acoustic absorption in HIFU).
# Matches MATLAB k-Wave example: example_diff_homogeneous_medium_diffusion.m

using KWave

# Grid — 1D tissue domain
Nx = 256; dx = 0.5e-3
kgrid = KWaveGrid(Nx, dx)

# Set time step for diffusion (larger than acoustic dt)
# For thermal diffusion, dt can be much larger
dt_thermal = 0.01  # 10 ms time step
kgrid.dt[] = dt_thermal
kgrid.Nt[] = 1000  # 10 seconds total

# Thermal medium — soft tissue
medium = ThermalMedium(
    thermal_conductivity=0.5,    # W/(m·K) — typical soft tissue
    density=1050.0,              # kg/m³
    specific_heat=3600.0,        # J/(kg·K)
    perfusion_rate=0.5,          # kg/(m³·s) — blood perfusion
    blood_temperature=37.0,      # °C
    blood_specific_heat=3617.0,  # J/(kg·K)
)

# Heat source — localized Gaussian heating (e.g., HIFU focus)
Q = zeros(Nx)
x_center = Nx ÷ 2
for i in 1:Nx
    Q[i] = 1e6 * exp(-((i - x_center) * dx)^2 / (2 * (2e-3)^2))
end

source = ThermalSource(Q=Q)

# Initial temperature — body temperature
T0 = 37.0

# Run thermal simulation
println("Running bioheat simulation...")
T_final, T_history = kwave_diffusion(kgrid, medium, source, T0, 1000; record_every=100)

println("Bioheat simulation complete!")
println("  Peak temperature: ", maximum(T_final), " °C")
println("  Temperature rise: ", maximum(T_final) - 37.0, " °C")
println("  History snapshots: ", size(T_history, 2))

# Compare with exact solution (no diffusion, just source + perfusion)
T_exact = bioheat_exact(37.0, maximum(Q), medium, 10.0)
println("  Exact peak (uniform): ", T_exact, " °C")
