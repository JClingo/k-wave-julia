# ============================================================================
# Example: Visualization with CairoMakie (Static / File Export)
# ============================================================================
# Demonstrates saving publication-quality figures from a 2D simulation.
# CairoMakie renders to PNG, PDF, and SVG without a display or GPU.
#
# Install: ] add CairoMakie

using KWave
using CairoMakie

# ----------------------------------------------------------------------------
# Simulation setup — 2D propagation with absorption
# ----------------------------------------------------------------------------
Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

medium = KWaveMedium(
    sound_speed=1500.0,
    density=1000.0,
    alpha_coeff=0.5,
    alpha_power=1.5,
)

make_time!(kgrid, 1500.0)

# Initial pressure — two discs
p0 = zeros(Nx, Ny)
p0[make_disc(Nx, Ny, 44, 64, 6)] .=  1.0
p0[make_disc(Nx, Ny, 84, 64, 6)] .= -1.0
source = KWaveSource(p0=p0)

# Record pressure and peak pressure at all grid points
sensor_mask = trues(Nx, Ny)
sensor = KWaveSensor(mask=sensor_mask, record=[:p_final, :p_max])

println("Running simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)
println("Done. Saving figures...")

# ----------------------------------------------------------------------------
# Figure 1: Beam pattern of the final pressure field (linear scale)
# ----------------------------------------------------------------------------
fig1 = beam_plot(result[:p_final]; db_scale=false)
save("cairomakie_beam_linear.png", fig1)
println("  Saved: cairomakie_beam_linear.png")

# ----------------------------------------------------------------------------
# Figure 2: Beam pattern in dB (40 dB dynamic range)
# ----------------------------------------------------------------------------
fig2 = beam_plot(result[:p_max]; db_scale=true, db_range=40)
save("cairomakie_beam_db.png", fig2)
save("cairomakie_beam_db.pdf", fig2)   # vector PDF for papers
save("cairomakie_beam_db.svg", fig2)   # scalable SVG for web
println("  Saved: cairomakie_beam_db.png / .pdf / .svg")

# ----------------------------------------------------------------------------
# Figure 3: Overlay — peak pressure on top of the sound speed map
# ----------------------------------------------------------------------------
fig3 = overlay_plot(medium.sound_speed, result[:p_max];
                    alpha=0.6,
                    background_cmap=:grays)
save("cairomakie_overlay.png", fig3)
println("  Saved: cairomakie_overlay.png")

# ----------------------------------------------------------------------------
# Figure 4: Stacked sensor signals along a horizontal line
# ----------------------------------------------------------------------------
# Extract signals from the top row of the grid
top_row_signals = result[:p][1:Ny, :]   # shape: (Ny, Nt)

fig4 = stacked_plot(top_row_signals;
                    dt=kgrid.dt[],
                    spacing=0.3)
save("cairomakie_stacked_signals.png", fig4)
println("  Saved: cairomakie_stacked_signals.png")

println("\nAll figures saved.")
