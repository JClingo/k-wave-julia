# ============================================================================
# Example: Visualization with CairoMakie (Static / File Export)
# ============================================================================
# Demonstrates saving publication-quality figures from a 2D simulation.
# CairoMakie renders to PNG, PDF, and SVG without a display or GPU.
#
# Install: ] add CairoMakie

using KWave
using CairoMakie

require_figure(name::AbstractString, fig) =
    fig === nothing ? error("$name returned nothing. Activate the local KWave project with `julia --project=KWave.jl` or `Pkg.activate(\"KWave.jl\"); Pkg.develop(path=\"KWave.jl\")`, then load CairoMakie again so the extension methods are available.") : fig

const OUTPUT_DIR = joinpath(@__DIR__, "generated")
mkpath(OUTPUT_DIR)

output_path(name::AbstractString) = joinpath(OUTPUT_DIR, name)

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

# Record the full time series as well as summary fields used below.
sensor_mask = trues(Nx, Ny)
sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_final, :p_max])

println("Running simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)
println("Done. Saving figures...")

# Full-field sensor data is returned flattened over the grid.
p_final = reshape(result[:p_final], Nx, Ny)
p_max = reshape(result[:p_max], Nx, Ny)
p = reshape(result[:p], Nx, Ny, :)

# ----------------------------------------------------------------------------
# Figure 1: Beam pattern of the final pressure field (linear scale)
# ----------------------------------------------------------------------------
fig1 = require_figure("beam_plot", beam_plot(p_final; db_scale=false))
save(output_path("cairomakie_beam_linear.png"), fig1)
println("  Saved: $(output_path("cairomakie_beam_linear.png"))")

# ----------------------------------------------------------------------------
# Figure 2: Beam pattern in dB (40 dB dynamic range)
# ----------------------------------------------------------------------------
fig2 = require_figure("beam_plot", beam_plot(p_max; db_scale=true, db_range=40))
save(output_path("cairomakie_beam_db.png"), fig2)
save(output_path("cairomakie_beam_db.pdf"), fig2)   # vector PDF for papers
save(output_path("cairomakie_beam_db.svg"), fig2)   # scalable SVG for web
println("  Saved: $(output_path("cairomakie_beam_db.png")) / .pdf / .svg")

# ----------------------------------------------------------------------------
# Figure 3: Overlay — peak pressure on top of the sound speed map
# ----------------------------------------------------------------------------
fig3 = require_figure(
    "overlay_plot",
    overlay_plot(fill(medium.sound_speed, Nx, Ny), p_max;
                 alpha=0.6,
                 background_cmap=:grays),
)
save(output_path("cairomakie_overlay.png"), fig3)
println("  Saved: $(output_path("cairomakie_overlay.png"))")

# ----------------------------------------------------------------------------
# Figure 4: Stacked sensor signals along a horizontal line
# ----------------------------------------------------------------------------
# Extract signals from the top row of the grid.
top_row_signals = p[1, :, :]

fig4 = require_figure(
    "stacked_plot",
    stacked_plot(top_row_signals;
                 dt=kgrid.dt[],
                 spacing=0.3),
)
save(output_path("cairomakie_stacked_signals.png"), fig4)
println("  Saved: $(output_path("cairomakie_stacked_signals.png"))")

println("\nAll figures saved.")
