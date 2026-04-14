# ============================================================================
# Example: Visualization with WGLMakie (Jupyter / Web)
# ============================================================================
# Demonstrates interactive plots that render inline in a Jupyter notebook
# or web browser. No GPU or display server required.
#
# To run in Jupyter:
#   1. Install: ] add WGLMakie Bonito IJulia
#   2. Start:   jupyter notebook
#   3. Open this file or paste cells into a notebook
#
# To run as a standalone web app (Bonito/Pluto):
#   ] add Bonito Pluto
#
# Install: ] add WGLMakie Bonito

using KWave
using Bonito
using WGLMakie

Page()
WGLMakie.activate!()

# ----------------------------------------------------------------------------
# Simulation — 2D heterogeneous medium
# ----------------------------------------------------------------------------
Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

# Heterogeneous sound speed — soft tissue over water
c_map = fill(1500.0, Nx, Ny)
c_map[Nx÷2:end, :] .= 1800.0      # faster upper half (e.g. bone layer)

rho_map = fill(1000.0, Nx, Ny)
rho_map[Nx÷2:end, :] .= 1200.0

medium = KWaveMedium(sound_speed=c_map, density=rho_map)
make_time!(kgrid, maximum(c_map))

# Source — disc
p0 = zeros(Nx, Ny)
p0[make_disc(Nx, Ny, 20, Ny÷2, 6)] .= 1.0
source = KWaveSource(p0=p0)

# Sensor — full field
sensor = KWaveSensor(mask=trues(Nx, Ny), record=[:p_final, :p_max])

println("Running simulation...")
result = kspace_first_order(kgrid, medium, source, sensor)
println("Done.")

# Full-field sensor data is returned flattened over the grid.
p_final = reshape(result[:p_final], Nx, Ny)
p_max = reshape(result[:p_max], Nx, Ny)

# ----------------------------------------------------------------------------
# Cell 1: Inline beam plot
# (In Jupyter this renders as an interactive figure in the output cell)
# ----------------------------------------------------------------------------
fig1 = beam_plot(p_final; db_scale=true, db_range=40)
fig1 === nothing && error("WGLMakie plotting extension did not load. Make sure both Bonito and WGLMakie are loaded in this session.")
display(fig1)

# ----------------------------------------------------------------------------
# Cell 2: Overlay — pressure field on the sound speed map
# Pan and zoom work directly in the notebook output cell
# ----------------------------------------------------------------------------
fig2 = overlay_plot(c_map, p_max;
                    alpha=0.65,
                    background_cmap=:thermal)
display(fig2)

# ----------------------------------------------------------------------------
# Cell 3: Stacked sensor signals
# Record along a vertical line and display all channels
# ----------------------------------------------------------------------------
line_mask = falses(Nx, Ny)
line_mask[Nx÷4, :] .= true
sensor_line = KWaveSensor(mask=line_mask, record=[:p])

println("Running second pass for sensor time series...")
result_line = kspace_first_order(kgrid, medium, source, sensor_line)

fig3 = stacked_plot(result_line[:p];
                    dt=kgrid.dt[],
                    spacing=0.25)
display(fig3)

# ----------------------------------------------------------------------------
# Note: 3D visualization
# WGLMakie supports 3D plots rendered via WebGL in the browser.
# Run a 3D simulation (see ivp_3d_simulation.jl) then:
#
#   using WGLMakie
#   voxel_plot(result3d[:p_max]; mode=:volume, db_scale=true)
#   fly_through(result3d[:p_final]; dim=3)
#
# These render as interactive WebGL canvases in the notebook output cell.
# ----------------------------------------------------------------------------

println("Notebook visualization complete.")
