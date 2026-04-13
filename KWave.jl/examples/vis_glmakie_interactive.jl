# ============================================================================
# Example: Interactive Visualization with GLMakie
# ============================================================================
# Demonstrates interactive 2D/3D visualization in native windows.
# GLMakie provides GPU-accelerated rendering with mouse pan/zoom/rotate
# and slider controls.
#
# Requires: OpenGL-capable GPU or integrated graphics (Linux/macOS/Windows)
# Install:  ] add GLMakie

using KWave
using GLMakie

# ============================================================================
# Part 1: 2D simulation — interactive beam plot
# ============================================================================
println("=== Part 1: 2D interactive beam plot ===")

Nx = 128; dx = 0.1e-3
Ny = 128; dy = 0.1e-3
kgrid = KWaveGrid(Nx, dx, Ny, dy)

medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
make_time!(kgrid, 1500.0)

p0 = zeros(Nx, Ny)
p0[make_disc(Nx, Ny, Nx÷2, Ny÷2, 8)] .= 1.0
source = KWaveSource(p0=p0)

sensor_mask = trues(Nx, Ny)
sensor = KWaveSensor(mask=sensor_mask, record=[:p_final, :p_max])

println("Running 2D simulation...")
result2d = kspace_first_order(kgrid, medium, source, sensor)

# Interactive beam pattern — opens a native window with pan/zoom
println("Opening beam plot (close window to continue)...")
fig1 = beam_plot(result2d[:p_final]; db_scale=true, db_range=40)
display(fig1)

# Overlay — pressure on the medium map
fig2 = overlay_plot(medium.sound_speed, result2d[:p_max]; alpha=0.6)
display(fig2)

# ============================================================================
# Part 2: 3D simulation — fly-through and voxel rendering
# ============================================================================
println("\n=== Part 2: 3D interactive visualization ===")

Nx3 = 64; dx3 = 0.1e-3
Ny3 = 64; dy3 = 0.1e-3
Nz3 = 64; dz3 = 0.1e-3
kgrid3d = KWaveGrid(Nx3, dx3, Ny3, dy3, Nz3, dz3)

medium3d = KWaveMedium(sound_speed=1500.0, density=1000.0)
make_time!(kgrid3d, 1500.0)

p0_3d = zeros(Nx3, Ny3, Nz3)
p0_3d[make_ball(Nx3, Ny3, Nz3, Nx3÷2, Ny3÷2, Nz3÷2, 5)] .= 1.0
source3d = KWaveSource(p0=p0_3d)

sensor3d = KWaveSensor(mask=trues(Nx3, Ny3, Nz3), record=[:p_final, :p_max])

println("Running 3D simulation (this may take a moment)...")
result3d = kspace_first_order(kgrid3d, medium3d, source3d, sensor3d)

# Slider-controlled slice viewer — drag the slider to move through Z planes
println("Opening fly-through viewer (drag slider to slice)...")
fig3 = fly_through(result3d[:p_final]; dim=3, db_scale=false)
display(fig3)

# Volume rendering — interactive 3D view with mouse rotation
println("Opening voxel volume render...")
fig4 = voxel_plot(result3d[:p_max];
                  mode=:volume,
                  db_scale=true, db_range=30,
                  dx=dx3, dy=dy3, dz=dz3)
display(fig4)

# Isosurface at 50% of peak
println("Opening isosurface plot...")
threshold = 0.5 * maximum(result3d[:p_max])
fig5 = isosurface_plot(result3d[:p_max], threshold;
                       alpha=0.8,
                       dx=dx3, dy=dy3, dz=dz3)
display(fig5)

# Maximum intensity projection
println("Opening max-intensity projection...")
fig6 = max_intensity_projection(result3d[:p_max])
display(fig6)

# ============================================================================
# Part 3: Real-time simulation monitoring
# ============================================================================
println("\n=== Part 3: Real-time simulation display ===")

kgrid_rt = KWaveGrid(Nx, dx, Ny, dy)
make_time!(kgrid_rt, 1500.0)

p0_rt = zeros(Nx, Ny)
p0_rt[make_disc(Nx, Ny, 40, 64, 5)] .= 1.0
source_rt = KWaveSource(p0=p0_rt)
sensor_rt = KWaveSensor(mask=trues(Nx, Ny), record=[:p_max])

# plot_sim=true opens a live window that updates each time step
println("Running with real-time display (watch the wavefront evolve)...")
result_rt = kspace_first_order(kgrid_rt, medium, source_rt, sensor_rt;
                               plot_sim=true)

println("\nDone. Close any open windows to exit.")
