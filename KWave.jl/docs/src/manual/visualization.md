# Visualization

## Loading a Backend

Visualization requires one of the Makie backend extensions:

```julia
using KWave, GLMakie      # desktop — interactive OpenGL window
using KWave, CairoMakie   # static PNG/SVG/PDF output
using KWave, WGLMakie     # Jupyter notebook (requires Bonito)
```

## Real-Time Simulation Display

Display the pressure field while the simulation runs:

```julia
output = kspace_first_order(kgrid, medium, source, sensor; plot_sim=true)
```

For manual control (e.g. embedding in a custom layout):

```julia
disp = create_sim_display(kgrid; layout=:default)
# inside a simulation loop:
update_sim_display!(disp, p_field, t_index)
# when done:
close_sim_display!(disp)
```

`NullDisplay` is a no-op display used when `plot_sim=false`.

## Movie Recording

```julia
rec = create_movie_recorder("simulation.mp4"; fps=30)
# inside a simulation loop:
record_frame!(rec, p_field)
# when done:
finalize_movie!(rec)
```

Pass `record_movie="simulation.mp4"` to `kspace_first_order` to record automatically.

## 2D Static Plots

```julia
using KWave, CairoMakie

# Overlay source mask on pressure field
overlay_plot(p_final, source.p_mask)

# Stack multiple line sensor signals
stacked_plot(output[:p])

# Beam pattern (log scale)
beam_plot(output[:p_max])
```

## 3D Visualization

```julia
using KWave, GLMakie

# Volume rendering
voxel_plot(p_volume)

# Isosurface at threshold
isosurface_plot(p_volume; level=0.5)

# Maximum intensity projection
max_intensity_projection(p_volume)           # plot (requires Makie)
mip = max_intensity_projection(p_volume; dims=:all)  # returns named tuple
mip.xy  # xy-plane projection
mip.xz  # xz-plane projection
```

## Colormap

`get_color_map()` returns the k-Wave blue-white-red colormap:

```julia
using GLMakie
fig, ax, hm = heatmap(p_field; colormap=get_color_map())
```

## See Also

[`get_color_map`](@ref), [`create_sim_display`](@ref), [`create_movie_recorder`](@ref),
[`voxel_plot`](@ref), [`max_intensity_projection`](@ref)
