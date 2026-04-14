# Visualization

Visualization functions require a Makie backend extension to be loaded.
Load one of the following before calling plot functions:

```julia
using GLMakie      # desktop interactive
using CairoMakie   # static / publication quality
using WGLMakie     # Jupyter notebook
```

## Colormap

```@docs
get_color_map
```

## Real-Time Simulation Display

```@docs
SimulationDisplay
NullDisplay
create_sim_display
update_sim_display!
close_sim_display!
```

## Movie Recording

```@docs
MovieRecorder
NullRecorder
create_movie_recorder
record_frame!
finalize_movie!
```

## Static Plots

```@docs
beam_plot
fly_through
overlay_plot
stacked_plot
```

## 3D Visualization

```@docs
voxel_plot
isosurface_plot
max_intensity_projection
```

## Notes

- `NullDisplay` and `NullRecorder` are no-op implementations used when `plot_sim=false` or `record_movie=nothing`. They allow solver internals to call display/recording functions unconditionally.
- `get_color_map` returns the k-Wave blue-white-red colormap as a vector of `RGBf` tuples, compatible with Makie's `colormap` keyword.
