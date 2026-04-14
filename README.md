# KWave.jl

A Julia port of the [k-Wave](http://www.k-wave.org/) acoustic simulation toolbox, providing time-domain simulation of acoustic wave propagation using the k-space pseudospectral method.

KWave.jl supports 1D, 2D, and 3D simulations of linear and nonlinear propagation in heterogeneous media with power-law absorption. It is a clean-room reimplementation of the [MATLAB k-Wave toolbox](https://github.com/ucl-bug/k-wave) (v1.4.1).

## Features

- **Fluid solvers** — k-space corrected pseudospectral time-domain (PSTD) solver in 1D, 2D, and 3D
- **Elastic solvers** — 2D/3D viscoelastic wave propagation with coupled compressional and shear wave fields
- **Axisymmetric solver** — efficient axially symmetric simulations
- **Continuous wave** — acoustic field propagator and angular spectrum methods
- **Thermal simulation** — bioheat diffusion solver with analytical references
- **Absorption** — power-law acoustic absorption and dispersion models
- **PML boundaries** — perfectly matched layer absorbing boundary conditions
- **Transducer arrays** — arc, bowl, disc, rect, and sphere element types
- **Reconstruction** — FFT-based line/plane reconstruction and delay-and-sum beamforming
- **Geometry** — shape primitives (disc, circle, ball, sphere, arc, bowl, line)
- **Signal processing** — tone burst, Gaussian pulse, filtering, envelope detection, noise injection
- **Material library** — water properties (sound speed, density, absorption, nonlinearity)
- **Visualization** — field display, beam plots, fly-through, voxel/isosurface rendering (via Makie extensions)
- **GPU acceleration** — CUDA, AMD ROCm, and Apple Metal via package extensions
- **I/O** — HDF5 read/write for simulation data
- **Python interop** — utilities for exchanging data with k-wave-python

## Requirements

- Julia >= 1.12

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/JClingo/k-wave-julia.git", subdir="KWave.jl")
```

Or in the Julia REPL package mode:

```
] add https://github.com/JClingo/k-wave-julia.git:KWave.jl
```

## Quick Start

```julia
using KWave

# Create a 2D grid
grid = KWaveGrid(128, 0.1e-3, 128, 0.1e-3)

# Define the medium
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

# Set up a source
source = KWaveSource(p_mask=make_disc(grid, [64, 64], 5),
                     p=tone_burst(1/grid.dt, 1e6, 5))

# Configure sensor
sensor = KWaveSensor(mask=make_circle(grid, [64, 64], 40))

# Run simulation
output = kspace_first_order(grid, medium, source, sensor)
```

## Visualization

KWave.jl uses [Makie](https://docs.makie.org/) for visualization via package extensions. Load a Makie backend alongside KWave to enable plotting. No backend is required to run simulations — visualization is fully optional.

### Choosing a Backend

There are three Makie backends, each suited to different workflows:

| Backend | Interactive | GPU | File Export | Best For |
|---------|:-----------:|:---:|:-----------:|----------|
| **GLMakie** | ✅ | ✅ | Screenshot | Real-time exploration, 3D, animations |
| **CairoMakie** | ❌ | ❌ | PNG / PDF / SVG | Publication figures, headless/CI |
| **WGLMakie** | ✅ | ❌ | Browser | Jupyter notebooks, web apps |

Install whichever backend(s) you need:

```julia
using Pkg
Pkg.add("GLMakie")      # interactive native windows
Pkg.add("CairoMakie")   # static publication-quality output
Pkg.add("WGLMakie")     # browser / Jupyter
```

Only one backend should be loaded per Julia session. Start a new session to switch backends.

#### GLMakie — Interactive Windows

Requires OpenGL drivers (works on Linux, macOS, and Windows with a GPU or integrated graphics).

```julia
using KWave, GLMakie

output = kspace_first_order(grid, medium, source, sensor)

beam_plot(output.p_final)       # Opens an interactive window
fly_through(output_3d.p_final)  # Slider-controlled 3D slicer
voxel_plot(output_3d.p_max)     # Real-time 3D volume rendering
```

#### CairoMakie — Static Figures & File Export

No GPU required. Best for headless servers, CI pipelines, and saving publication-quality figures.

```julia
using KWave, CairoMakie

output = kspace_first_order(grid, medium, source, sensor)

fig = beam_plot(output.p_final; db_scale=true)
save("beam.png", fig)   # PNG
save("beam.pdf", fig)   # Vector PDF (for papers)
save("beam.svg", fig)   # Scalable SVG
```

#### WGLMakie — Jupyter / Web

Renders interactive plots directly inside Jupyter notebooks or web applications.

```julia
using KWave, Bonito, WGLMakie   # run inside a Jupyter notebook

Page()
WGLMakie.activate!()

output = kspace_first_order(grid, medium, source, sensor)
beam_plot(output.p_final)   # Renders inline, with pan/zoom
```

---

### 2D Beam Pattern

```julia
using KWave, CairoMakie  # or GLMakie for interactive windows

# Run a 2D simulation...
output = kspace_first_order(grid, medium, source, sensor)

# Plot the final pressure field as a beam pattern
beam_plot(output.p_final; db_scale=true, db_range=40)
```

### Overlay Plot (Pressure on Medium)

```julia
using KWave, CairoMakie

# Overlay pressure field on the sound speed map
overlay_plot(medium.sound_speed, output.p_max;
             alpha=0.6,
             background_cmap=:grays)
```

### Stacked Sensor Signals

```julia
using KWave, CairoMakie

# Display time series recorded at multiple sensor points
stacked_plot(output.p; dt=grid.dt, labels=["Sensor $i" for i in 1:size(output.p, 1)])
```

### Interactive 3D Fly-Through

```julia
using KWave, GLMakie  # GLMakie required for interactive sliders

output_3d = kspace_first_order(grid3d, medium, source, sensor)

# Slider-controlled slice viewer through a 3D volume
fly_through(output_3d.p_final; dim=3)
```

### 3D Voxel and Isosurface Rendering

```julia
using KWave, GLMakie

# Volume rendering
voxel_plot(output_3d.p_max; mode=:volume, db_scale=true, db_range=30)

# Isosurface at 50% of peak pressure
isosurface_plot(output_3d.p_max, 0.5 * maximum(output_3d.p_max); alpha=0.8)

# Maximum intensity projection
max_intensity_projection(output_3d.p_max)
```

### Real-Time Simulation Display

```julia
using KWave, GLMakie

# Pass plot_sim=true to watch the simulation evolve in real time
output = kspace_first_order(grid, medium, source, sensor; plot_sim=true)
```

### Saving Figures

All Makie plots can be saved to file:

```julia
using CairoMakie

fig = beam_plot(output.p_final; db_scale=true)
save("beam_pattern.png", fig)   # PNG, PDF, SVG supported with CairoMakie
```

## GPU Acceleration

Load the appropriate GPU extension for your hardware:

```julia
using KWave
using CUDA        # NVIDIA GPUs
# using AMDGPU    # AMD GPUs
# using Metal     # Apple Silicon
```

## Examples

See the [`examples/`](KWave.jl/examples/) directory for complete simulation scripts:

| Example | Description |
|---------|-------------|
| `ivp_homogeneous_medium_2d.jl` | Basic 2D acoustic propagation |
| `ivp_heterogeneous_medium_2d.jl` | Heterogeneous medium simulation |
| `ivp_absorption_2d.jl` | Absorption modeling |
| `ivp_1d_simulation.jl` | 1D simulation |
| `ivp_3d_simulation.jl` | 3D simulation |
| `tvsp_homogeneous_medium_2d.jl` | Time-varying source problem |
| `tvsp_nonlinear_propagation_1d.jl` | Nonlinear propagation |
| `pr_time_reversal_2d.jl` | Photoacoustic time reversal |
| `pr_fft_reconstruction_2d.jl` | FFT-based reconstruction |
| `us_beamforming_2d.jl` | Ultrasound beamforming |
| `us_phased_array_3d.jl` | 3D phased array |
| `ewp_elastic_2d.jl` | 2D elastic wave propagation |
| `sd_sensor_directivity_2d.jl` | Sensor directivity |
| `diff_bioheat_1d.jl` | Bioheat diffusion |
| `vis_cairomakie_2d.jl` | Static figures saved to PNG/PDF/SVG (CairoMakie) |
| `vis_glmakie_interactive.jl` | Interactive 2D/3D windows, real-time display (GLMakie) |
| `vis_wglmakie_notebook.jl` | Inline Jupyter / web browser plots (WGLMakie) |

## Running Tests

```bash
cd KWave.jl
julia --project -e 'using Pkg; Pkg.test()'
```

## Project Structure

```
KWave.jl/
  src/
    KWave.jl          # Main module
    solver/            # Acoustic, elastic, thermal, CW solvers
    geometry/          # Shape primitives
    signal/            # Signal generation and processing
    filter/            # Spatial and temporal filters
    material/          # Material property functions
    reconstruction/    # Image reconstruction algorithms
    visualization/     # Plotting and display
    io/                # HDF5 I/O
    interop/           # Python interoperability
    reference/         # Analytical reference solutions
  ext/                 # GPU and visualization extensions
  test/                # Test suite
  examples/            # Example simulation scripts
```

## License

Apache-2.0

## Acknowledgments

Based on the [k-Wave MATLAB toolbox](http://www.k-wave.org/) developed at University College London. This is a clean-room port — no MATLAB source code was used.

## References

- B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave fields," *J. Biomed. Opt.*, vol. 15, no. 2, 021314, 2010.
