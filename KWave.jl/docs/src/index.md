# KWave.jl

**KWave.jl** is a Julia port of the [k-Wave MATLAB toolbox](http://www.k-wave.org/) for time-domain acoustic and elastic wave simulation using k-space pseudospectral methods.

## Features

| Capability | Details |
|---|---|
| Dimensionality | 1D, 2D, 3D, axisymmetric (2D-AS) |
| Wave types | Acoustic (linear/nonlinear), elastic (P+S waves), thermal diffusion |
| Absorption | Power-law, Stokes, no-dispersion modes |
| GPU backends | CUDA, Metal, AMDGPU via package extensions |
| Visualization | GLMakie, CairoMakie, WGLMakie (Jupyter) |
| Units | Unitful.jl integration via extension |
| I/O | HDF5 (compatible with k-Wave C++ binary solver) |

## Quick Start

```julia
using KWave

# 2D grid: 256×256 points, 0.1 mm spacing
kgrid = KWaveGrid(256, 0.1e-3, 256, 0.1e-3)
make_time!(kgrid, 1500.0)

# Homogeneous water medium
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

# Gaussian initial pressure at grid centre
p0 = zeros(256, 256)
p0[128, 128] = 1.0
source = KWaveSource(p0=p0)

# Record pressure at all grid points on the boundary
sensor = KWaveSensor()

# Run simulation
output = kspace_first_order(kgrid, medium, source, sensor)
```

## Navigation

- **[Getting Started](getting_started.md)** — installation, first simulation, common pitfalls
- **[Theory](theory.md)** — k-space pseudospectral method background
- **[Manual](manual/grid.md)** — in-depth guides for each subsystem
- **[Examples](examples/index.md)** — 17 annotated example scripts
- **[API Reference](api/index.md)** — complete function and type documentation

## Citing

If you use KWave.jl, please cite the original k-Wave paper:

> B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave fields," *Journal of Biomedical Optics*, 15(2), 021314, 2010.
