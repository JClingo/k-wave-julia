# Getting Started

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/JClingo/k-wave-julia.git", subdir="KWave.jl")
```

Or from a Julia REPL:

```
pkg> add https://github.com/JClingo/k-wave-julia.git:KWave.jl
```

### Optional Dependencies

| Package | Purpose | Load pattern |
|---|---|---|
| `GLMakie` | Interactive visualization | `using KWave, GLMakie` |
| `CairoMakie` | Static / publication plots | `using KWave, CairoMakie` |
| `WGLMakie` + `Bonito` | Jupyter notebook visualization | `using KWave, WGLMakie, Bonito` |
| `CUDA` | NVIDIA GPU acceleration | `using KWave, CUDA` |
| `Metal` | Apple Silicon GPU acceleration | `using KWave, Metal` |
| `AMDGPU` | AMD GPU acceleration | `using KWave, AMDGPU` |
| `Unitful` | Physical units in medium/grid | `using KWave, Unitful` |

## Your First Simulation

The example below simulates a 2D initial value problem (photoacoustic source) propagating in a homogeneous medium.

### 1. Create the Grid

```julia
using KWave

Nx, Ny = 256, 256          # grid points
dx, dy = 0.1e-3, 0.1e-3   # grid spacing [m]
kgrid = KWaveGrid(Nx, dx, Ny, dy)
```

`KWaveGrid` pre-computes wavenumber vectors, staggered shift operators, and spatial coordinate arrays.

### 2. Set the Time Step

```julia
c_max = 1500.0             # maximum sound speed [m/s]
make_time!(kgrid, c_max)   # sets kgrid.dt and kgrid.Nt using CFL = 0.3
```

`make_time!` must be called before running any solver. It modifies `kgrid.dt[]`, `kgrid.Nt[]`, and `kgrid.t_array` in-place.

### 3. Define the Medium

```julia
medium = KWaveMedium(
    sound_speed = 1500.0,   # [m/s]
    density     = 1000.0,   # [kg/m³]
)
```

For absorbing media, add `alpha_coeff` and `alpha_power`:

```julia
medium = KWaveMedium(
    sound_speed = 1500.0,
    density     = 1000.0,
    alpha_coeff = 0.75,     # [dB/(MHz^y cm)]
    alpha_power = 1.5,
)
```

### 4. Set Up the Source

```julia
p0 = zeros(Nx, Ny)
p0[128, 128] = 1.0                 # point source at grid centre
source = KWaveSource(p0=p0)
```

### 5. Configure the Sensor

```julia
# Record pressure time series at all grid points on all four edges
mask = falses(Nx, Ny)
mask[1, :]   .= true
mask[end, :] .= true
mask[:, 1]   .= true
mask[:, end] .= true

sensor = KWaveSensor(mask=mask, record=[:p])
```

### 6. Run the Simulation

```julia
output = kspace_first_order(kgrid, medium, source, sensor)
```

### 7. Access Results

```julia
p_data = output[:p]    # Matrix (n_sensor_points × Nt)
```

## Common Pitfalls

| Problem | Cause | Fix |
|---|---|---|
| `UndefRefError: access to undefined reference` | `make_time!` not called | Call `make_time!(kgrid, c_max)` before solver |
| PML artefacts at boundaries | PML too thin | Increase `pml_size` (default 20); try 40 for large models |
| Memory error on large 3D grid | Full field sensor | Set `sensor.mask` to record only points of interest |
| Slow first run | JIT compilation | Subsequent runs of same grid size are fast |
| GPU simulation gives wrong results | `Float64` on GPU | Set `data_cast=Float32` when using GPU backends |
