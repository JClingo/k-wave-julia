# Python Interop

KWave.jl can construct Python-compatible k-Wave objects for use with the Python k-Wave package
via [`juliacall`](https://github.com/JuliaPy/PythonCall.jl).
This allows running the Python k-Wave simulation engine while building the simulation setup in Julia.

## Functions

```@docs
python_kwave_grid
python_medium
python_source
python_sensor
python_run_simulation
python_output_to_dict
```

## Workflow

```julia
using KWave

# Build Julia-side objects
kgrid  = KWaveGrid(128, 0.1e-3, 128, 0.1e-3)
make_time!(kgrid, 1500.0)
medium = KWaveMedium(sound_speed=1500.0)
source = KWaveSource(p0=rand(128, 128))
sensor = KWaveSensor()

# Convert to Python k-Wave objects
py_grid   = python_kwave_grid(kgrid)
py_medium = python_medium(medium)
py_source = python_source(source)
py_sensor = python_sensor(sensor)

# Run via Python k-Wave engine
output = python_run_simulation(py_grid, py_medium, py_source, py_sensor)

# Convert output back to Julia Dict
result = python_output_to_dict(output)
```
