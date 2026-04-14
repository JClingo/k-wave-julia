# Python Interop

KWave.jl can bridge to the [Python k-Wave package](https://github.com/waltsims/k-wave-python)
using [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl) / `juliacall`.

## Setup

Install prerequisites:

```
pkg> add PythonCall
pip install k-wave-python
```

## Workflow

Build simulation objects in Julia, convert to Python, run, convert results back:

```julia
using KWave

kgrid  = KWaveGrid(128, 0.1e-3, 128, 0.1e-3)
make_time!(kgrid, 1500.0)
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
source = KWaveSource(p0=rand(128, 128))
sensor = KWaveSensor()

py_grid   = python_kwave_grid(kgrid)
py_medium = python_medium(medium)
py_source = python_source(source)
py_sensor = python_sensor(sensor)

py_output = python_run_simulation(py_grid, py_medium, py_source, py_sensor)
result    = python_output_to_dict(py_output)

# result is a Dict{String, Array} of sensor data
p_data = result["p"]
```

## Use Case

Useful when you want to:
- Compare KWave.jl results against the Python k-Wave engine
- Use Julia for setup/post-processing and Python k-Wave's GPU engine for the simulation
- Gradually migrate an existing Python k-Wave workflow to Julia

## See Also

[`python_kwave_grid`](@ref), [`python_run_simulation`](@ref)
