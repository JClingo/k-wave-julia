# ============================================================================
# KWave.jl — Python interop utilities via PythonCall/juliacall
# ============================================================================
# Provides helper functions for Python users accessing KWave.jl
# via juliacall. These utilities handle numpy array conversion and
# provide a Python-friendly API.
# ============================================================================

"""
    python_kwave_grid(nx, dx, ny=nothing, dy=nothing, nz=nothing, dz=nothing)

Python-friendly grid constructor. Dispatches to the appropriate KWaveGrid
based on which dimensions are provided.

# Usage from Python:
```python
from juliacall import Main as jl
jl.seval("using KWave")
grid = jl.KWave.python_kwave_grid(128, 1e-4, 128, 1e-4)
```
"""
function python_kwave_grid(nx::Int, dx::Real,
                           ny::Union{Nothing, Int}=nothing,
                           dy::Union{Nothing, Real}=nothing,
                           nz::Union{Nothing, Int}=nothing,
                           dz::Union{Nothing, Real}=nothing)
    if ny === nothing
        return KWaveGrid(nx, Float64(dx))
    elseif nz === nothing
        return KWaveGrid(nx, Float64(dx), ny, Float64(dy))
    else
        return KWaveGrid(nx, Float64(dx), ny, Float64(dy), nz, Float64(dz))
    end
end

"""
    python_medium(; sound_speed, density=1.0, kwargs...)

Python-friendly medium constructor. Accepts plain numbers and numpy arrays
(which juliacall auto-converts to Julia arrays).
"""
function python_medium(; sound_speed, density=1.0,
                       alpha_coeff=nothing, alpha_power=nothing,
                       alpha_mode=:no_absorption, BonA=nothing)
    return KWaveMedium(;
        sound_speed=sound_speed,
        density=density,
        alpha_coeff=alpha_coeff,
        alpha_power=alpha_power,
        alpha_mode=alpha_mode,
        BonA=BonA,
    )
end

"""
    python_source(; kwargs...)

Python-friendly source constructor.
"""
function python_source(; p0=nothing, p_mask=nothing, p=nothing,
                       p_mode=Additive, u_mask=nothing,
                       ux=nothing, uy=nothing, uz=nothing,
                       u_mode=Additive)
    return KWaveSource(;
        p0=p0, p_mask=p_mask, p=p, p_mode=p_mode,
        u_mask=u_mask, ux=ux, uy=uy, uz=uz, u_mode=u_mode,
    )
end

"""
    python_sensor(; mask=nothing, record=[:p], kwargs...)

Python-friendly sensor constructor.
"""
function python_sensor(; mask=nothing, record=[:p],
                       time_reversal_boundary_data=nothing)
    # Convert Python list of strings to Vector{Symbol} if needed
    rec = record isa Vector{Symbol} ? record : Symbol.(record)
    return KWaveSensor(;
        mask=mask,
        record=rec,
        time_reversal_boundary_data=time_reversal_boundary_data,
    )
end

"""
    python_run_simulation(grid, medium, source, sensor; kwargs...)

Python-friendly simulation runner. Returns results as a Dict
that juliacall maps to a Python dict.
"""
function python_run_simulation(grid, medium, source, sensor;
                               pml_size=20, pml_alpha=2.0,
                               smooth_p0=true, data_cast=Float64)
    result = kspace_first_order(grid, medium, source, sensor;
                                pml_size=pml_size,
                                pml_alpha=pml_alpha,
                                smooth_p0=smooth_p0,
                                data_cast=data_cast)

    # Convert to plain Dict for Python friendliness
    return Dict(string(k) => v for (k, v) in result.data)
end

"""
    python_output_to_dict(output::SimulationOutput)

Convert SimulationOutput to a Dict{String, Array} for Python consumption.
"""
function python_output_to_dict(output::SimulationOutput)
    return Dict(string(k) => v for (k, v) in output.data)
end
