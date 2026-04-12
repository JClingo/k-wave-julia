# ============================================================================
# KWave.jl — KWaveSensor and SimulationOutput
# ============================================================================

"""
    KWaveSensor

Sensor/detector definitions for k-Wave simulations.

# Fields
- `mask`: Binary mask (record at grid points where true), Cartesian coordinates,
  or `nothing` (record entire field)
- `record`: Vector of symbols specifying which fields to record:
  `:p`, `:p_max`, `:p_min`, `:p_rms`, `:p_final`,
  `:ux`, `:uy`, `:uz`, `:u_max`, `:u_rms`, `:u_final`,
  `:I_avg`, `:I_max`
- `time_reversal_boundary_data`: Sensor data for time-reversal reconstruction
"""
Base.@kwdef struct KWaveSensor
    mask::Union{Nothing, AbstractArray{Bool}, AbstractMatrix{Float64}} = nothing
    record::Vector{Symbol} = [:p]
    time_reversal_boundary_data::Union{Nothing, AbstractArray} = nothing
    directivity_angle::Union{Nothing, AbstractArray} = nothing
    directivity_size::Union{Nothing, Float64} = nothing
    frequency_response::Union{Nothing, Tuple{Float64, Float64}} = nothing
end

# Valid record field names
const VALID_RECORD_FIELDS = Set([
    :p, :p_max, :p_min, :p_rms, :p_final,
    :ux, :uy, :uz, :u_max, :u_rms, :u_final,
    :I_avg, :I_max,
])

"""
    validate_record_fields(sensor)

Check that all requested record fields are valid.
"""
function validate_record_fields(sensor::KWaveSensor)
    for field in sensor.record
        if field ∉ VALID_RECORD_FIELDS
            error("Invalid sensor record field: :$field. Valid fields: $(join(sort(collect(VALID_RECORD_FIELDS)), ", "))")
        end
    end
end

"""
    is_time_reversal(sensor)

Check if the sensor is configured for time-reversal reconstruction.
"""
is_time_reversal(s::KWaveSensor) = s.time_reversal_boundary_data !== nothing

# ============================================================================
# SimulationOutput
# ============================================================================

"""
    SimulationOutput

Container for simulation results. Access fields via `output[:p]`, `output[:p_max]`, etc.
"""
struct SimulationOutput
    data::Dict{Symbol, AbstractArray}
end

Base.getindex(o::SimulationOutput, key::Symbol) = o.data[key]
Base.haskey(o::SimulationOutput, key::Symbol) = haskey(o.data, key)
Base.keys(o::SimulationOutput) = keys(o.data)

function Base.show(io::IO, o::SimulationOutput)
    fields = sort(collect(keys(o.data)))
    print(io, "SimulationOutput with fields: ", join(fields, ", "))
end

"""
    _create_sensor_data(sensor, kgrid, Nt)

Pre-allocate storage arrays for sensor data based on requested record fields.
"""
function _create_sensor_data(sensor::KWaveSensor, grid_size::Tuple, Nt::Int)
    data = Dict{Symbol, AbstractArray}()

    if sensor.mask === nothing
        return data
    end

    # Determine number of sensor points
    if sensor.mask isa AbstractArray{Bool}
        n_sensor = count(sensor.mask)
    else
        # Cartesian coordinates: columns are points
        n_sensor = size(sensor.mask, 2)
    end

    for field in sensor.record
        if field == :p
            data[:p] = zeros(Float64, n_sensor, Nt)
        elseif field == :p_max
            data[:p_max] = fill(-Inf, n_sensor)
        elseif field == :p_min
            data[:p_min] = fill(Inf, n_sensor)
        elseif field == :p_rms
            data[:p_rms] = zeros(Float64, n_sensor)
        elseif field == :p_final
            data[:p_final] = zeros(Float64, n_sensor)
        elseif field in (:ux, :uy, :uz)
            data[field] = zeros(Float64, n_sensor, Nt)
        elseif field == :u_max
            data[:u_max] = fill(-Inf, n_sensor)
        elseif field == :u_rms
            data[:u_rms] = zeros(Float64, n_sensor)
        elseif field == :u_final
            data[:u_final] = zeros(Float64, n_sensor)
        elseif field == :I_avg
            data[:I_avg] = zeros(Float64, n_sensor)
        elseif field == :I_max
            data[:I_max] = fill(-Inf, n_sensor)
        end
    end

    return data
end
