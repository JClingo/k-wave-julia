# ============================================================================
# KWave.jl — KWaveSource
# ============================================================================

"""
    KWaveSource{T}

Source definitions for k-Wave simulations.

# Fields
## Initial value (photoacoustic) source
- `p0`: Initial pressure distribution

## Time-varying pressure source
- `p_mask`: Binary mask defining pressure source locations
- `p`: Pressure time series (n_sources × Nt) or (n_sources × 1) for steady
- `p_mode`: Source injection mode (default: `Additive`)

## Time-varying velocity source
- `u_mask`: Binary mask defining velocity source locations
- `ux`, `uy`, `uz`: Velocity component time series
- `u_mode`: Source injection mode (default: `Additive`)
"""
struct KWaveSource{T<:AbstractFloat}
    p0::Union{Nothing, AbstractArray{T}}
    p_mask::Union{Nothing, AbstractArray{Bool}}
    p::Union{Nothing, AbstractArray{T}}
    p_mode::SourceMode
    u_mask::Union{Nothing, AbstractArray{Bool}}
    ux::Union{Nothing, AbstractArray{T}}
    uy::Union{Nothing, AbstractArray{T}}
    uz::Union{Nothing, AbstractArray{T}}
    u_mode::SourceMode
end

"""
    KWaveSource(; p0=nothing, kwargs...)

Create a KWaveSource. All array inputs are converted to `Float64`.

# Keyword Arguments
- `p0`: Initial pressure distribution (photoacoustic / IVP source)
- `p_mask`: Binary `BitArray` mask for time-varying pressure source locations
- `p`: Pressure time series — `(n_sources × Nt)` or `(n_sources × 1)` for steady-state
- `p_mode`: Source injection mode (default: `Additive`)
- `u_mask`: Binary mask for velocity source locations
- `ux`, `uy`, `uz`: Velocity component time series — `(n_sources × Nt)`
- `u_mode`: Velocity source injection mode (default: `Additive`)

# Notes
- `p0` and time-varying sources (`p_mask`/`p`) can coexist in the same simulation.
- For a single source point, `p` should be a `(1 × Nt)` matrix (not a vector).
- Use [`KWaveArray`](@ref) with [`get_distributed_source_signal`](@ref) for multi-element arrays.

# See Also
[`KWaveSensor`](@ref), [`SourceMode`](@ref), [`tone_burst`](@ref),
[`KWaveArray`](@ref), [`kspace_first_order`](@ref)
"""
function KWaveSource(; p0::Union{Nothing, AbstractArray{<:Real}}=nothing,
                     p_mask::Union{Nothing, AbstractArray{Bool}}=nothing,
                     p::Union{Nothing, AbstractArray{<:Real}}=nothing,
                     p_mode::SourceMode=Additive,
                     u_mask::Union{Nothing, AbstractArray{Bool}}=nothing,
                     ux::Union{Nothing, AbstractArray{<:Real}}=nothing,
                     uy::Union{Nothing, AbstractArray{<:Real}}=nothing,
                     uz::Union{Nothing, AbstractArray{<:Real}}=nothing,
                     u_mode::SourceMode=Additive)
    _f64 = x -> x === nothing ? nothing : Float64.(x)
    return KWaveSource{Float64}(
        _f64(p0), p_mask, _f64(p), p_mode,
        u_mask, _f64(ux), _f64(uy), _f64(uz), u_mode
    )
end

"""
    has_p0(source)

Check if the source has an initial pressure distribution.
"""
has_p0(s::KWaveSource) = s.p0 !== nothing

"""
    has_pressure_source(source)

Check if the source has a time-varying pressure source.
"""
has_pressure_source(s::KWaveSource) = s.p_mask !== nothing && s.p !== nothing

"""
    has_velocity_source(source)

Check if the source has a time-varying velocity source.
"""
has_velocity_source(s::KWaveSource) = s.u_mask !== nothing
