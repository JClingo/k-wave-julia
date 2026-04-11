# ============================================================================
# KWave.jl — KWaveMedium
# ============================================================================

"""
    KWaveMedium{T}

Acoustic medium properties for k-Wave simulations.

# Fields
- `sound_speed`: Sound speed [m/s] — scalar (homogeneous) or array (heterogeneous)
- `density`: Mass density [kg/m³] — scalar or array (default: 1.0)
- `alpha_coeff`: Power-law absorption coefficient [dB/(MHz^y cm)] (optional)
- `alpha_power`: Power-law absorption exponent (optional)
- `alpha_mode`: Absorption mode — `:no_absorption`, `:no_dispersion`, `:stokes`
- `BonA`: Nonlinearity parameter B/A (optional)
"""
struct KWaveMedium{T<:AbstractFloat}
    sound_speed::Union{T, AbstractArray{T}}
    density::Union{T, AbstractArray{T}}
    alpha_coeff::Union{Nothing, T, AbstractArray{T}}
    alpha_power::Union{Nothing, T}
    alpha_mode::Symbol
    BonA::Union{Nothing, T, AbstractArray{T}}
end

"""
    KWaveMedium(; sound_speed, density=1.0, kwargs...)

Create a KWaveMedium with Float64 values by default.
"""
function KWaveMedium(; sound_speed::Union{Real, AbstractArray{<:Real}},
                     density::Union{Real, AbstractArray{<:Real}}=1.0,
                     alpha_coeff::Union{Nothing, Real, AbstractArray{<:Real}}=nothing,
                     alpha_power::Union{Nothing, Real}=nothing,
                     alpha_mode::Symbol=:no_absorption,
                     BonA::Union{Nothing, Real, AbstractArray{<:Real}}=nothing)
    _f64 = x -> x === nothing ? nothing : (x isa AbstractArray ? Float64.(x) : Float64(x))
    return KWaveMedium{Float64}(
        _f64(sound_speed), _f64(density),
        _f64(alpha_coeff),
        alpha_power === nothing ? nothing : Float64(alpha_power),
        alpha_mode, _f64(BonA)
    )
end

"""
    is_homogeneous(medium)

Check if the medium has uniform (scalar) sound speed and density.
"""
is_homogeneous(m::KWaveMedium) = m.sound_speed isa Real && m.density isa Real

"""
    is_lossless(medium)

Check if the medium has no absorption.
"""
is_lossless(m::KWaveMedium) = m.alpha_coeff === nothing || m.alpha_mode == :no_absorption

"""
    is_nonlinear(medium)

Check if the medium has nonlinear propagation enabled.
"""
is_nonlinear(m::KWaveMedium) = m.BonA !== nothing

"""
    max_sound_speed(medium)

Return the maximum sound speed in the medium.
"""
max_sound_speed(m::KWaveMedium) = m.sound_speed isa Real ? m.sound_speed : maximum(m.sound_speed)

"""
    min_sound_speed(medium)

Return the minimum sound speed in the medium.
"""
min_sound_speed(m::KWaveMedium) = m.sound_speed isa Real ? m.sound_speed : minimum(m.sound_speed)
