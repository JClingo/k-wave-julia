# ============================================================================
# KWave.jl Unitful Extension — Physical units integration
# ============================================================================

module KWaveUnitfulExt

using KWave
using Unitful
using Unitful: m, s, Hz, Pa, kg, K, W, J

"""
    strip_units(x)

Remove Unitful units from a value, converting to SI base units first.
"""
strip_units(x::Unitful.Quantity) = ustrip(uconvert(Unitful.NoUnits, x))
strip_units(x::Real) = Float64(x)
strip_units(x::Nothing) = nothing

"""
    strip_to_SI(x, unit)

Convert a Unitful quantity to the specified SI unit and strip.
"""
strip_to_SI(x::Unitful.Quantity, unit) = Float64(ustrip(uconvert(unit, x)))
strip_to_SI(x::Real, unit) = Float64(x)
strip_to_SI(x::Nothing, unit) = nothing

function strip_to_SI(x::AbstractArray{<:Unitful.Quantity}, unit)
    return Float64[ustrip(uconvert(unit, xi)) for xi in x]
end
strip_to_SI(x::AbstractArray{<:Real}, unit) = Float64.(x)

"""
    KWaveGrid(nx, dx::Unitful.Length)
    KWaveGrid(nx, dx::Unitful.Length, ny, dy::Unitful.Length)
    KWaveGrid(nx, dx::Unitful.Length, ny, dy::Unitful.Length, nz, dz::Unitful.Length)

Construct a KWaveGrid from dimensions with physical units.
Grid spacings are automatically converted to meters.
"""
function KWave.KWaveGrid(nx::Int, dx::Unitful.Length)
    return KWave.KWaveGrid(nx, strip_to_SI(dx, m))
end

function KWave.KWaveGrid(nx::Int, dx::Unitful.Length, ny::Int, dy::Unitful.Length)
    return KWave.KWaveGrid(nx, strip_to_SI(dx, m), ny, strip_to_SI(dy, m))
end

function KWave.KWaveGrid(nx::Int, dx::Unitful.Length, ny::Int, dy::Unitful.Length,
                         nz::Int, dz::Unitful.Length)
    return KWave.KWaveGrid(nx, strip_to_SI(dx, m), ny, strip_to_SI(dy, m),
                           nz, strip_to_SI(dz, m))
end

"""
    KWaveMedium(; sound_speed::Unitful.Velocity, density::Unitful.Density, ...)

Construct a KWaveMedium with physical units. Values are automatically
converted to SI base units (m/s, kg/m³, etc.).
"""
function KWave.KWaveMedium(;
    sound_speed::Union{Unitful.Velocity, AbstractArray{<:Unitful.Velocity}},
    density::Union{Unitful.Density, AbstractArray{<:Unitful.Density}}=1000.0u"kg/m^3",
    kwargs...)

    c = sound_speed isa AbstractArray ?
        Float64[ustrip(uconvert(m/s, ci)) for ci in sound_speed] :
        Float64(ustrip(uconvert(m/s, sound_speed)))

    rho = density isa AbstractArray ?
        Float64[ustrip(uconvert(kg/m^3, di)) for di in density] :
        Float64(ustrip(uconvert(kg/m^3, density)))

    return KWave.KWaveMedium(; sound_speed=c, density=rho, kwargs...)
end

"""
    make_time!(kgrid, sound_speed::Unitful.Velocity; kwargs...)

Set the time array with a sound speed that has physical units.
"""
function KWave.make_time!(kgrid::KWave.AbstractKWaveGrid,
                          sound_speed::Union{Unitful.Velocity, AbstractArray{<:Unitful.Velocity}};
                          kwargs...)
    c = sound_speed isa AbstractArray ?
        Float64[ustrip(uconvert(m/s, ci)) for ci in sound_speed] :
        Float64(ustrip(uconvert(m/s, sound_speed)))
    return KWave.make_time!(kgrid, c; kwargs...)
end

"""
    water_sound_speed(T::Unitful.Temperature)

Compute water sound speed from a temperature with units.
"""
function KWave.water_sound_speed(T::Unitful.Temperature)
    T_celsius = ustrip(uconvert(u"°C", T))
    return KWave.water_sound_speed(T_celsius)
end

"""
    water_density(T::Unitful.Temperature)

Compute water density from a temperature with units.
"""
function KWave.water_density(T::Unitful.Temperature)
    T_celsius = ustrip(uconvert(u"°C", T))
    return KWave.water_density(T_celsius)
end

"""
    water_absorption(f::Unitful.Frequency, T::Unitful.Temperature)

Compute water absorption from frequency and temperature with units.
"""
function KWave.water_absorption(f::Unitful.Frequency, T::Unitful.Temperature)
    f_hz = ustrip(uconvert(Hz, f))
    T_celsius = ustrip(uconvert(u"°C", T))
    return KWave.water_absorption(f_hz, T_celsius)
end

end # module KWaveUnitfulExt
