# ============================================================================
# KWave.jl Metal Extension — Apple GPU compute backend
# ============================================================================

module KWaveMetalExt

using KWave
using Metal

"""
    gpu_metal(x)

Move an array to Metal GPU (MtlArray). Scalars and `nothing` pass through.
"""
gpu_metal(x::AbstractArray) = MtlArray(x)
gpu_metal(x::Real) = x
gpu_metal(x::Nothing) = nothing

"""
    cpu(x::MtlArray)

Move a MtlArray back to CPU.
"""
cpu(x::MtlArray) = Array(x)
cpu(x::AbstractArray) = x
cpu(x::Real) = x
cpu(x::Nothing) = nothing

"""
    MetalFFTPlans

GPU FFT plans using Metal.
Note: Metal.jl FFT support may be limited compared to CUDA/AMDGPU.
"""
struct MetalFFTPlans
    forward::Any
    inverse::Any
end

"""
    create_metal_fft_plans(dims)

Create FFT/IFFT plans for MtlArrays.
"""
function create_metal_fft_plans(dims::Tuple)
    # Metal.jl FFT support — create plans for MtlArrays
    tmp = Metal.zeros(ComplexF32, dims...)
    fwd = plan_fft!(tmp)
    inv = plan_ifft!(tmp)
    return MetalFFTPlans(fwd, inv)
end

"""
    to_metal_medium(medium)

Convert a KWaveMedium's array fields to MtlArrays for Metal GPU computation.
Metal typically uses Float32 for best performance.
"""
function to_metal_medium(medium::KWave.KWaveMedium{T}) where T
    _g = x -> x isa AbstractArray ? MtlArray{Float32}(x) : (x === nothing ? nothing : Float32(x))
    return KWave.KWaveMedium{Float32}(
        _g(medium.sound_speed),
        _g(medium.density),
        _g(medium.alpha_coeff),
        medium.alpha_power === nothing ? nothing : Float32(medium.alpha_power),
        medium.alpha_mode,
        _g(medium.BonA),
    )
end

"""
    to_metal_source(source)

Convert a KWaveSource's array fields to MtlArrays.
"""
function to_metal_source(source::KWave.KWaveSource{T}) where T
    _g = x -> x === nothing ? nothing : MtlArray{Float32}(x)
    _gb = x -> x === nothing ? nothing : MtlArray(x)
    return KWave.KWaveSource{Float32}(
        _g(source.p0),
        _gb(source.p_mask),
        _g(source.p),
        source.p_mode,
        _gb(source.u_mask),
        _g(source.ux),
        _g(source.uy),
        _g(source.uz),
        source.u_mode,
    )
end

end # module KWaveMetalExt
