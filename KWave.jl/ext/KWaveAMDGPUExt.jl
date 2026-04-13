# ============================================================================
# KWave.jl AMDGPU Extension — AMD GPU compute backend
# ============================================================================

module KWaveAMDGPUExt

using KWave
using AMDGPU

"""
    gpu_amd(x)

Move an array to AMD GPU (ROCArray). Scalars and `nothing` pass through.
"""
gpu_amd(x::AbstractArray) = ROCArray(x)
gpu_amd(x::Real) = x
gpu_amd(x::Nothing) = nothing

"""
    cpu(x::ROCArray)

Move a ROCArray back to CPU.
"""
cpu(x::ROCArray) = Array(x)
cpu(x::AbstractArray) = x
cpu(x::Real) = x
cpu(x::Nothing) = nothing

"""
    AMDGPUFFTPlans

GPU FFT plans using rocFFT via AMDGPU.jl.
"""
struct AMDGPUFFTPlans
    forward::Any
    inverse::Any
end

"""
    create_amdgpu_fft_plans(dims)

Create in-place FFT/IFFT plans for ROCArrays.
"""
function create_amdgpu_fft_plans(dims::Tuple)
    tmp = AMDGPU.zeros(ComplexF64, dims...)
    fwd = plan_fft!(tmp)
    inv = plan_ifft!(tmp)
    return AMDGPUFFTPlans(fwd, inv)
end

"""
    to_gpu_medium(medium::KWaveMedium; T=Float64)

Convert a KWaveMedium's array fields to ROCArrays for AMD GPU computation.
"""
function to_amdgpu_medium(medium::KWave.KWaveMedium{T}) where T
    _g = x -> x isa AbstractArray ? ROCArray{T}(x) : (x === nothing ? nothing : T(x))
    return KWave.KWaveMedium{T}(
        _g(medium.sound_speed),
        _g(medium.density),
        _g(medium.alpha_coeff),
        medium.alpha_power === nothing ? nothing : T(medium.alpha_power),
        medium.alpha_mode,
        _g(medium.BonA),
    )
end

"""
    to_amdgpu_source(source)

Convert a KWaveSource's array fields to ROCArrays.
"""
function to_amdgpu_source(source::KWave.KWaveSource{T}) where T
    _g = x -> x === nothing ? nothing : ROCArray{T}(x)
    _gb = x -> x === nothing ? nothing : ROCArray(x)
    return KWave.KWaveSource{T}(
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

end # module KWaveAMDGPUExt
