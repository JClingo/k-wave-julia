# ============================================================================
# KWave.jl CUDA Extension — GPU compute backend
# ============================================================================

module KWaveCUDAExt

using KWave
using CUDA
using CUDA.CUFFT

"""
    gpu(x)

Move an array to the GPU (CuArray). Scalars and `nothing` pass through.
"""
gpu(x::AbstractArray) = CuArray(x)
gpu(x::Real) = x
gpu(x::Nothing) = nothing

"""
    cpu(x)

Move a CuArray back to the CPU. Scalars and `nothing` pass through.
"""
cpu(x::CuArray) = Array(x)
cpu(x::AbstractArray) = x
cpu(x::Real) = x
cpu(x::Nothing) = nothing

"""
    CUDAFFTPlans

GPU FFT plans using CUFFT.
"""
struct CUDAFFTPlans
    forward::Any
    inverse::Any
end

"""
    create_cuda_fft_plans(dims)

Create in-place FFT/IFFT plans for CuArrays.
"""
function create_cuda_fft_plans(dims::Tuple)
    tmp = CUDA.zeros(ComplexF32, dims...)
    fwd = plan_fft!(tmp)
    inv = plan_ifft!(tmp)
    return CUDAFFTPlans(fwd, inv)
end

"""
    to_gpu_medium(medium)

Convert a KWaveMedium's array fields to CuArrays for GPU computation.
"""
function to_gpu_medium(medium::KWave.KWaveMedium{T}) where T
    _g = x -> x isa AbstractArray ? CuArray{Float32}(x) : (x === nothing ? nothing : Float32(x))
    return KWave.KWaveMedium{Float32}(
        _g(medium.sound_speed),
        _g(medium.density),
        _g(medium.alpha_coeff),
        medium.alpha_power === nothing ? nothing : Float32(medium.alpha_power),
        medium.alpha_mode,
        _g(medium.BonA)
    )
end

"""
    to_gpu_source(source)

Convert a KWaveSource's array fields to CuArrays.
"""
function to_gpu_source(source::KWave.KWaveSource{T}) where T
    _g = x -> x === nothing ? nothing : CuArray{Float32}(x)
    _gb = x -> x === nothing ? nothing : CuArray(x)
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

end # module KWaveCUDAExt
