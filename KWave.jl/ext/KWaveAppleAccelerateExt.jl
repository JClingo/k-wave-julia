module KWaveAppleAccelerateExt

# Apple's Accelerate framework (vDSP) provides FFT routines optimized for
# Apple Silicon (NEON/AMX).  On aarch64 macOS, FFTW_jll is compiled without
# ARM SIMD flags, so vDSP can be significantly faster for 1D sizes.
#
# AppleAccelerate.jl v0.6 does NOT implement a usable plan interface for rfft
# (FFTSetup supports neither `*` nor `LinearAlgebra.mul!`), so we use the
# functional API (rfft / irfft) wrapped in thin structs that satisfy the
# mul!(out, plan, x) contract expected by spectral_gradient!.
#
# Allocation note: AppleAccelerate.rfft / irfft allocate their output array.
# This is unavoidable with the current library API.  On Apple Silicon the
# allocator is fast and vDSP compute savings dominate for the small 1D sizes
# used in k-Wave 1D simulations.
#
# 2D and 3D grids continue to use FFTW — AppleAccelerate does not expose a
# multi-dimensional rfft.
#
# Usage:
#   using AppleAccelerate   # must be loaded before or alongside KWave
#   using KWave

import AppleAccelerate
import FFTW
import KWave
import KWave: FFTPlans
import LinearAlgebra

# ── Thin plan wrappers around the functional vDSP API ───────────────────────

struct AccelForwardPlan{T<:AbstractFloat}
    n_real::Int
end

struct AccelInversePlan{T<:AbstractFloat}
    n_real::Int
end

function LinearAlgebra.mul!(out::Vector{Complex{T}},
                            plan::AccelForwardPlan{T},
                            x::Vector{T}) where T<:AbstractFloat
    copyto!(out, AppleAccelerate.rfft(x))
    return out
end

function LinearAlgebra.mul!(out::Vector{T},
                            plan::AccelInversePlan{T},
                            x::Vector{Complex{T}}) where T<:AbstractFloat
    copyto!(out, AppleAccelerate.irfft(x, plan.n_real))
    return out
end

# ── __init__ ─────────────────────────────────────────────────────────────────

function __init__()
    if !(Sys.isapple() && Sys.ARCH === :aarch64)
        @warn "KWave: KWaveAppleAccelerateExt loaded on non-Apple-Silicon platform — no effect."
        return
    end

    # Verify normalization: forward=none, inverse=1/N (matches FFTW convention).
    ok = true
    for (T, atol) in ((Float32, 1f-4), (Float64, 1e-10))
        x   = T[1, 2, 3, 4]
        fwd = AppleAccelerate.rfft(x)
        inv = AppleAccelerate.irfft(fwd, 4)
        if !isapprox(inv, x; atol=atol)
            ok = false
            break
        end
    end

    if !ok
        @warn "KWave: AppleAccelerate irfft normalization mismatch — " *
              "falling back to FFTW for all grids."
        return
    end

    # AppleAccelerate.jl v0.6 does not expose a plan-based rfft interface with
    # mul! support.  The functional rfft(x) creates a new vDSP setup per call,
    # making it ~200x slower than planned FFTW even without ARM SIMD.
    # FFTW with a single thread (the default for small 1D grids — see
    # fft_utils.jl) already outperforms MATLAB for 1D sizes.
    #
    # _CREATE_1D_PLANS is intentionally left as nothing so that the base
    # create_fft_plans fallback (single-threaded FFTW) is used.
    # This hook is reserved for a future AppleAccelerate release that provides
    # an allocation-free plan interface for rfft.
    @info "KWave: KWaveAppleAccelerateExt loaded. " *
          "AppleAccelerate.jl v0.6 lacks a plan-based rfft interface; " *
          "using single-threaded FFTW for 1D (already faster than MATLAB on Apple Silicon)."
end

end # module KWaveAppleAccelerateExt
