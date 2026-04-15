module KWaveMKLExt

# MKL.jl must be loaded BEFORE FFTW so that FFTW picks up Intel MKL as its
# FFT backend instead of the default open-source FFTW library.
#
# On x86 systems this replaces FFTW's pthreads-based parallelism with Intel
# MKL's OpenMP-based FFT, which typically achieves 3–5× higher throughput on
# large 2D/3D transforms — bringing Julia's FFT performance in line with MATLAB.
#
# Usage (add MKL to your environment, then load before KWave):
#
#   using MKL        # must be first
#   using KWave
#
# Or equivalently: `julia -e 'import Pkg; Pkg.add("MKL")' && julia -e 'using MKL, KWave'`
#
# Note: MKL.jl is Intel-only and has no effect on Apple Silicon (aarch64) or
# other non-x86 architectures.  On those platforms, multi-threaded FFTW
# (`julia -t auto`) is the best available option.

import MKL
import KWave
import FFTW

function __init__()
    if Sys.ARCH === :x86_64 || Sys.ARCH === :i686
        # MKL was loaded before FFTW; nothing extra to do — FFTW has already
        # switched its backend to MKL at import time via MKL.jl's __init__.
        @info "KWave: MKL FFT backend active ($(Sys.ARCH)). " *
              "Large 2D/3D FFT throughput should match MATLAB."
    else
        @warn "KWave: MKL.jl loaded but architecture is $(Sys.ARCH) — " *
              "MKL is x86-only, no speedup applies. " *
              "On Apple Silicon, use `julia -t auto` for multi-threaded FFTW."
    end
end

end # module KWaveMKLExt
