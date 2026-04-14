# GPU Acceleration

KWave.jl supports GPU acceleration via package extensions.

## Available Backends

| Backend | Extension | Hardware |
|---|---|---|
| `CUDA` | `KWaveCUDAExt` | NVIDIA GPUs |
| `Metal` | `KWaveMetalExt` | Apple Silicon (M-series) |
| `AMDGPU` | `KWaveAMDGPUExt` | AMD GPUs (ROCm) |

## Usage

Load the GPU package alongside KWave:

```julia
using KWave, CUDA
# or: using KWave, Metal
# or: using KWave, AMDGPU
```

Then pass `data_cast=Float32` to the solver:

```julia
output = kspace_first_order(kgrid, medium, source, sensor; data_cast=Float32)
```

The extension automatically moves arrays to the GPU before the time loop and returns results as CPU arrays.

## Requirements

- `data_cast=Float32` is required — most consumer GPUs do not support `Float64` at useful speeds.
- GPU must have sufficient VRAM for the simulation arrays:
  - A 512×512×512 grid at `Float32` requires ~4 GB for pressure + velocity fields alone.
- FFTW is replaced by cuFFT (CUDA), Metal FFT, or rocFFT automatically.

## Performance Notes

- First call is slow due to kernel compilation. Subsequent calls with the same grid size are fast.
- 2D simulations on GPU show less speedup than 3D — the problem size must be large enough to saturate the GPU.
- Heterogeneous media (array-valued `sound_speed`, `density`) benefit more from GPU than homogeneous.

## See Also

[`kspace_first_order`](@ref), [Running Simulations](simulation.md)
