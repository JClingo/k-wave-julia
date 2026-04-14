# Examples

All example scripts are in the [`examples/`](https://github.com/JClingo/k-wave-julia/tree/main/KWave.jl/examples) directory.

## Initial Value Problems (IVP)

Simulate photoacoustic wave propagation from an initial pressure distribution.

| Example | Description |
|---|---|
| [`ivp_homogeneous_medium_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/ivp_homogeneous_medium_2d.jl) | 2D IVP in a lossless homogeneous medium |
| [`ivp_heterogeneous_medium_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/ivp_heterogeneous_medium_2d.jl) | 2D IVP with spatially varying sound speed |
| [`ivp_absorption_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/ivp_absorption_2d.jl) | 2D IVP with power-law absorption |
| [`ivp_1d_simulation.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/ivp_1d_simulation.jl) | 1D propagation — comparison with analytical solution |
| [`ivp_3d_simulation.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/ivp_3d_simulation.jl) | 3D IVP with full pressure field output |

## Time-Varying Source Problems (TVSP)

Simulate propagation from a driven transducer source.

| Example | Description |
|---|---|
| [`tvsp_homogeneous_medium_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/tvsp_homogeneous_medium_2d.jl) | 2D CW source in homogeneous medium |
| [`tvsp_nonlinear_propagation_1d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/tvsp_nonlinear_propagation_1d.jl) | 1D nonlinear wave distortion (comparison with Mendousse solution) |

## Photoacoustic Reconstruction (PR)

Reconstruct initial pressure distributions from time-series sensor data.

| Example | Description |
|---|---|
| [`pr_fft_reconstruction_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/pr_fft_reconstruction_2d.jl) | 2D k-space FFT reconstruction from a line sensor |
| [`pr_time_reversal_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/pr_time_reversal_2d.jl) | 2D time-reversal reconstruction |

## Sensor Directivity (SD)

| Example | Description |
|---|---|
| [`sd_sensor_directivity_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/sd_sensor_directivity_2d.jl) | 2D directional sensor frequency response |

## Ultrasound (US)

Pulse-echo and beamforming simulations.

| Example | Description |
|---|---|
| [`us_beamforming_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/us_beamforming_2d.jl) | 2D delay-and-sum beamforming with `KWaveTransducer` |
| [`us_phased_array_3d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/us_phased_array_3d.jl) | 3D phased array with `KWaveArray` |

## Elastic Wave Propagation (EWP)

| Example | Description |
|---|---|
| [`ewp_elastic_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/ewp_elastic_2d.jl) | 2D P and S wave propagation |

## Diffusion / Bioheat (DIFF)

| Example | Description |
|---|---|
| [`diff_bioheat_1d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/diff_bioheat_1d.jl) | 1D bioheat with comparison to analytical solution |

## Visualization

| Example | Description |
|---|---|
| [`vis_cairomakie_2d.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/vis_cairomakie_2d.jl) | Static 2D plots with CairoMakie |
| [`vis_glmakie_interactive.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/vis_glmakie_interactive.jl) | Interactive 3D visualization with GLMakie |
| [`vis_wglmakie_notebook.jl`](https://github.com/JClingo/k-wave-julia/blob/main/KWave.jl/examples/vis_wglmakie_notebook.jl) | WGLMakie Jupyter notebook visualization |
