module KWave

using FFTW
using HDF5
using LinearAlgebra
using Statistics
using DSP
using Interpolations
using ProgressMeter

# Types and enums
include("types.jl")
export SourceMode, Dirichlet, Additive, AdditiveNoCorrection

# Data structures
include("grid.jl")
include("medium.jl")
include("source.jl")
include("sensor.jl")
export AbstractKWaveGrid, KWaveGrid1D, KWaveGrid2D, KWaveGrid3D
export KWaveGrid, make_time!
export total_grid_points, grid_size, grid_spacing, k_max
export KWaveMedium, KWaveSource, KWaveSensor, SimulationOutput
export is_homogeneous, is_lossless, is_nonlinear
export is_time_reversal

# Core utilities
include("pml.jl")
include("fft_utils.jl")
include("utils.jl")
export get_pml, get_optimal_pml_size
export cart2grid, grid2cart, expand_matrix, resize_array

# Geometry
include("geometry/shapes.jl")
export make_disc, make_circle, make_ball, make_sphere

# Signal & filtering
include("signal/generation.jl")
include("filter/filters.jl")
include("material/conversion.jl")
export tone_burst, gaussian_pulse
export smooth, apply_filter, gaussian_filter, get_win
export db2neper, neper2db

# I/O
include("io/hdf5.jl")
export write_matrix, write_grid, write_flags, write_attributes
export read_matrix, read_output

# Visualization
include("visualization/colormap.jl")
include("visualization/field_display.jl")
export get_color_map
export SimulationDisplay, NullDisplay
export create_sim_display, update_sim_display!, close_sim_display!
export MovieRecorder, NullRecorder
export create_movie_recorder, record_frame!, finalize_movie!

# Solver
include("solver/first_order_steps.jl")
include("solver/first_order.jl")
export kspace_first_order
export AbsorptionParams

end # module KWave
