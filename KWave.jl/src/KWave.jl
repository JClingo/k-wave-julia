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
include("geometry/cartesian.jl")
export make_disc, make_circle, make_ball, make_sphere
export make_arc, make_line, make_bowl
export make_multi_arc, make_multi_bowl, make_spherical_section
export make_cart_circle, make_cart_sphere, make_cart_arc
export make_cart_bowl, make_cart_disc, make_cart_rect

# Signal & filtering
include("signal/generation.jl")
include("signal/processing.jl")
include("filter/filters.jl")
include("material/conversion.jl")
include("material/properties.jl")
export tone_burst, gaussian_pulse
export smooth, apply_filter, gaussian_filter, get_win
export db2neper, neper2db
# Phase 3: extended signal processing
export add_noise, create_cw_signals, extract_amp_phase
export log_compression, envelope_detection
export gradient_fd, gradient_spect, filter_time_series, spect
# Phase 3: material properties
export water_sound_speed, water_density, water_absorption, water_nonlinearity
export fit_power_law_params

# I/O
include("io/hdf5.jl")
export write_matrix, write_grid, write_flags, write_attributes
export read_matrix, read_output

# Transducer arrays
include("array.jl")
include("transducer.jl")
export KWaveArray, ArrayElement
export ArcElement, BowlElement, DiscElement, RectElement, SphereElement
export add_arc_element!, add_bowl_element!, add_disc_element!, add_rect_element!
export get_element_binary_mask, get_array_binary_mask
export get_distributed_source_signal, combine_sensor_data
export KWaveTransducer
export get_transducer_binary_mask, get_transducer_source
export combine_transducer_sensor_data

# Reconstruction
include("reconstruction/fft_recon.jl")
include("reconstruction/beamform.jl")
export kspace_line_recon, kspace_plane_recon
# Phase 4: beamforming and scan conversion
export beamform_delay_and_sum, scan_conversion

# Visualization
include("visualization/colormap.jl")
include("visualization/field_display.jl")
include("visualization/plots.jl")
include("visualization/voxel.jl")
export get_color_map
export SimulationDisplay, NullDisplay
export create_sim_display, update_sim_display!, close_sim_display!
export MovieRecorder, NullRecorder
export create_movie_recorder, record_frame!, finalize_movie!
export beam_plot, fly_through, overlay_plot, stacked_plot
# Phase 4: 3D visualization
export voxel_plot, isosurface_plot, max_intensity_projection

# Solver — acoustic (fluid)
include("solver/first_order_steps.jl")
include("solver/first_order.jl")
include("solver/axisymmetric.jl")
include("solver/cw.jl")
export kspace_first_order
export kspace_first_order_as
export acoustic_field_propagator, angular_spectrum_cw
export AbsorptionParams

# Phase 4: Elastic wave solvers
include("solver/elastic.jl")
export ElasticMedium, ElasticSource
export pstd_elastic_2d, pstd_elastic_3d

# Phase 4: Thermal simulation
include("solver/diffusion.jl")
export ThermalMedium, ThermalSource
export kwave_diffusion, bioheat_exact

# Phase 4: Analytical reference solutions
include("reference/analytical.jl")
export focused_annulus_oneil, focused_bowl_oneil, mendousse

# Phase 4: Python interop utilities
include("interop/python.jl")
export python_kwave_grid, python_medium, python_source, python_sensor
export python_run_simulation, python_output_to_dict

# CLI
include("cli.jl")

end # module KWave
