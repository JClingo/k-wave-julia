using KWave
using Test
using LinearAlgebra
using FFTW
using Statistics
using HDF5

@testset "KWave.jl" begin
    # Phase 1 & 2 tests
    include("test_grid.jl")
    include("test_pml.jl")
    include("test_fft.jl")
    include("test_utils.jl")
    include("test_geometry.jl")
    include("test_signal.jl")
    include("test_filter.jl")
    include("test_io.jl")
    include("test_solver_2d.jl")
    include("test_solver_1d.jl")
    include("test_solver_3d.jl")
    include("test_absorption.jl")
    include("test_time_reversal.jl")

    # Phase 3 tests
    include("test_geometry_extended.jl")
    include("test_cartesian.jl")
    include("test_signal_processing.jl")
    include("test_materials.jl")
    include("test_array.jl")
    include("test_transducer.jl")
    include("test_reconstruction.jl")
    include("test_cw.jl")
    include("test_axisymmetric.jl")
    include("test_visualization.jl")

    # Phase 4 tests
    include("test_elastic.jl")
    include("test_thermal.jl")
    include("test_analytical.jl")
    include("test_beamform.jl")
    include("test_voxel.jl")
    include("test_python_interop.jl")
end
