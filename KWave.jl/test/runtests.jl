using KWave
using Test
using LinearAlgebra
using FFTW
using Statistics
using HDF5

@testset "KWave.jl" begin
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
end
