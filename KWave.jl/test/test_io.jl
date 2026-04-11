@testset "HDF5 I/O" begin
    test_file = tempname() * ".h5"

    @testset "write_matrix / read_matrix round-trip" begin
        data = randn(16, 16)
        write_matrix(test_file, "test_data", data)

        loaded = read_matrix(test_file, "test_data")
        @test loaded ≈ data

        rm(test_file; force=true)
    end

    @testset "write_matrix scalar" begin
        write_matrix(test_file, "scalar_val", 42.0)
        loaded = read_matrix(test_file, "scalar_val")
        @test loaded[1] ≈ 42.0

        rm(test_file; force=true)
    end

    @testset "write_grid" begin
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        make_time!(kgrid, 1500.0)

        write_grid(test_file, kgrid; pml_size=10, pml_alpha=2.0)

        # Verify contents
        nx = read_matrix(test_file, "Nx")
        @test nx[1] == 64
        ny = read_matrix(test_file, "Ny")
        @test ny[1] == 64
        dt_val = read_matrix(test_file, "dt")
        @test dt_val[1] ≈ kgrid.dt[]

        rm(test_file; force=true)
    end

    @testset "write_flags" begin
        flags = Dict("ux_source_flag" => 0, "p_source_flag" => 1)
        write_flags(test_file, flags)

        loaded = read_matrix(test_file, "ux_source_flag")
        @test loaded[1] == 0

        rm(test_file; force=true)
    end

    @testset "write_attributes" begin
        write_attributes(test_file)

        HDF5.h5open(test_file, "r") do fid
            @test HDF5.attrs(fid)["created_by"] == "KWave.jl"
            @test haskey(HDF5.attrs(fid), "creation_date")
        end

        rm(test_file; force=true)
    end
end
