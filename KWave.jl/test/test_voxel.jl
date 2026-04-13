@testset "Voxel and 3D Visualization" begin
    @testset "voxel_plot fallback" begin
        vol = randn(10, 10, 10)
        # Without a Makie backend, should return nothing with a warning
        result = @test_logs (:warn,) voxel_plot(vol)
        @test result === nothing
    end

    @testset "isosurface_plot fallback" begin
        vol = randn(10, 10, 10)
        result = @test_logs (:warn,) isosurface_plot(vol, 0.5)
        @test result === nothing
    end

    @testset "max_intensity_projection single dim" begin
        vol = zeros(10, 12, 14)
        vol[5, 6, 7] = 1.0

        mip = max_intensity_projection(vol; dims=3)
        @test size(mip) == (10, 12)
        @test maximum(mip) == 1.0

        mip2 = max_intensity_projection(vol; dims=2)
        @test size(mip2) == (10, 14)

        mip1 = max_intensity_projection(vol; dims=1)
        @test size(mip1) == (12, 14)
    end

    @testset "max_intensity_projection all dims" begin
        vol = zeros(8, 10, 12)
        vol[4, 5, 6] = 2.0

        result = max_intensity_projection(vol; dims=:all)
        @test haskey(result, :xy)
        @test haskey(result, :xz)
        @test haskey(result, :yz)
        @test size(result.xy) == (8, 10)
        @test size(result.xz) == (8, 12)
        @test size(result.yz) == (10, 12)
    end

    @testset "max_intensity_projection dB scale" begin
        vol = zeros(8, 8, 8)
        vol[4, 4, 4] = 1.0
        vol[2, 2, 2] = 0.01

        mip = max_intensity_projection(vol; dims=3, db_scale=true, db_range=40)
        @test all(mip .>= -40)
        @test maximum(mip) ≈ 0.0 atol=1e-10
    end

    @testset "_apply_db_scale" begin
        field = [1.0 0.1; 0.01 0.0]
        db = KWave._apply_db_scale(field, 60)
        @test db[1, 1] ≈ 0.0 atol=1e-10
        @test db[1, 2] ≈ -20.0 atol=0.1
        @test db[2, 1] ≈ -40.0 atol=0.1
        @test db[2, 2] >= -60.0
    end
end
