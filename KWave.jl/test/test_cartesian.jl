@testset "Cartesian Geometry (Phase 3)" begin
    @testset "make_cart_circle" begin
        coords = make_cart_circle(0.01, 100, (0.0, 0.0))
        @test size(coords) == (2, 100)
        # All points should be approximately on the circle
        radii = sqrt.(coords[1, :].^2 .+ coords[2, :].^2)
        @test all(isapprox.(radii, 0.01; atol=1e-10))
    end

    @testset "make_cart_circle partial arc" begin
        coords = make_cart_circle(0.01, 50, (0.0, 0.0); arc_angle=Float64(π))
        @test size(coords) == (2, 50)
        # All points should have y >= 0 (within tolerance)
        @test all(coords[2, :] .>= -1e-8)
    end

    @testset "make_cart_sphere" begin
        coords = make_cart_sphere(0.01, 200, (0.0, 0.0, 0.0))
        @test size(coords) == (3, 200)
        radii = sqrt.(coords[1, :].^2 .+ coords[2, :].^2 .+ coords[3, :].^2)
        @test all(isapprox.(radii, 0.01; atol=1e-10))
    end

    @testset "make_cart_arc" begin
        coords = make_cart_arc((0.0, 0.0), 0.01, 0.005, (0.01, 0.0), 50)
        @test size(coords) == (2, 50)
        radii = sqrt.(coords[1, :].^2 .+ coords[2, :].^2)
        @test all(isapprox.(radii, 0.01; atol=1e-10))
    end

    @testset "make_cart_bowl" begin
        coords = make_cart_bowl((0.0, 0.0, 0.0), 0.01, 0.005,
                                (0.01, 0.0, 0.0), 100)
        @test size(coords) == (3, 100)
        radii = sqrt.(coords[1, :].^2 .+ coords[2, :].^2 .+ coords[3, :].^2)
        @test all(isapprox.(radii, 0.01; atol=1e-10))
    end

    @testset "make_cart_disc" begin
        coords = make_cart_disc((0.0, 0.0, 0.0), 0.01, (0.0, 0.0, 1.0), 100)
        @test size(coords) == (3, 100)
        # All points should be on the z=0 plane
        @test all(isapprox.(coords[3, :], 0.0; atol=1e-10))
        # All points within radius
        radii = sqrt.(coords[1, :].^2 .+ coords[2, :].^2)
        @test all(radii .<= 0.01 + 1e-10)
    end

    @testset "make_cart_rect" begin
        coords = make_cart_rect((0.0, 0.0, 0.0), 0.02, 0.01,
                                (0.0, 0.0, 1.0), 100)
        @test size(coords, 1) == 3
        @test size(coords, 2) > 0
        # All points on z=0 plane
        @test all(isapprox.(coords[3, :], 0.0; atol=1e-10))
    end
end
