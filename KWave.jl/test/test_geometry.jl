@testset "Geometry" begin
    @testset "make_disc" begin
        disc = make_disc(64, 64, 32, 32, 10)
        @test disc isa BitMatrix
        @test size(disc) == (64, 64)
        # Center should be filled
        @test disc[32, 32] == true
        # Far corner should be empty
        @test disc[1, 1] == false
        # Should be roughly circular: count should be ≈ π*r²
        area = count(disc)
        @test abs(area - π * 10^2) / (π * 10^2) < 0.1  # within 10% of analytic

        # Symmetry
        @test disc[32-5, 32] == disc[32+5, 32]
        @test disc[32, 32-5] == disc[32, 32+5]
    end

    @testset "make_circle" begin
        circle = make_circle(64, 64, 32, 32, 15)
        @test circle isa BitMatrix
        @test size(circle) == (64, 64)
        # Center should be empty (it's a perimeter)
        @test circle[32, 32] == false
        # Points on the circle should be set
        @test circle[32+15, 32] == true || circle[32-15, 32] == true
        # Total points should be roughly 2πr
        n_points = count(circle)
        @test n_points > 0
        @test abs(n_points - 2π * 15) / (2π * 15) < 0.3
    end

    @testset "make_ball" begin
        ball = make_ball(32, 32, 32, 16, 16, 16, 8)
        @test ball isa BitArray{3}
        @test size(ball) == (32, 32, 32)
        @test ball[16, 16, 16] == true
        @test ball[1, 1, 1] == false
        # Volume should be ≈ (4/3)π*r³
        vol = count(ball)
        expected_vol = (4/3) * π * 8^3
        @test abs(vol - expected_vol) / expected_vol < 0.15
    end

    @testset "make_sphere" begin
        sphere = make_sphere(32, 32, 32, 16, 16, 16, 10)
        @test sphere isa BitArray{3}
        @test size(sphere) == (32, 32, 32)
        # Center should be empty (it's a shell)
        @test sphere[16, 16, 16] == false
        # Surface points should exist
        @test count(sphere) > 0
        # Surface area should be ≈ 4π*r²
        n_points = count(sphere)
        expected_area = 4π * 10^2
        @test abs(n_points - expected_area) / expected_area < 0.3
    end
end
