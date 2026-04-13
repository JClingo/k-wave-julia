@testset "Geometry Extended (Phase 3)" begin
    @testset "make_arc" begin
        Nx, Ny = 64, 64
        cx, cy = 32, 32
        radius = 20.0
        diameter = 30.0
        focus = (32, 10)
        arc = make_arc(Nx, Ny, cx, cy, radius, diameter, focus)
        @test size(arc) == (Nx, Ny)
        @test arc isa BitMatrix
        @test count(arc) > 0
    end

    @testset "make_line 2D" begin
        Nx, Ny = 32, 32
        line = make_line(Nx, Ny, (5, 5), (25, 25))
        @test size(line) == (Nx, Ny)
        @test line[5, 5] == true
        @test line[25, 25] == true
        @test count(line) >= 20  # diagonal line has at least 21 points
    end

    @testset "make_line 3D" begin
        Nx, Ny, Nz = 16, 16, 16
        line = make_line(Nx, Ny, Nz, (2, 2, 2), (14, 14, 14))
        @test size(line) == (Nx, Ny, Nz)
        @test line[2, 2, 2] == true
        @test line[14, 14, 14] == true
        @test count(line) >= 12
    end

    @testset "make_bowl" begin
        Nx, Ny, Nz = 32, 32, 32
        center = (16, 16, 16)
        radius = 12.0
        diameter = 16.0
        focus = (16, 16, 5)
        bowl = make_bowl(Nx, Ny, Nz, center, radius, diameter, focus)
        @test size(bowl) == (Nx, Ny, Nz)
        @test count(bowl) > 0
    end

    @testset "make_multi_arc" begin
        Nx, Ny = 64, 64
        arcs = [
            (32, 32, 15.0, 20.0, (32, 10)),
            (32, 32, 15.0, 20.0, (32, 54)),
        ]
        mask = make_multi_arc(Nx, Ny, arcs)
        @test count(mask) > 0
    end

    @testset "make_multi_bowl" begin
        Nx, Ny, Nz = 32, 32, 32
        bowls = [
            ((16, 16, 16), 10.0, 12.0, (16, 16, 5)),
            ((16, 16, 16), 10.0, 12.0, (16, 16, 27)),
        ]
        mask = make_multi_bowl(Nx, Ny, Nz, bowls)
        @test count(mask) > 0
    end

    @testset "make_spherical_section" begin
        Nx, Ny, Nz = 32, 32, 32
        section = make_spherical_section(Nx, Ny, Nz, (16, 16, 16), 10.0, 5.0)
        @test size(section) == (Nx, Ny, Nz)
        @test count(section) > 0
    end
end
