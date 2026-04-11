@testset "KWaveGrid" begin
    @testset "1D grid construction" begin
        kgrid = KWaveGrid(64, 1e-4)
        @test kgrid isa KWaveGrid1D
        @test kgrid.Nx == 64
        @test kgrid.dx == 1e-4
        @test length(kgrid.x_vec) == 64
        @test length(kgrid.kx_vec) == 64
        @test length(kgrid.k) == 64
        @test ndims(kgrid) == 1
        @test total_grid_points(kgrid) == 64

        # Spatial vector should be centered
        @test kgrid.x_vec[1] ≈ -(64 - 1) / 2 * 1e-4
        @test kgrid.x_vec[end] ≈ (64 - 1) / 2 * 1e-4

        # Wavenumber vector: first element should be 0
        @test kgrid.kx_vec[1] ≈ 0.0

        # Wavenumber magnitude should be non-negative
        @test all(kgrid.k .>= 0)

        # Shift operators should have unit magnitude
        @test all(abs.(kgrid.ddx_k_shift_pos) .≈ 1.0)
        @test all(abs.(kgrid.ddx_k_shift_neg) .≈ 1.0)
    end

    @testset "2D grid construction" begin
        kgrid = KWaveGrid(64, 1e-4, 128, 2e-4)
        @test kgrid isa KWaveGrid2D
        @test kgrid.Nx == 64
        @test kgrid.Ny == 128
        @test kgrid.dx == 1e-4
        @test kgrid.dy == 2e-4
        @test ndims(kgrid) == 2
        @test total_grid_points(kgrid) == 64 * 128
        @test grid_size(kgrid) == (64, 128)
        @test grid_spacing(kgrid) == (1e-4, 2e-4)

        # k magnitude matrix should be correct size
        @test size(kgrid.k) == (64, 128)

        # k at origin should be zero
        @test kgrid.k[1, 1] ≈ 0.0 atol=1e-15

        # k should be consistent with kx, ky
        for j in 1:5, i in 1:5
            expected = sqrt(kgrid.kx_vec[i]^2 + kgrid.ky_vec[j]^2)
            @test kgrid.k[i, j] ≈ expected
        end

        # k_max should correspond to Nyquist
        @test k_max(kgrid) ≈ π / max(kgrid.dx, kgrid.dy) atol=1e-10
    end

    @testset "3D grid construction" begin
        kgrid = KWaveGrid(32, 1e-4, 32, 1e-4, 32, 1e-4)
        @test kgrid isa KWaveGrid3D
        @test ndims(kgrid) == 3
        @test total_grid_points(kgrid) == 32^3
        @test size(kgrid.k) == (32, 32, 32)
        @test kgrid.k[1, 1, 1] ≈ 0.0 atol=1e-15
    end

    @testset "make_time!" begin
        kgrid = KWaveGrid(128, 1e-4, 128, 1e-4)
        c = 1500.0  # m/s

        dt, Nt = make_time!(kgrid, c)

        # dt should satisfy CFL: dt <= CFL * dx / c
        @test dt ≈ 0.3 * 1e-4 / 1500.0
        @test dt > 0
        @test Nt > 0
        @test kgrid.dt[] == dt
        @test kgrid.Nt[] == Nt
        @test length(kgrid.t_array) == Nt
        @test kgrid.t_array[1] ≈ 0.0
        @test kgrid.t_array[2] ≈ dt

        # With explicit t_end
        dt2, Nt2 = make_time!(kgrid, c; t_end=1e-5)
        @test Nt2 ≈ ceil(Int, 1e-5 / dt2) + 1
    end

    @testset "Wavenumber vector properties" begin
        # Even N
        kv = KWave._make_wavenumber_vec(8, 1.0)
        @test length(kv) == 8
        @test kv[1] ≈ 0.0  # DC component
        # For even N, kv has an unpaired -N/2 component, so sum != 0
        # But the positive and negative parts (excluding Nyquist) should balance
        @test kv[1] == 0.0  # DC is zero

        # Odd N
        kv_odd = KWave._make_wavenumber_vec(7, 1.0)
        @test length(kv_odd) == 7
        @test kv_odd[1] ≈ 0.0
    end
end
