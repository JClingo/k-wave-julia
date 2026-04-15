@testset "FFT Utilities" begin
    @testset "FFT plans creation" begin
        plans = KWave.create_fft_plans((32, 32))
        @test plans isa KWave.FFTPlans

        # Plans for a grid
        kgrid = KWaveGrid(32, 1e-4, 32, 1e-4)
        plans2 = KWave.create_fft_plans(kgrid)
        @test plans2 isa KWave.FFTPlans
    end

    @testset "FFT round-trip" begin
        N = 32
        plans = KWave.create_fft_plans((N, N))
        data = randn(Float64, N, N)
        original = copy(data)

        scratch = plans.forward * data
        result = plans.inverse * scratch

        @test result ≈ original atol=1e-12
    end

    @testset "Spectral gradient - 2D sinusoidal" begin
        # Test: gradient of sin(kx * x) should be kx * cos(kx * x)
        Nx, Ny = 64, 64
        dx, dy = 1.0, 1.0
        kgrid = KWaveGrid(Nx, dx, Ny, dy)

        plans = KWave.create_fft_plans(kgrid)
        scratch = zeros(ComplexF64, Nx÷2+1, Ny)

        # Create a sinusoidal field: f(x,y) = sin(2π * x / L)
        Lx = Nx * dx
        kx_test = 2π / Lx
        f = zeros(Float64, Nx, Ny)
        for j in 1:Ny, i in 1:Nx
            f[i, j] = sin(kx_test * kgrid.x_vec[i])
        end

        # Expected: df/dx = kx * cos(kx * x)
        expected = zeros(Float64, Nx, Ny)
        for j in 1:Ny, i in 1:Nx
            expected[i, j] = kx_test * cos(kx_test * kgrid.x_vec[i])
        end

        # Compute spectral gradient (without staggered shift)
        df = similar(f)
        # Use identity shift (no staggering) for this test
        no_shift = ones(ComplexF64, Nx)
        KWave.spectral_gradient!(df, f, kgrid.kx_vec, no_shift, scratch, plans, 1, 2)

        # Should match analytical gradient (away from boundaries)
        interior = 5:Nx-4
        @test df[interior, Ny÷2] ≈ expected[interior, Ny÷2] atol=1e-10
    end
end
