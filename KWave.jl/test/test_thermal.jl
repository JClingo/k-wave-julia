@testset "Thermal Simulation" begin
    @testset "ThermalMedium construction" begin
        tm = ThermalMedium(
            thermal_conductivity=0.5,
            density=1050.0,
            specific_heat=3600.0,
        )
        @test tm.thermal_conductivity == 0.5
        @test tm.density == 1050.0
        @test tm.specific_heat == 3600.0
        @test tm.perfusion_rate === nothing
        @test tm.blood_temperature == 37.0

        # With perfusion
        tm2 = ThermalMedium(
            thermal_conductivity=0.5,
            density=1050.0,
            specific_heat=3600.0,
            perfusion_rate=0.5,
            metabolic_rate=420.0,
        )
        @test tm2.perfusion_rate == 0.5
        @test tm2.metabolic_rate == 420.0
    end

    @testset "ThermalSource construction" begin
        ts = ThermalSource()
        @test ts.Q === nothing
        @test ts.mode == :steady

        Q = ones(32)
        ts2 = ThermalSource(Q=Q, mode=:steady)
        @test ts2.Q !== nothing
    end

    @testset "kwave_diffusion 1D basic" begin
        Nx = 64; dx = 1e-3
        kgrid = KWaveGrid(Nx, dx)
        kgrid.dt[] = 0.01
        kgrid.Nt[] = 100

        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
        )

        source = ThermalSource()
        T0 = 37.0

        T_final, T_hist = kwave_diffusion(kgrid, medium, source, T0, 100)

        # Without any source, temperature should stay at initial value
        @test all(isapprox.(T_final, 37.0, atol=1e-6))
        @test size(T_hist, 2) >= 2
    end

    @testset "kwave_diffusion 1D with source" begin
        Nx = 64; dx = 1e-3
        kgrid = KWaveGrid(Nx, dx)
        kgrid.dt[] = 0.001
        kgrid.Nt[] = 100

        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
        )

        Q = zeros(Nx)
        Q[Nx÷2] = 1e6
        source = ThermalSource(Q=Q)

        T_final, _ = kwave_diffusion(kgrid, medium, source, 37.0, 100)

        # Temperature at center should increase
        @test T_final[Nx÷2] > 37.0
        # Temperature at edges should stay near initial
        @test T_final[1] ≈ 37.0 atol=0.1
    end

    @testset "kwave_diffusion 1D with perfusion" begin
        Nx = 64; dx = 1e-3
        kgrid = KWaveGrid(Nx, dx)
        kgrid.dt[] = 0.01
        kgrid.Nt[] = 100

        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
            perfusion_rate=1.0,
        )

        # Start above blood temperature — perfusion should cool it down
        T0 = fill(45.0, Nx)
        source = ThermalSource()

        T_final, _ = kwave_diffusion(kgrid, medium, source, T0, 100)

        # Should cool toward blood temperature (37°C)
        @test maximum(T_final) < 45.0
    end

    @testset "kwave_diffusion 2D" begin
        Nx, Ny = 32, 32; dx, dy = 1e-3, 1e-3
        kgrid = KWaveGrid(Nx, dx, Ny, dy)
        kgrid.dt[] = 0.001
        kgrid.Nt[] = 50

        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
        )

        Q = zeros(Nx, Ny)
        Q[Nx÷2, Ny÷2] = 1e6
        source = ThermalSource(Q=Q)

        T_final, _ = kwave_diffusion(kgrid, medium, source, 37.0, 50)

        @test T_final[Nx÷2, Ny÷2] > 37.0
        @test size(T_final) == (Nx, Ny)
    end

    @testset "kwave_diffusion 3D" begin
        Nx, Ny, Nz = 16, 16, 16; dx, dy, dz = 1e-3, 1e-3, 1e-3
        kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz)
        kgrid.dt[] = 0.001
        kgrid.Nt[] = 20

        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
        )

        Q = zeros(Nx, Ny, Nz)
        Q[Nx÷2, Ny÷2, Nz÷2] = 1e6
        source = ThermalSource(Q=Q)

        T_final, _ = kwave_diffusion(kgrid, medium, source, 37.0, 20)

        @test T_final[Nx÷2, Ny÷2, Nz÷2] > 37.0
        @test size(T_final) == (Nx, Ny, Nz)
    end

    @testset "bioheat_exact scalar" begin
        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
            perfusion_rate=0.5,
        )

        # At t=0, should return initial temperature
        T = bioheat_exact(37.0, 0.0, medium, 0.0)
        @test T ≈ 37.0

        # With source, temperature should increase
        T = bioheat_exact(37.0, 1e6, medium, 10.0)
        @test T > 37.0

        # Vector of times
        T_vec = bioheat_exact(37.0, 1e6, medium, [0.0, 5.0, 10.0])
        @test length(T_vec) == 3
        @test T_vec[1] ≈ 37.0
        @test T_vec[3] > T_vec[2] > T_vec[1]
    end

    @testset "bioheat_exact no perfusion" begin
        medium = ThermalMedium(
            thermal_conductivity=0.5,
            density=1000.0,
            specific_heat=4000.0,
        )

        # Linear temperature rise with constant source
        T = bioheat_exact(37.0, 4e6, medium, 1.0)
        expected = 37.0 + 4e6 / (1000.0 * 4000.0) * 1.0
        @test T ≈ expected
    end
end
