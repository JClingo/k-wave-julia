@testset "Absorption and Nonlinear" begin
    @testset "Power-law db2neper conversion" begin
        # Test two-argument form
        alpha_coeff = 0.75  # dB/(MHz^y cm)
        alpha_power = 1.5
        alpha_nepers = db2neper(alpha_coeff, alpha_power)
        @test alpha_nepers > 0
        @test isfinite(alpha_nepers)

        # Verify dimensional consistency
        # Higher power should give smaller neper coefficient (at standard frequencies)
        a1 = db2neper(1.0, 1.0)
        a2 = db2neper(1.0, 2.0)
        @test a1 > a2
    end

    @testset "Absorbing 2D simulation runs without error" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(
            sound_speed=1500.0,
            density=1000.0,
            alpha_coeff=0.75,
            alpha_power=1.5,
            alpha_mode=:no_dispersion,
        )

        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2, Ny÷2] = true
        sensor_mask[Nx÷2+15, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test output isa SimulationOutput
        @test all(isfinite.(output[:p]))
    end

    @testset "Absorption reduces total field energy" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        # Lossless run
        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium_lossless = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5))
        source = KWaveSource(p0=p0)

        # Record entire field to compare total energy
        sensor = KWaveSensor(mask=trues(Nx, Ny), record=[:p_rms])

        out_lossless = kspace_first_order(kgrid, medium_lossless, source, sensor; smooth_p0=false)

        # Absorbing run
        kgrid2 = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid2, 1500.0; t_end=5e-6)

        medium_absorbing = KWaveMedium(
            sound_speed=1500.0,
            density=1000.0,
            alpha_coeff=2.0,
            alpha_power=1.5,
            alpha_mode=:no_dispersion,
        )

        source2 = KWaveSource(p0=p0)
        out_absorbing = kspace_first_order(kgrid2, medium_absorbing, source2, sensor; smooth_p0=false)

        # Absorbing should have lower total RMS energy
        energy_lossless = sum(out_lossless[:p_rms].^2)
        energy_absorbing = sum(out_absorbing[:p_rms].^2)
        @test energy_absorbing < energy_lossless
    end

    @testset "Absorbing with dispersion" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(
            sound_speed=1500.0,
            density=1000.0,
            alpha_coeff=0.75,
            alpha_power=1.5,
            alpha_mode=:stokes,
        )

        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2+15, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test all(isfinite.(output[:p]))
    end

    @testset "Absorbing 1D" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(
            sound_speed=1500.0,
            density=1000.0,
            alpha_coeff=0.75,
            alpha_power=1.5,
            alpha_mode=:no_dispersion,
        )

        p0 = zeros(Float64, Nx)
        p0[Nx÷2-3:Nx÷2+3] .= 1.0
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx)
        sensor_mask[Nx÷2+20] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        @test all(isfinite.(output[:p]))
    end

    @testset "Nonlinear 2D simulation" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(
            sound_speed=1500.0,
            density=1000.0,
            BonA=5.0,  # typical water-like nonlinearity
        )

        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5)) .* 1e6  # high amplitude for nonlinear effects
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2, Ny÷2] = true
        sensor_mask[Nx÷2+15, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test output isa SimulationOutput
        @test all(isfinite.(output[:p]))
    end

    @testset "Nonlinear + absorbing combined" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(
            sound_speed=1500.0,
            density=1000.0,
            alpha_coeff=0.75,
            alpha_power=1.5,
            alpha_mode=:no_dispersion,
            BonA=5.0,
        )

        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5)) .* 1e6
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2+15, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test all(isfinite.(output[:p]))
    end

    @testset "Heterogeneous sound speed" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)

        # Heterogeneous sound speed: two-layer medium
        c0 = fill(1500.0, Nx, Ny)
        c0[Nx÷2+1:end, :] .= 1600.0
        make_time!(kgrid, c0; t_end=5e-6)

        medium = KWaveMedium(sound_speed=c0, density=1000.0)

        p0 = Float64.(make_disc(Nx, Ny, Nx÷4, Ny÷2, 5))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[3*Nx÷4, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test all(isfinite.(output[:p]))
        @test maximum(abs, output[:p]) > 0
    end

    @testset "Heterogeneous density" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        # Heterogeneous density
        rho0 = fill(1000.0, Nx, Ny)
        rho0[Nx÷2+1:end, :] .= 1200.0

        medium = KWaveMedium(sound_speed=1500.0, density=rho0)

        p0 = Float64.(make_disc(Nx, Ny, Nx÷4, Ny÷2, 5))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[3*Nx÷4, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test all(isfinite.(output[:p]))
    end
end
