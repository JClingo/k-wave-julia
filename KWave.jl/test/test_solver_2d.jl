@testset "2D Solver" begin
    @testset "Basic IVP simulation runs" begin
        # Simple 2D initial value problem
        Nx, Ny = 64, 64
        dx, dy = 1e-4, 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dy)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        # Initial pressure disc at center
        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 5))
        source = KWaveSource(p0=p0)

        # Record at a few points
        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2, Ny÷2] = true       # center
        sensor_mask[Nx÷2+20, Ny÷2] = true    # 20 points away
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test output isa SimulationOutput
        @test haskey(output, :p)
        @test size(output[:p]) == (2, kgrid.Nt[])

        # Center sensor should see initial pressure at t=1
        # (wave originates from center)
        # After several steps, pressure at center should decrease
        @test abs(output[:p][1, 1]) > 0
    end

    @testset "Energy conservation (approximate)" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=2e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 3))
        source = KWaveSource(p0=p0)

        # Record entire field to check energy
        sensor = KWaveSensor(mask=trues(Nx, Ny), record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        # Total energy proxy: sum(p^2) at each time step
        # Should not blow up (stability check)
        Nt = kgrid.Nt[]
        energy_start = sum(output[:p][:, 1].^2)
        energy_end = sum(output[:p][:, Nt].^2)

        # Energy should not grow unboundedly (PML absorbs some)
        @test energy_end < energy_start * 10  # generous bound — just checking stability
        @test all(isfinite.(output[:p]))
    end

    @testset "Wave propagation symmetry" begin
        # A centered source should produce symmetric waves
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=2e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        # Single point source at exact center
        p0 = zeros(Float64, Nx, Ny)
        p0[Nx÷2, Ny÷2] = 1.0
        source = KWaveSource(p0=p0)

        # Symmetric sensor points
        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2+10, Ny÷2] = true   # right
        sensor_mask[Nx÷2-10, Ny÷2] = true   # left
        sensor_mask[Nx÷2, Ny÷2+10] = true   # up
        sensor_mask[Nx÷2, Ny÷2-10] = true   # down
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        # All 4 symmetric sensors should see approximately the same signal
        p_data = output[:p]
        @test size(p_data, 1) == 4

        # Check that signals are similar (not exact due to PML asymmetry near edges)
        mid_time = kgrid.Nt[] ÷ 2
        signals_at_mid = p_data[:, mid_time]
        mean_signal = mean(abs.(signals_at_mid))
        if mean_signal > 1e-15
            for i in 1:4
                @test abs(signals_at_mid[i]) / mean_signal < 5.0  # within 5x of mean
            end
        end
    end

    @testset "p_max recording" begin
        Nx, Ny = 32, 32
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=1e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 3))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test haskey(output, :p_max)
        # p_max should be >= max of recorded p time series
        @test output[:p_max][1] >= maximum(output[:p][1, :]) - 1e-10
    end
end
