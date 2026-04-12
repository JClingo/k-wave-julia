@testset "Time Reversal" begin
    @testset "2D time-reversal reconstruction" begin
        Nx, Ny = 64, 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)
        Nt = kgrid.Nt[]

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0_original = Float64.(make_disc(Nx, Ny, Nx÷2, Ny÷2, 3))
        source = KWaveSource(p0=p0_original)

        # Sensor around the boundary
        sensor_mask = falses(Nx, Ny)
        sensor_mask[1, :] .= true
        sensor_mask[Nx, :] .= true
        sensor_mask[:, 1] .= true
        sensor_mask[:, Ny] .= true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        # Forward run
        output_fwd = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)
        recorded_p = output_fwd[:p]

        # Time-reversal reconstruction
        # sensor.mask = injection mask (boundary), p_final = full field output
        kgrid_tr = KWaveGrid(Nx, dx, Ny, dx)
        make_time!(kgrid_tr, 1500.0; t_end=5e-6)

        sensor_tr = KWaveSensor(
            mask=sensor_mask,
            time_reversal_boundary_data=recorded_p,
            record=[:p_final],
        )
        source_tr = KWaveSource()

        output_tr = kspace_first_order(kgrid_tr, medium, source_tr, sensor_tr; smooth_p0=false)

        @test haskey(output_tr, :p_final)
        p_recon = output_tr[:p_final]

        # Full-field output: reshaped to grid
        p_grid = reshape(p_recon, Nx, Ny)

        # Find peak location in reconstructed field
        peak_idx = argmax(abs.(p_grid))
        peak_i, peak_j = Tuple(peak_idx)

        # Peak should be near center (within a few grid points)
        @test abs(peak_i - Nx÷2) <= 5
        @test abs(peak_j - Ny÷2) <= 5
    end

    @testset "1D time-reversal" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0 = zeros(Float64, Nx)
        p0[Nx÷2-2:Nx÷2+2] .= 1.0
        source = KWaveSource(p0=p0)

        # Sensors at boundaries
        sensor_mask = falses(Nx)
        sensor_mask[1] = true
        sensor_mask[Nx] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output_fwd = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        # Time-reversal
        kgrid_tr = KWaveGrid(Nx, dx)
        make_time!(kgrid_tr, 1500.0; t_end=5e-6)

        sensor_tr = KWaveSensor(
            mask=sensor_mask,
            time_reversal_boundary_data=output_fwd[:p],
            record=[:p_final],
        )
        source_tr = KWaveSource()

        output_tr = kspace_first_order(kgrid_tr, medium, source_tr, sensor_tr; smooth_p0=false)

        @test haskey(output_tr, :p_final)
        p_recon = output_tr[:p_final]
        @test all(isfinite.(p_recon))

        # Peak should be near center
        peak_idx = argmax(abs.(p_recon))
        @test abs(peak_idx - Nx÷2) <= 5
    end
end
