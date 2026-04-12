@testset "1D Solver" begin
    @testset "Basic IVP simulation runs" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        # Gaussian initial pressure
        p0 = zeros(Float64, Nx)
        p0[Nx÷2-5:Nx÷2+5] .= 1.0
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx)
        sensor_mask[Nx÷2] = true
        sensor_mask[Nx÷2+20] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test output isa SimulationOutput
        @test haskey(output, :p)
        @test size(output[:p]) == (2, kgrid.Nt[])
        @test abs(output[:p][1, 1]) > 0
    end

    @testset "Energy stability" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=2e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0 = zeros(Float64, Nx)
        p0[Nx÷2-2:Nx÷2+2] .= 1.0
        source = KWaveSource(p0=p0)

        sensor = KWaveSensor(mask=trues(Nx), record=[:p])
        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        Nt = kgrid.Nt[]
        energy_start = sum(output[:p][:, 1].^2)
        energy_end = sum(output[:p][:, Nt].^2)

        @test energy_end < energy_start * 10
        @test all(isfinite.(output[:p]))
    end

    @testset "Wave propagation symmetry" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=2e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0 = zeros(Float64, Nx)
        p0[Nx÷2] = 1.0
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx)
        sensor_mask[Nx÷2+10] = true  # right
        sensor_mask[Nx÷2-10] = true  # left
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        p_data = output[:p]
        @test size(p_data, 1) == 2

        # Both sensors should see similar amplitude
        mid_time = kgrid.Nt[] ÷ 2
        @test abs(abs(p_data[1, mid_time]) - abs(p_data[2, mid_time])) < 0.5 * max(abs(p_data[1, mid_time]), abs(p_data[2, mid_time]), 1e-15)
    end

    @testset "Velocity source 1D" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)
        Nt = kgrid.Nt[]

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        source_mask = falses(Nx)
        source_mask[Nx÷2] = true
        signal = sin.(2π * 1e6 .* kgrid.t_array)'
        source = KWaveSource(u_mask=source_mask, ux=signal)

        sensor_mask = falses(Nx)
        sensor_mask[Nx÷2+20] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor)
        @test all(isfinite.(output[:p]))
        @test maximum(abs, output[:p]) > 0
    end

    @testset "Pressure source 1D" begin
        Nx = 128
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=5e-6)
        Nt = kgrid.Nt[]

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        source_mask = falses(Nx)
        source_mask[Nx÷2] = true
        signal = sin.(2π * 1e6 .* kgrid.t_array)'
        source = KWaveSource(p_mask=source_mask, p=signal)

        sensor_mask = falses(Nx)
        sensor_mask[Nx÷2+20] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor)
        @test all(isfinite.(output[:p]))
        @test maximum(abs, output[:p]) > 0
    end

    @testset "p_max and p_rms recording" begin
        Nx = 64
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx)
        make_time!(kgrid, 1500.0; t_end=1e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = zeros(Float64, Nx)
        p0[Nx÷2-2:Nx÷2+2] .= 1.0
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx)
        sensor_mask[Nx÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max, :p_rms])

        output = kspace_first_order(kgrid, medium, source, sensor; smooth_p0=false)

        @test haskey(output, :p_max)
        @test haskey(output, :p_rms)
        @test output[:p_max][1] >= maximum(output[:p][1, :]) - 1e-10
        @test output[:p_rms][1] > 0
    end
end
