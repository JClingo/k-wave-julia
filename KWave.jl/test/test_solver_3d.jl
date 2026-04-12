@testset "3D Solver" begin
    @testset "Basic IVP simulation runs" begin
        Nx, Ny, Nz = 32, 32, 32
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx, Nz, dx)
        make_time!(kgrid, 1500.0; t_end=2e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0 = Float64.(make_ball(Nx, Ny, Nz, Nx÷2, Ny÷2, Nz÷2, 3))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny, Nz)
        sensor_mask[Nx÷2, Ny÷2, Nz÷2] = true
        sensor_mask[Nx÷2+8, Ny÷2, Nz÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        output = kspace_first_order(kgrid, medium, source, sensor;
                                     smooth_p0=false, pml_size=10)

        @test output isa SimulationOutput
        @test haskey(output, :p)
        @test size(output[:p]) == (2, kgrid.Nt[])
        @test abs(output[:p][1, 1]) > 0
    end

    @testset "Energy stability 3D" begin
        Nx, Ny, Nz = 24, 24, 24
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx, Nz, dx)
        make_time!(kgrid, 1500.0; t_end=1e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

        p0 = Float64.(make_ball(Nx, Ny, Nz, Nx÷2, Ny÷2, Nz÷2, 2))
        source = KWaveSource(p0=p0)

        sensor = KWaveSensor(mask=trues(Nx, Ny, Nz), record=[:p])
        output = kspace_first_order(kgrid, medium, source, sensor;
                                     smooth_p0=false, pml_size=8)

        Nt = kgrid.Nt[]
        energy_start = sum(output[:p][:, 1].^2)
        energy_end = sum(output[:p][:, Nt].^2)

        @test energy_end < energy_start * 10
        @test all(isfinite.(output[:p]))
    end

    @testset "p_max recording 3D" begin
        Nx, Ny, Nz = 24, 24, 24
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx, Nz, dx)
        make_time!(kgrid, 1500.0; t_end=1e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_ball(Nx, Ny, Nz, Nx÷2, Ny÷2, Nz÷2, 2))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny, Nz)
        sensor_mask[Nx÷2, Ny÷2, Nz÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p, :p_max])

        output = kspace_first_order(kgrid, medium, source, sensor;
                                     smooth_p0=false, pml_size=8)

        @test haskey(output, :p_max)
        @test output[:p_max][1] >= maximum(output[:p][1, :]) - 1e-10
    end

    @testset "Velocity recording 3D" begin
        Nx, Ny, Nz = 24, 24, 24
        dx = 1e-4

        kgrid = KWaveGrid(Nx, dx, Ny, dx, Nz, dx)
        make_time!(kgrid, 1500.0; t_end=1e-6)

        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        p0 = Float64.(make_ball(Nx, Ny, Nz, Nx÷2, Ny÷2, Nz÷2, 2))
        source = KWaveSource(p0=p0)

        sensor_mask = falses(Nx, Ny, Nz)
        sensor_mask[Nx÷2, Ny÷2, Nz÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p, :ux, :uy, :uz, :u_max])

        output = kspace_first_order(kgrid, medium, source, sensor;
                                     smooth_p0=false, pml_size=8)

        @test haskey(output, :ux)
        @test haskey(output, :uy)
        @test haskey(output, :uz)
        @test haskey(output, :u_max)
        @test all(isfinite.(output[:ux]))
        @test output[:u_max][1] >= 0
    end
end
