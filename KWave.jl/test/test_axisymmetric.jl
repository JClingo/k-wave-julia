@testset "Axisymmetric Solver (Phase 3)" begin
    @testset "basic run" begin
        # Create a 2D grid (axial × radial)
        Nx, Nr = 64, 32
        dx = 1e-4
        kgrid = KWaveGrid(Nx, dx, Nr, dx)

        # Homogeneous medium
        medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
        make_time!(kgrid, medium.sound_speed; t_end=2e-6)

        # Initial pressure: Gaussian blob
        p0 = zeros(Nx, Nr)
        for j in 1:Nr, i in 1:Nx
            r2 = ((i - Nx÷2) * dx)^2 + ((j - Nr÷2) * dx)^2
            p0[i, j] = exp(-r2 / (5 * dx)^2)
        end

        source = KWaveSource(p0=p0)
        sensor_mask = falses(Nx, Nr)
        sensor_mask[Nx÷2, :] .= true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        result = kspace_first_order_as(kgrid, medium, source, sensor)
        @test haskey(result, :p)
        @test size(result[:p], 2) == kgrid.Nt[]
    end

    @testset "with absorption" begin
        Nx, Nr = 64, 64
        dx = 1e-4
        kgrid = KWaveGrid(Nx, dx, Nr, dx)
        medium = KWaveMedium(
            sound_speed=1500.0, density=1000.0,
            alpha_coeff=0.5, alpha_power=1.5,
            alpha_mode=:no_dispersion
        )
        make_time!(kgrid, medium.sound_speed; t_end=1e-6)

        p0 = zeros(Nx, Nr)
        p0[Nx÷2, Nr÷2] = 1.0
        source = KWaveSource(p0=p0)
        sensor_mask = falses(Nx, Nr)
        sensor_mask[Nx÷2, 1] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p_max])

        result = kspace_first_order_as(kgrid, medium, source, sensor; pml_size=10)
        @test haskey(result, :p_max)
    end
end
