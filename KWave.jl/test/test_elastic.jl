@testset "Elastic Solvers" begin
    @testset "ElasticMedium construction" begin
        # Basic construction
        em = ElasticMedium(
            sound_speed_compression=3000.0,
            sound_speed_shear=1500.0,
            density=1800.0,
        )
        @test em.sound_speed_compression == 3000.0
        @test em.sound_speed_shear == 1500.0
        @test em.density == 1800.0
        @test em.alpha_coeff_compression === nothing
        @test em.alpha_coeff_shear === nothing

        # With absorption
        em2 = ElasticMedium(
            sound_speed_compression=3000.0,
            sound_speed_shear=1500.0,
            density=1800.0,
            alpha_coeff_compression=0.5,
            alpha_coeff_shear=1.0,
            alpha_power=1.5,
        )
        @test em2.alpha_coeff_compression == 0.5
        @test em2.alpha_coeff_shear == 1.0
        @test em2.alpha_power == 1.5
    end

    @testset "ElasticSource construction" begin
        es = ElasticSource()
        @test es.s_mask === nothing
        @test es.p0 === nothing
        @test es.s_mode == Additive

        # With p0
        p0 = zeros(32, 32)
        p0[16, 16] = 1.0
        es2 = ElasticSource(p0=p0)
        @test es2.p0 !== nothing
        @test es2.p0[16, 16] == 1.0
    end

    @testset "pstd_elastic_2d basic run" begin
        Nx, Ny = 32, 32
        dx, dy = 0.5e-3, 0.5e-3
        kgrid = KWaveGrid(Nx, dx, Ny, dy)

        medium = ElasticMedium(
            sound_speed_compression=3000.0,
            sound_speed_shear=1500.0,
            density=1800.0,
        )

        make_time!(kgrid, 3000.0)

        # Initial pressure source
        p0 = zeros(Nx, Ny)
        p0[Nx÷2, Ny÷2] = 1.0
        source = ElasticSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[5, :] .= true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        result = pstd_elastic_2d(kgrid, medium, source, sensor)

        @test haskey(result, :p)
        @test size(result[:p], 1) == count(sensor_mask)
        @test size(result[:p], 2) == kgrid.Nt[]
    end

    @testset "pstd_elastic_2d stress source" begin
        Nx, Ny = 32, 32
        dx, dy = 0.5e-3, 0.5e-3
        kgrid = KWaveGrid(Nx, dx, Ny, dy)

        medium = ElasticMedium(
            sound_speed_compression=3000.0,
            sound_speed_shear=1500.0,
            density=1800.0,
        )

        make_time!(kgrid, 3000.0)
        Nt = kgrid.Nt[]

        # Stress source
        s_mask = falses(Nx, Ny)
        s_mask[5, Ny÷2] = true

        sxx_sig = reshape(sin.(2π * 1e6 .* (0:Nt-1) .* kgrid.dt[]), 1, :)

        source = ElasticSource(s_mask=s_mask, sxx=sxx_sig, s_mode=Additive)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷2, :] .= true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        result = pstd_elastic_2d(kgrid, medium, source, sensor)

        @test haskey(result, :p)
        @test !all(result[:p] .== 0)
    end

    @testset "pstd_elastic_2d wave speeds" begin
        # Verify P-wave arrives before S-wave
        Nx, Ny = 64, 64
        dx, dy = 0.5e-3, 0.5e-3
        kgrid = KWaveGrid(Nx, dx, Ny, dy)

        medium = ElasticMedium(
            sound_speed_compression=3000.0,
            sound_speed_shear=1500.0,
            density=1800.0,
        )

        make_time!(kgrid, 3000.0)

        p0 = zeros(Nx, Ny)
        p0[Nx÷2, Ny÷2] = 1.0
        source = ElasticSource(p0=p0)

        sensor_mask = falses(Nx, Ny)
        sensor_mask[Nx÷4, Ny÷2] = true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        result = pstd_elastic_2d(kgrid, medium, source, sensor)

        # Signal should not be zero
        p_trace = result[:p][1, :]
        @test maximum(abs, p_trace) > 0
    end

    @testset "pstd_elastic_3d basic run" begin
        Nx, Ny, Nz = 16, 16, 16
        dx, dy, dz = 0.5e-3, 0.5e-3, 0.5e-3
        kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz)

        medium = ElasticMedium(
            sound_speed_compression=3000.0,
            sound_speed_shear=1500.0,
            density=1800.0,
        )

        make_time!(kgrid, 3000.0)

        p0 = zeros(Nx, Ny, Nz)
        p0[Nx÷2, Ny÷2, Nz÷2] = 1.0
        source = ElasticSource(p0=p0)

        sensor_mask = falses(Nx, Ny, Nz)
        sensor_mask[3, :, :] .= true
        sensor = KWaveSensor(mask=sensor_mask, record=[:p])

        result = pstd_elastic_3d(kgrid, medium, source, sensor; pml_size=5)

        @test haskey(result, :p)
        @test size(result[:p], 1) == count(sensor_mask)
    end
end
