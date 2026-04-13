@testset "Python Interop Utilities" begin
    @testset "python_kwave_grid" begin
        # 1D
        g1 = python_kwave_grid(64, 1e-4)
        @test g1 isa KWaveGrid1D
        @test g1.Nx == 64

        # 2D
        g2 = python_kwave_grid(32, 1e-4, 32, 1e-4)
        @test g2 isa KWaveGrid2D
        @test g2.Nx == 32

        # 3D
        g3 = python_kwave_grid(16, 1e-4, 16, 1e-4, 16, 1e-4)
        @test g3 isa KWaveGrid3D
    end

    @testset "python_medium" begin
        m = python_medium(sound_speed=1500.0, density=1000.0)
        @test m isa KWaveMedium
        @test m.sound_speed == 1500.0
        @test m.density == 1000.0

        m2 = python_medium(sound_speed=1500.0, BonA=5.0)
        @test m2.BonA == 5.0
    end

    @testset "python_source" begin
        s = python_source()
        @test s isa KWaveSource
        @test s.p0 === nothing

        p0 = zeros(32, 32)
        p0[16, 16] = 1.0
        s2 = python_source(p0=p0)
        @test s2.p0 !== nothing
    end

    @testset "python_sensor" begin
        s = python_sensor()
        @test s isa KWaveSensor
        @test s.record == [:p]

        mask = falses(32, 32)
        mask[1, :] .= true
        s2 = python_sensor(mask=mask, record=[:p, :p_max])
        @test :p_max in s2.record
    end

    @testset "python_output_to_dict" begin
        data = Dict{Symbol, AbstractArray}(:p => ones(5, 10), :p_max => ones(5))
        output = SimulationOutput(data)

        d = python_output_to_dict(output)
        @test d isa Dict{String}
        @test haskey(d, "p")
        @test haskey(d, "p_max")
    end
end
