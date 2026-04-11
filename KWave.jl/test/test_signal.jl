@testset "Signal Generation" begin
    @testset "tone_burst" begin
        Fs = 40e6
        f0 = 1e6
        cycles = 5
        signal = tone_burst(Fs, f0, cycles)

        @test signal isa Vector{Float64}
        @test length(signal) > 0

        # Signal should start and end near zero (windowed)
        @test abs(signal[1]) < 0.1
        @test abs(signal[end]) < 0.1

        # Maximum should be near 1.0 (unit amplitude carrier with envelope)
        @test maximum(abs.(signal)) > 0.5
        @test maximum(abs.(signal)) <= 1.0 + 1e-10

        # With offset
        signal_offset = tone_burst(Fs, f0, cycles; signal_offset=100)
        @test length(signal_offset) == length(signal) + 100
        @test all(signal_offset[1:100] .== 0.0)
    end

    @testset "gaussian_pulse" begin
        # Scalar
        @test gaussian_pulse(0.0) ≈ 1.0
        @test gaussian_pulse(0.0; magnitude=2.0) ≈ 2.0
        @test gaussian_pulse(1.0; mean=1.0) ≈ 1.0

        # Vector
        x = collect(-3:0.1:3)
        g = gaussian_pulse(x)
        @test length(g) == length(x)
        @test g[31] ≈ 1.0 atol=1e-10  # x = 0
        @test all(g .>= 0)
        @test all(g .<= 1.0 + 1e-10)

        # Symmetry
        @test gaussian_pulse(-1.0) ≈ gaussian_pulse(1.0)
    end
end
