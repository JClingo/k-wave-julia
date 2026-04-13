@testset "CW Propagation (Phase 3)" begin
    @testset "acoustic_field_propagator 2D" begin
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        medium = KWaveMedium(sound_speed=1500.0)

        source_amp = zeros(64, 64)
        source_phase = zeros(64, 64)
        source_amp[32, 32] = 1.0  # Point source at center

        amp, phase = acoustic_field_propagator(kgrid, medium, source_amp, source_phase, 1e6)
        @test size(amp) == (64, 64)
        @test size(phase) == (64, 64)
        @test maximum(amp) > 0
    end

    @testset "angular_spectrum_cw 1D" begin
        N = 64
        dx = 1e-4
        freq = 1e6
        c0 = 1500.0

        # Gaussian source
        x = ((1:N) .- N÷2) .* dx
        input = exp.(-x.^2 ./ (5 * dx)^2)

        result = angular_spectrum_cw(input, dx, freq, c0, 10 * dx)
        @test length(result) == N
    end

    @testset "angular_spectrum_cw 2D" begin
        Nx, Ny = 32, 32
        dx, dy = 1e-4, 1e-4
        freq = 1e6
        c0 = 1500.0

        input = zeros(Nx, Ny)
        input[16, 16] = 1.0

        result = angular_spectrum_cw(input, dx, dy, freq, c0, 5 * dx)
        @test size(result) == (Nx, Ny)
    end

    @testset "angular_spectrum_cw multiple planes" begin
        N = 32
        dx = 1e-4
        input = zeros(N)
        input[16] = 1.0

        z_positions = [1e-4, 2e-4, 5e-4]
        result = angular_spectrum_cw(input, dx, 1e6, 1500.0, z_positions)
        @test size(result) == (N, 3)
    end
end
