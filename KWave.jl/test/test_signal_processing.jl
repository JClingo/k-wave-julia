@testset "Signal Processing (Phase 3)" begin
    @testset "add_noise" begin
        signal = sin.(2π .* (0:999) ./ 100)
        noisy = add_noise(signal, 20.0)
        @test length(noisy) == length(signal)
        # Signal should be changed
        @test noisy != signal
    end

    @testset "create_cw_signals" begin
        kgrid = KWaveGrid(128, 1e-4)
        make_time!(kgrid, 1500.0; t_end=1e-5)
        signals = create_cw_signals(kgrid, 1e6, 1.0, 0.0)
        @test size(signals) == (1, kgrid.Nt[])
    end

    @testset "create_cw_signals with ramp" begin
        kgrid = KWaveGrid(128, 1e-4)
        make_time!(kgrid, 1500.0; t_end=1e-5)
        signals = create_cw_signals(kgrid, 1e6, [1.0, 2.0], [0.0, π/4]; num_cycles=3)
        @test size(signals, 1) == 2
        # Signal should start near zero (ramped)
        @test abs(signals[1, 1]) < 0.1
    end

    @testset "log_compression" begin
        signal = rand(100)
        compressed = log_compression(signal, 10.0)
        @test all(0 .<= compressed .<= 1)
    end

    @testset "envelope_detection" begin
        t = (0:999) ./ 1000.0
        signal = sin.(2π * 10 .* t) .* exp.(-5 .* (t .- 0.5).^2)
        env = envelope_detection(signal)
        @test length(env) == length(signal)
        # Envelope should be non-negative
        @test all(real.(env) .>= -1e-10)
    end

    @testset "gradient_fd" begin
        x = collect(0.0:0.1:10.0)
        f = x.^2
        df = gradient_fd(f, 0.1)
        # At x=5 (index 51), gradient should be approximately 2*5=10
        @test isapprox(df[51], 10.0; atol=0.1)
    end

    @testset "gradient_spect" begin
        kgrid = KWaveGrid(64, 2π/64)
        f = sin.(kgrid.x_vec)
        df = gradient_spect(f, kgrid)
        # Derivative of sin(x) = cos(x)
        expected = cos.(kgrid.x_vec)
        @test isapprox(df, expected; atol=1e-10)
    end

    @testset "spect" begin
        Fs = 1000.0
        t = (0:999) ./ Fs
        f_sig = 100.0
        signal = sin.(2π * f_sig .* t)
        freqs, amp, phase = spect(signal, Fs)
        # Peak should be at f_sig
        peak_idx = argmax(amp[2:end]) + 1  # skip DC
        @test isapprox(freqs[peak_idx], f_sig; atol=2.0)
    end
end
