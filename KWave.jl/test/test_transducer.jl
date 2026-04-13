@testset "KWaveTransducer (Phase 3)" begin
    @testset "construction" begin
        t = KWaveTransducer(
            number_elements=32,
            element_width=2,
            element_length=10,
        )
        @test t.number_elements == 32
        @test t.element_width == 2
        @test t.element_spacing == 0
        @test t.radius == Inf
        @test t.transmit_apodization == :rectangular
    end

    @testset "binary mask" begin
        kgrid = KWaveGrid(32, 1e-4, 64, 1e-4, 32, 1e-4)
        t = KWaveTransducer(
            number_elements=8,
            element_width=4,
            element_length=16,
            position=(16, 5, 8),
        )
        mask = get_transducer_binary_mask(t, kgrid)
        @test size(mask) == (32, 64, 32)
        @test count(mask) == 8 * 4 * 16  # 8 elements × 4 width × 16 length
    end

    @testset "binary mask with active elements" begin
        kgrid = KWaveGrid(32, 1e-4, 64, 1e-4, 32, 1e-4)
        t = KWaveTransducer(
            number_elements=8,
            element_width=4,
            element_length=16,
            position=(16, 5, 8),
            active_elements=[1, 3, 5],
        )
        mask = get_transducer_binary_mask(t, kgrid)
        @test count(mask) == 3 * 4 * 16
    end

    @testset "source generation" begin
        kgrid = KWaveGrid(32, 1e-4, 64, 1e-4, 32, 1e-4)
        make_time!(kgrid, 1500.0; t_end=2e-6)

        sig = sin.(2π * 1e6 .* kgrid.t_array)
        t = KWaveTransducer(
            number_elements=4,
            element_width=2,
            element_length=8,
            position=(16, 10, 12),
            input_signal=sig,
        )
        medium = KWaveMedium(sound_speed=1500.0)
        mask, source_signal = get_transducer_source(t, kgrid, medium)
        @test count(mask) > 0
        @test size(source_signal, 2) == kgrid.Nt[]
    end
end
