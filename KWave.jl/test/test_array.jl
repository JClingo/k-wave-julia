@testset "KWaveArray (Phase 3)" begin
    @testset "construction" begin
        array = KWaveArray()
        @test length(array) == 0

        add_arc_element!(array, (0.0, 0.0), 0.01, 0.005, (0.01, 0.0))
        @test length(array) == 1
        @test array.elements[1] isa ArcElement
    end

    @testset "add elements" begin
        array = KWaveArray()
        add_arc_element!(array, (0.0, 0.0), 0.01, 0.005, (0.01, 0.0))
        add_bowl_element!(array, (0.0, 0.0, 0.0), 0.01, 0.005, (0.01, 0.0, 0.0))
        add_disc_element!(array, (0.0, 0.0, 0.0), 0.005, (0.01, 0.0, 0.0))
        add_rect_element!(array, (0.0, 0.0, 0.0), 0.005, 0.003, (0.01, 0.0, 0.0))
        @test length(array) == 4
    end

    @testset "get_element_binary_mask 2D" begin
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        array = KWaveArray()
        add_arc_element!(array, (0.0, 0.0), 15 * 1e-4, 20 * 1e-4, (15 * 1e-4, 0.0))
        mask = get_element_binary_mask(array, kgrid, 1)
        @test size(mask) == (64, 64)
        @test count(mask) > 0
    end

    @testset "get_array_binary_mask" begin
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        array = KWaveArray()
        add_arc_element!(array, (0.0, 0.0), 15 * 1e-4, 10 * 1e-4, (15 * 1e-4, 0.0))
        add_arc_element!(array, (0.0, 0.0), 15 * 1e-4, 10 * 1e-4, (-15 * 1e-4, 0.0))
        mask = get_array_binary_mask(array, kgrid)
        @test count(mask) > 0
    end

    @testset "get_distributed_source_signal" begin
        kgrid = KWaveGrid(64, 1e-4, 64, 1e-4)
        make_time!(kgrid, 1500.0; t_end=1e-5)
        Nt = kgrid.Nt[]

        array = KWaveArray()
        add_arc_element!(array, (0.0, 0.0), 15 * 1e-4, 10 * 1e-4, (15 * 1e-4, 0.0))

        source_signal = sin.(2π .* 1e6 .* kgrid.t_array')  # 1 × Nt
        mask, distributed = get_distributed_source_signal(array, kgrid, source_signal)
        @test count(mask) > 0
        @test size(distributed, 2) == Nt
    end
end
