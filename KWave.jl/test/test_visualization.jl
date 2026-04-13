@testset "Visualization (Phase 3)" begin
    @testset "plot function interfaces exist" begin
        # Without a Makie backend, these should warn and return nothing
        @test beam_plot(rand(10, 10)) === nothing
        @test fly_through(rand(10, 10, 10)) === nothing
        @test overlay_plot(rand(10, 10), rand(10, 10)) === nothing
        @test stacked_plot(rand(5, 100)) === nothing
    end

    @testset "get_color_map" begin
        cmap = get_color_map()
        @test length(cmap) == 256

        cmap64 = get_color_map(num_colors=64)
        @test length(cmap64) == 64

        # Check center is white
        r, g, b = cmap[128]
        @test r > 0.9
        @test g > 0.9
        @test b > 0.9
    end
end
