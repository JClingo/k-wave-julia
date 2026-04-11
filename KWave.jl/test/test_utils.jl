@testset "Grid Utilities" begin
    @testset "expand_matrix" begin
        m = ones(4, 4)
        expanded = expand_matrix(m, 2)
        @test size(expanded) == (8, 8)
        # Center should be ones
        @test all(expanded[3:6, 3:6] .== 1.0)
        # Edges should be zero (default)
        @test expanded[1, 1] == 0.0

        # Non-uniform expansion
        expanded2 = expand_matrix(m, ((1, 2), (3, 4)))
        @test size(expanded2) == (4 + 1 + 2, 4 + 3 + 4)

        # Custom edge value
        expanded3 = expand_matrix(m, 1; edge_val=5.0)
        @test expanded3[1, 1] == 5.0
    end

    @testset "resize_array 1D" begin
        v = Float64[1, 2, 3, 4, 5]
        resized = resize_array(v, (9,))
        @test length(resized) == 9
        @test resized[1] ≈ 1.0
        @test resized[end] ≈ 5.0

        # Identity resize
        same = resize_array(v, (5,))
        @test same ≈ v
    end

    @testset "resize_array 2D" begin
        m = [1.0 2.0; 3.0 4.0]
        resized = resize_array(m, (4, 4))
        @test size(resized) == (4, 4)
        # Corners should preserve original values
        @test resized[1, 1] ≈ 1.0
        @test resized[4, 4] ≈ 4.0
    end

    @testset "cart2grid / grid2cart round-trip" begin
        kgrid = KWaveGrid(32, 1e-4, 32, 1e-4)

        # Create some Cartesian points
        cart = [kgrid.x_vec[10] kgrid.x_vec[20];
                kgrid.y_vec[15] kgrid.y_vec[25]]

        mask, order_idx, reorder_idx = cart2grid(kgrid, cart)
        @test mask isa BitMatrix
        @test size(mask) == (32, 32)
        @test count(mask) == 2
        @test mask[10, 15] == true
        @test mask[20, 25] == true

        # Convert back
        cart_back = grid2cart(kgrid, mask)
        @test size(cart_back, 2) == 2
    end
end
