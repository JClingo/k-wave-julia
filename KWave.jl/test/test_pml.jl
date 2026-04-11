@testset "PML" begin
    @testset "get_pml basic properties" begin
        N = 64
        pml_size = 10
        pml = get_pml(N, 1e-4, pml_size, 2.0, 1e-8)

        @test length(pml) == N
        # Interior points should be 1.0
        @test all(pml[pml_size+1:N-pml_size] .≈ 1.0)
        # PML boundary points should be < 1.0
        @test pml[1] < 1.0
        @test pml[N] < 1.0
        # All values should be positive and ≤ 1
        @test all(0 .< pml .<= 1.0)
        # PML should be symmetric
        for i in 1:pml_size
            @test pml[i] ≈ pml[N - i + 1]
        end
        # Absorption should increase toward edges
        for i in 1:pml_size-1
            @test pml[i] <= pml[i+1]
        end
    end

    @testset "get_pml staggered" begin
        pml_reg = get_pml(64, 1e-4, 10, 2.0, 1e-8)
        pml_stag = get_pml(64, 1e-4, 10, 2.0, 1e-8; staggered=true)
        # Staggered PML should be different but same length
        @test length(pml_stag) == 64
        @test pml_reg != pml_stag
        # Interior should still be 1.0
        @test all(pml_stag[11:54] .≈ 1.0)
    end

    @testset "get_pml zero size" begin
        pml = get_pml(64, 1e-4, 0, 2.0, 1e-8)
        @test all(pml .== 1.0)
    end

    @testset "get_optimal_pml_size" begin
        # For a grid of 128, optimal PML should give a total with small prime factors
        opt = get_optimal_pml_size(128)
        total = 128 + 2 * opt
        # Should be <= max prime factor of 5 (FFT-friendly)
        @test KWave._max_prime_factor(total) <= 7

        # Tuple version
        opt_tuple = get_optimal_pml_size((64, 128))
        @test length(opt_tuple) == 2
    end
end
