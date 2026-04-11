@testset "Filtering" begin
    @testset "get_win" begin
        # Hann window
        w = get_win(64, :hann)
        @test length(w) == 64
        @test w[1] ≈ 0.0 atol=1e-10
        @test maximum(w) > 0.99  # near-unity peak for even-length window
        @test all(0 .<= w .<= 1.0 + 1e-10)

        # Hamming window
        w_ham = get_win(64, :hamming)
        @test length(w_ham) == 64
        @test w_ham[33] > 0.9  # peak near center

        # Blackman window
        w_blk = get_win(64, :blackman)
        @test length(w_blk) == 64
        @test w_blk[1] ≈ 0.0 atol=1e-10

        # Rectangular window
        w_rect = get_win(64, :rectangular)
        @test all(w_rect .== 1.0)

        # Kaiser window
        w_kai = get_win(64, :kaiser; param=5.0)
        @test length(w_kai) == 64
        @test all(w_kai .>= 0)

        # Edge cases
        @test get_win(1, :hann) == [1.0]
        @test get_win(0, :hann) == Float64[]

        # Periodic window
        w_per = get_win(64, :hann; symmetric=false)
        @test length(w_per) == 64
    end

    @testset "smooth" begin
        # Smooth a step function
        m = zeros(32, 32)
        m[10:22, 10:22] .= 1.0
        smoothed = smooth(m)
        @test size(smoothed) == (32, 32)
        # Smoothed version should have less extreme values
        @test maximum(smoothed) <= maximum(m) + 0.1

        # With restore_max
        smoothed_rm = smooth(m; restore_max=true)
        @test maximum(smoothed_rm) ≈ maximum(m) atol=0.1
    end

    @testset "db2neper / neper2db" begin
        # Round-trip
        db_val = 10.0
        @test neper2db(db2neper(db_val)) ≈ db_val
        @test db2neper(neper2db(1.0)) ≈ 1.0

        # Known conversion: 1 Np ≈ 8.686 dB
        @test neper2db(1.0) ≈ 8.685889638 atol=1e-6

        # Vectorized
        db_vec = [0.0, 10.0, 20.0]
        np_vec = db2neper(db_vec)
        @test length(np_vec) == 3
        @test np_vec[1] ≈ 0.0
    end
end
