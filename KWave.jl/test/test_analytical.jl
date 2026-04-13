@testset "Analytical Reference Solutions" begin
    @testset "focused_annulus_oneil" begin
        # O'Neil solution for unfocused annulus
        r_outer = 10e-3
        r_inner = 5e-3
        amp = 1.0
        phase = 0.0
        freq = 1e6
        c0 = 1500.0
        rho = 1000.0
        z_pos = range(1e-3, 50e-3, length=50)

        p = focused_annulus_oneil(r_outer, r_inner, amp, phase, freq, c0, rho, collect(z_pos))

        @test length(p) == 50
        @test all(isfinite, p)
        # Pressure at z=0 should be zero
        p_near = focused_annulus_oneil(r_outer, r_inner, amp, phase, freq, c0, rho, [0.0])
        @test abs(p_near[1]) ≈ 0.0 atol=1e-10
    end

    @testset "focused_annulus_oneil with focus" begin
        r_outer = 10e-3
        r_inner = 0.0
        amp = 1.0
        phase = 0.0
        freq = 1e6
        c0 = 1500.0
        rho = 1000.0
        focus = 30e-3

        z_pos = range(1e-3, 50e-3, length=100)
        p = focused_annulus_oneil(r_outer, r_inner, amp, phase, freq, c0, rho,
                                  collect(z_pos); focus_distance=focus)

        @test length(p) == 100
        @test all(isfinite, p)
    end

    @testset "focused_bowl_oneil" begin
        radius = 50e-3
        diameter = 40e-3
        amp = 1.0
        phase = 0.0
        freq = 1e6
        c0 = 1500.0
        rho = 1000.0

        z_pos = range(1e-3, 80e-3, length=100)
        p = focused_bowl_oneil(radius, diameter, amp, phase, freq, c0, rho, collect(z_pos))

        @test length(p) == 100
        @test all(isfinite, p)

        # Should have a focus region — pressure should vary along axis
        @test maximum(abs, p) > 0
    end

    @testset "mendousse" begin
        # Mendousse solution for nonlinear 1D propagation
        x = 5e-3  # 5 mm propagation distance
        freq = 1e6
        c0 = 1500.0
        rho = 1000.0
        BonA = 5.0
        alpha_0 = 0.5  # Np/m
        u0 = 1.0  # m/s source amplitude

        t = range(0, 2/freq, length=200)
        p = mendousse(x, collect(t), u0, freq, c0, rho, BonA, alpha_0; num_harmonics=20)

        @test length(p) == 200
        @test all(isfinite, p)
        # Should oscillate
        @test maximum(p) > 0
        @test minimum(p) < 0
    end

    @testset "mendousse vector x" begin
        x = [1e-3, 5e-3, 10e-3]
        freq = 1e6
        c0 = 1500.0
        rho = 1000.0
        BonA = 5.0
        alpha_0 = 0.5
        u0 = 0.5

        t = range(0, 2/freq, length=100)
        p = mendousse(x, collect(t), u0, freq, c0, rho, BonA, alpha_0)

        @test size(p) == (3, 100)
        # Amplitude should decrease with distance (due to absorption)
        @test maximum(abs, p[3, :]) < maximum(abs, p[1, :])
    end
end
