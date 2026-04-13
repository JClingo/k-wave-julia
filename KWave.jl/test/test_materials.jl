@testset "Material Properties (Phase 3)" begin
    @testset "water_sound_speed" begin
        # At 20°C, water sound speed ≈ 1482 m/s
        c = water_sound_speed(20.0)
        @test isapprox(c, 1482.0; atol=5.0)

        # At 37°C (body temperature), ≈ 1524 m/s
        c37 = water_sound_speed(37.0)
        @test isapprox(c37, 1524.0; atol=5.0)

        # Sound speed increases with temperature (up to ~74°C)
        @test water_sound_speed(30.0) > water_sound_speed(20.0)
    end

    @testset "water_density" begin
        # At 20°C, water density ≈ 998 kg/m³
        rho = water_density(20.0)
        @test isapprox(rho, 998.2; atol=1.0)

        # Maximum density near 4°C
        @test water_density(4.0) > water_density(20.0)
    end

    @testset "water_absorption" begin
        # Absorption should be positive
        alpha = water_absorption(20.0, 1e6)
        @test alpha > 0

        # Absorption increases with frequency
        @test water_absorption(20.0, 5e6) > water_absorption(20.0, 1e6)
    end

    @testset "water_nonlinearity" begin
        # B/A for water at 20°C ≈ 5
        ba = water_nonlinearity(20.0)
        @test isapprox(ba, 5.6; atol=0.5)
    end

    @testset "fit_power_law_params" begin
        # Generate test data with known power law
        alpha0 = 0.5
        y = 1.5
        f = [1e6, 2e6, 5e6, 10e6]
        alpha = alpha0 .* (f ./ 1e6).^y
        alpha_coeff, alpha_power = fit_power_law_params(f, alpha, 1500.0)
        @test isapprox(alpha_coeff, alpha0; rtol=0.01)
        @test isapprox(alpha_power, y; rtol=0.01)
    end
end
