@testset "Reconstruction (Phase 3)" begin
    @testset "kspace_line_recon" begin
        # Create synthetic sensor data
        Ny = 64
        Nt = 128
        dy = 1e-4
        dt = 1e-7

        # Simple test: Gaussian pulse arriving at different times
        sensor_data = zeros(Float64, Ny, Nt)
        for j in 1:Ny
            t0 = 40 + abs(j - Ny÷2) * 0.5
            t_idx = round(Int, t0)
            if 1 <= t_idx <= Nt
                sensor_data[j, t_idx] = 1.0
            end
        end

        recon = kspace_line_recon(sensor_data, dy, dt; c=1500.0)
        @test size(recon) == (Ny, Nt)
        @test maximum(recon) > 0
    end

    @testset "kspace_line_recon positivity" begin
        Ny, Nt = 32, 64
        sensor_data = randn(Ny, Nt)
        recon = kspace_line_recon(sensor_data, 1e-4, 1e-7; pos_cond=true)
        @test all(recon .>= 0)
    end

    @testset "kspace_plane_recon" begin
        Ny, Nz, Nt = 16, 16, 32
        sensor_data = randn(Ny, Nz, Nt) .* 0.01
        # Add a signal
        sensor_data[8, 8, 10] = 1.0

        recon = kspace_plane_recon(sensor_data, 1e-4, 1e-4, 1e-7; c=1500.0)
        @test size(recon) == (Ny, Nz, Nt)
    end
end
