@testset "Beamforming and Scan Conversion" begin
    @testset "beamform_delay_and_sum basic" begin
        # Create simple synthetic data — point source at center
        n_sensors = 32
        n_samples = 100
        c0 = 1500.0
        dt = 1e-7

        # Sensor positions — linear array
        sensor_pos = range(-5e-3, 5e-3, length=n_sensors)

        # Point source at (0, 5mm) — create time-delayed signals
        xp, zp = 0.0, 5e-3
        sensor_data = zeros(n_sensors, n_samples)
        for (i, sx) in enumerate(sensor_pos)
            dist = sqrt((sx - xp)^2 + zp^2)
            sample = round(Int, dist / (c0 * dt)) + 1
            if 1 <= sample <= n_samples
                sensor_data[i, sample] = 1.0
            end
        end

        # Reconstruct
        grid_x = range(-5e-3, 5e-3, length=20)
        grid_z = range(1e-3, 10e-3, length=20)

        image = beamform_delay_and_sum(
            sensor_data, collect(sensor_pos), c0, dt,
            collect(grid_x), collect(grid_z),
        )

        @test size(image) == (20, 20)
        # Maximum should be near the true source location
        @test maximum(image) > 0
    end

    @testset "beamform_delay_and_sum apodization" begin
        n_sensors = 16
        n_samples = 50
        sensor_data = randn(n_sensors, n_samples)
        sensor_pos = range(-2e-3, 2e-3, length=n_sensors)
        grid_x = range(-2e-3, 2e-3, length=10)
        grid_z = range(0.5e-3, 5e-3, length=10)

        for apod in [:rectangular, :hamming, :hann, :blackman]
            image = beamform_delay_and_sum(
                sensor_data, collect(sensor_pos), 1500.0, 1e-7,
                collect(grid_x), collect(grid_z); apodization=apod,
            )
            @test size(image) == (10, 10)
            @test all(isfinite, image)
        end
    end

    @testset "beamform_delay_and_sum f_number" begin
        n_sensors = 16
        n_samples = 50
        sensor_data = randn(n_sensors, n_samples)
        sensor_pos = range(-2e-3, 2e-3, length=n_sensors)
        grid_x = range(-2e-3, 2e-3, length=10)
        grid_z = range(0.5e-3, 5e-3, length=10)

        image = beamform_delay_and_sum(
            sensor_data, collect(sensor_pos), 1500.0, 1e-7,
            collect(grid_x), collect(grid_z); f_number=2.0,
        )
        @test size(image) == (10, 10)
    end

    @testset "scan_conversion basic" begin
        n_lines = 16
        n_samples = 100

        data = randn(n_lines, n_samples)
        angles = range(-π/6, π/6, length=n_lines)
        positions = zeros(n_lines)

        image, x_vec, z_vec = scan_conversion(
            data, collect(angles), positions, 1500.0, 1e-7;
            image_size=(32, 32),
        )

        @test size(image) == (32, 32)
        @test length(x_vec) == 32
        @test length(z_vec) == 32
    end

    @testset "scan_conversion custom range" begin
        data = randn(8, 50)
        angles = range(-π/8, π/8, length=8)
        positions = zeros(8)

        image, x_vec, z_vec = scan_conversion(
            data, collect(angles), positions, 1500.0, 1e-7;
            image_size=(20, 20),
            x_range=(-3e-3, 3e-3),
            z_range=(0.0, 5e-3),
        )

        @test size(image) == (20, 20)
        @test x_vec[1] ≈ -3e-3
        @test x_vec[end] ≈ 3e-3
        @test z_vec[end] ≈ 5e-3
    end
end
