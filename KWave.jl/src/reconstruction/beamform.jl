# ============================================================================
# KWave.jl — Beamforming reconstruction and scan conversion
# ============================================================================

"""
    scan_conversion(data, steering_angles, element_positions, c0, dt;
                    image_size=(256, 256), x_range=nothing, z_range=nothing)

Convert scan-line data from polar/sector format to a Cartesian image.

Performs scan conversion of ultrasound B-mode data from a phased array
or sector scan into a rectangular image grid.

# Arguments
- `data`: Scan-line data matrix (num_elements × num_time_steps)
- `steering_angles`: Vector of steering angles [rad] for each scan line
- `element_positions`: Vector of element positions along the array [m]
- `c0`: Sound speed [m/s]
- `dt`: Time step [s]

# Keyword Arguments
- `image_size`: Output image size as (Nx, Nz) tuple (default: (256, 256))
- `x_range`: Lateral range [m] as (xmin, xmax) (default: auto from geometry)
- `z_range`: Depth range [m] as (zmin, zmax) (default: auto from data)

# Returns
`(image, x_vec, z_vec)` — scan-converted image and coordinate vectors.
"""
function scan_conversion(
    data::AbstractMatrix,
    steering_angles::AbstractVector,
    element_positions::AbstractVector,
    c0::Real,
    dt::Real;
    image_size::Tuple{Int, Int}=(256, 256),
    x_range::Union{Nothing, Tuple{Real, Real}}=nothing,
    z_range::Union{Nothing, Tuple{Real, Real}}=nothing,
)
    n_lines, n_samples = size(data)
    c0 = Float64(c0)
    dt = Float64(dt)

    # Depth vector from time samples
    depth = collect(0:n_samples-1) .* (c0 * dt / 2)  # divide by 2 for pulse-echo

    # Auto-compute ranges
    if z_range === nothing
        z_range = (0.0, depth[end])
    end

    max_angle = maximum(abs, steering_angles)
    max_depth = z_range[2]
    if x_range === nothing
        lateral_extent = max_depth * sin(max_angle) + maximum(abs, element_positions)
        x_range = (-lateral_extent, lateral_extent)
    end

    Nx_out, Nz_out = image_size
    x_vec = range(x_range[1], x_range[2], length=Nx_out)
    z_vec = range(z_range[1], z_range[2], length=Nz_out)

    image = zeros(Float64, Nx_out, Nz_out)

    for (line_idx, angle) in enumerate(steering_angles)
        for (iz, z) in enumerate(z_vec)
            for (ix, x) in enumerate(x_vec)
                # Convert (x, z) to range and angle relative to this scan line
                r = sqrt(x^2 + z^2)
                theta = atan(x, z)

                # Angular distance from this scan line
                d_angle = abs(theta - angle)

                # Only include points within ±half beam width
                if d_angle < π / (2 * n_lines) || n_lines == 1
                    # Interpolate along the scan line depth
                    sample_idx = 2 * r / (c0 * dt) + 1  # +1 for Julia indexing

                    if 1 <= sample_idx <= n_samples
                        # Linear interpolation
                        idx_low = floor(Int, sample_idx)
                        idx_high = min(idx_low + 1, n_samples)
                        frac = sample_idx - idx_low

                        val = (1 - frac) * data[line_idx, idx_low] + frac * data[line_idx, idx_high]
                        image[ix, iz] += val
                    end
                end
            end
        end
    end

    return image, collect(x_vec), collect(z_vec)
end

"""
    beamform_delay_and_sum(sensor_data, sensor_positions, c0, dt, grid_x, grid_z;
                           apodization=:rectangular, f_number=nothing)

Delay-and-sum beamforming reconstruction.

Reconstructs an image from sensor data by computing time delays
from each pixel to each sensor element and summing the coherently
delayed signals.

# Arguments
- `sensor_data`: Recorded pressure data (num_sensors × num_time_steps)
- `sensor_positions`: Sensor coordinates — matrix (2 × num_sensors) for 2D [x; z],
  or vector of x-positions for line array at z=0
- `c0`: Sound speed [m/s]
- `dt`: Time step [s]
- `grid_x`: Reconstruction x-coordinates [m] (vector)
- `grid_z`: Reconstruction z-coordinates [m] (vector)

# Keyword Arguments
- `apodization`: Window function for element weighting (default: :rectangular)
- `f_number`: F-number for dynamic aperture (default: nothing = use all elements)

# Returns
Reconstructed image (length(grid_x) × length(grid_z)).
"""
function beamform_delay_and_sum(
    sensor_data::AbstractMatrix,
    sensor_positions::Union{AbstractVector, AbstractMatrix},
    c0::Real,
    dt::Real,
    grid_x::AbstractVector,
    grid_z::AbstractVector;
    apodization::Symbol=:rectangular,
    f_number::Union{Nothing, Real}=nothing,
)
    n_sensors, n_samples = size(sensor_data)
    c0 = Float64(c0)
    dt = Float64(dt)

    # Parse sensor positions
    if sensor_positions isa AbstractVector
        sx = Float64.(sensor_positions)
        sz = zeros(Float64, n_sensors)
    else
        sx = Float64.(sensor_positions[1, :])
        sz = Float64.(sensor_positions[2, :])
    end

    Nx = length(grid_x)
    Nz = length(grid_z)
    image = zeros(Float64, Nx, Nz)

    # Apodization weights
    weights = _compute_beamform_apodization(n_sensors, apodization)

    for iz in 1:Nz
        z = Float64(grid_z[iz])
        for ix in 1:Nx
            x = Float64(grid_x[ix])

            val = 0.0
            n_active = 0

            for s in 1:n_sensors
                # Dynamic aperture check
                if f_number !== nothing
                    half_aperture = z / (2 * f_number)
                    if abs(sx[s] - x) > half_aperture
                        continue
                    end
                end

                # Distance from pixel to sensor
                dist = sqrt((x - sx[s])^2 + (z - sz[s])^2)

                # Time delay → sample index
                delay = dist / c0
                sample_idx = delay / dt + 1  # +1 for Julia indexing

                if 1 <= sample_idx <= n_samples
                    idx_low = floor(Int, sample_idx)
                    idx_high = min(idx_low + 1, n_samples)
                    frac = sample_idx - idx_low

                    val += weights[s] * ((1 - frac) * sensor_data[s, idx_low] +
                                          frac * sensor_data[s, idx_high])
                    n_active += 1
                end
            end

            if n_active > 0
                image[ix, iz] = val / n_active
            end
        end
    end

    return image
end

"""
Compute apodization weights for beamforming.
"""
function _compute_beamform_apodization(n::Int, type::Symbol)
    if type == :rectangular
        return ones(Float64, n)
    elseif type == :hamming
        return [0.54 - 0.46 * cos(2π * i / (n - 1)) for i in 0:n-1]
    elseif type == :hann
        return [0.5 * (1 - cos(2π * i / (n - 1))) for i in 0:n-1]
    elseif type == :blackman
        return [0.42 - 0.5 * cos(2π * i / (n - 1)) + 0.08 * cos(4π * i / (n - 1)) for i in 0:n-1]
    else
        error("Unknown apodization type: $type. Use :rectangular, :hamming, :hann, or :blackman.")
    end
end
