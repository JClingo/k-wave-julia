# ============================================================================
# KWave.jl — KWaveTransducer: Linear array transducer model
# ============================================================================

"""
    KWaveTransducer

Linear array transducer model for 3D simulations.

Models a phased/linear array with configurable number of elements,
element geometry, focus, steering, and apodization.
"""
Base.@kwdef struct KWaveTransducer
    number_elements::Int
    element_width::Int            # grid points
    element_length::Int           # grid points
    element_spacing::Int = 0      # grid points (kerf)
    position::NTuple{3,Int} = (1, 1, 1)
    radius::Float64 = Inf        # focus radius (Inf = flat)
    focus_distance::Float64 = Inf
    steering_angle::Float64 = 0.0  # [rad]
    transmit_apodization::Symbol = :rectangular
    receive_apodization::Symbol = :rectangular
    active_elements::Union{Nothing, Vector{Int}} = nothing
    input_signal::Union{Nothing, AbstractVector{Float64}} = nothing
end

"""
    get_transducer_binary_mask(transducer, kgrid)

Generate a binary mask for the transducer mapped onto a 3D grid.

The transducer elements are laid out along the y-axis at the position specified
by `transducer.position`. Each element spans `element_width` points in y and
`element_length` points in z.

# Arguments
- `transducer`: KWaveTransducer
- `kgrid`: KWaveGrid3D

# Returns
`BitArray{3}` of grid size with `true` at transducer element locations.
"""
function get_transducer_binary_mask(transducer::KWaveTransducer, kgrid::KWaveGrid3D)
    mask = falses(kgrid.Nx, kgrid.Ny, kgrid.Nz)
    px, py, pz = transducer.position

    active = transducer.active_elements
    if active === nothing
        active = collect(1:transducer.number_elements)
    end

    for elem_idx in active
        # Y-offset for this element
        y_start = py + (elem_idx - 1) * (transducer.element_width + transducer.element_spacing)
        y_end = y_start + transducer.element_width - 1

        z_start = pz
        z_end = pz + transducer.element_length - 1

        for k in max(1, z_start):min(kgrid.Nz, z_end)
            for j in max(1, y_start):min(kgrid.Ny, y_end)
                if 1 <= px <= kgrid.Nx
                    mask[px, j, k] = true
                end
            end
        end
    end

    return mask
end

"""
    get_transducer_source(transducer, kgrid, medium)

Generate source signals for the transducer with appropriate delays
for focusing and steering.

# Arguments
- `transducer`: KWaveTransducer
- `kgrid`: KWaveGrid3D
- `medium`: KWaveMedium

# Returns
`(mask, source_signal)` — binary mask and time-delayed source signals.
"""
function get_transducer_source(transducer::KWaveTransducer, kgrid::KWaveGrid3D,
                               medium::KWaveMedium)
    if transducer.input_signal === nothing
        error("Transducer input_signal must be set before generating source")
    end

    mask = get_transducer_binary_mask(transducer, kgrid)
    mask_indices = findall(mask)
    n_points = length(mask_indices)
    Nt = kgrid.Nt[]
    dt = kgrid.dt[]

    c0 = medium.sound_speed isa Real ? medium.sound_speed : mean(medium.sound_speed)
    px, py, pz = transducer.position

    # Compute per-element delays
    active = transducer.active_elements
    if active === nothing
        active = collect(1:transducer.number_elements)
    end

    # Element center positions (in grid coordinates)
    element_centers_y = Float64[]
    for elem_idx in active
        y_center = py + (elem_idx - 1) * (transducer.element_width + transducer.element_spacing) +
                   transducer.element_width / 2
        push!(element_centers_y, y_center)
    end

    # Compute transmit apodization
    tx_apod = _compute_apodization(length(active), transducer.transmit_apodization)

    # Compute focusing delays
    delays = zeros(Float64, length(active))
    if isfinite(transducer.focus_distance)
        # Focus point
        array_center_y = mean(element_centers_y)
        for (i, y_c) in enumerate(element_centers_y)
            dy = (y_c - array_center_y) * kgrid.dy
            dist = sqrt(transducer.focus_distance^2 + dy^2)
            delays[i] = (dist - transducer.focus_distance) / c0
        end
    end

    # Add steering delays
    if transducer.steering_angle != 0.0
        array_center_y = mean(element_centers_y)
        for (i, y_c) in enumerate(element_centers_y)
            dy = (y_c - array_center_y) * kgrid.dy
            delays[i] += dy * sin(transducer.steering_angle) / c0
        end
    end

    # Normalize delays (make minimum zero)
    delays .-= minimum(delays)

    # Generate time-delayed source signal for each grid point
    source_signal = zeros(Float64, n_points, Nt)
    sig = transducer.input_signal
    sig_len = length(sig)

    point_idx = 0
    for (elem_i, elem_idx) in enumerate(active)
        y_start = py + (elem_idx - 1) * (transducer.element_width + transducer.element_spacing)
        y_end = y_start + transducer.element_width - 1
        z_start = pz
        z_end = pz + transducer.element_length - 1

        delay_samples = round(Int, delays[elem_i] / dt)

        for k in max(1, z_start):min(kgrid.Nz, z_end)
            for j in max(1, y_start):min(kgrid.Ny, y_end)
                if 1 <= px <= kgrid.Nx
                    point_idx += 1
                    for t in 1:min(sig_len, Nt - delay_samples)
                        source_signal[point_idx, t + delay_samples] = sig[t] * tx_apod[elem_i]
                    end
                end
            end
        end
    end

    return mask, source_signal
end

"""
    combine_transducer_sensor_data(transducer, kgrid, sensor_data)

Combine sensor data from grid points back into per-element signals,
applying receive apodization.

# Arguments
- `transducer`: KWaveTransducer
- `kgrid`: KWaveGrid3D
- `sensor_data`: Matrix of sensor data (n_sensor_points × Nt)

# Returns
Matrix of per-element sensor data (num_active_elements × Nt).
"""
function combine_transducer_sensor_data(transducer::KWaveTransducer, kgrid::KWaveGrid3D,
                                        sensor_data::AbstractMatrix)
    active = transducer.active_elements
    if active === nothing
        active = collect(1:transducer.number_elements)
    end

    n_active = length(active)
    Nt = size(sensor_data, 2)
    combined = zeros(Float64, n_active, Nt)
    rx_apod = _compute_apodization(n_active, transducer.receive_apodization)

    px, py, pz = transducer.position
    point_idx = 0

    for (elem_i, elem_idx) in enumerate(active)
        y_start = py + (elem_idx - 1) * (transducer.element_width + transducer.element_spacing)
        y_end = y_start + transducer.element_width - 1
        z_start = pz
        z_end = pz + transducer.element_length - 1

        n_pts = 0
        for k in max(1, z_start):min(kgrid.Nz, z_end)
            for j in max(1, y_start):min(kgrid.Ny, y_end)
                if 1 <= px <= kgrid.Nx
                    point_idx += 1
                    if point_idx <= size(sensor_data, 1)
                        combined[elem_i, :] .+= sensor_data[point_idx, :]
                        n_pts += 1
                    end
                end
            end
        end
        if n_pts > 0
            combined[elem_i, :] .*= rx_apod[elem_i] / n_pts
        end
    end

    return combined
end

# ============================================================================
# Apodization
# ============================================================================

function _compute_apodization(n_elements::Int, type::Symbol)
    if type == :rectangular
        return ones(Float64, n_elements)
    elseif type == :hann || type == :hanning
        return get_win(n_elements, :hann)
    elseif type == :hamming
        return get_win(n_elements, :hamming)
    elseif type == :blackman
        return get_win(n_elements, :blackman)
    elseif type == :triangular || type == :bartlett
        return get_win(n_elements, :bartlett)
    else
        error("Unknown apodization type: $type")
    end
end
