# ============================================================================
# KWave.jl — KWaveArray: Off-grid transducer array modeling
# ============================================================================

# ============================================================================
# Array element types
# ============================================================================

abstract type ArrayElement end

"""Arc element in 2D."""
struct ArcElement <: ArrayElement
    position::NTuple{2,Float64}   # center of curvature
    radius::Float64               # radius of curvature [m]
    diameter::Float64             # aperture diameter [m]
    focus::NTuple{2,Float64}      # focus position [m]
end

"""Bowl element in 3D."""
struct BowlElement <: ArrayElement
    position::NTuple{3,Float64}
    radius::Float64
    diameter::Float64
    focus::NTuple{3,Float64}
end

"""Disc element in 3D (flat circular piston)."""
struct DiscElement <: ArrayElement
    position::NTuple{3,Float64}
    diameter::Float64
    focus::NTuple{3,Float64}
end

"""Rectangular element in 3D."""
struct RectElement <: ArrayElement
    position::NTuple{3,Float64}
    width::Float64
    height::Float64
    focus::NTuple{3,Float64}
end

"""Sphere segment element in 3D."""
struct SphereElement <: ArrayElement
    position::NTuple{3,Float64}
    radius::Float64
    diameter::Float64
end

# ============================================================================
# KWaveArray
# ============================================================================

"""
    KWaveArray

Off-grid transducer array model for multi-element transducers.

Elements are defined in physical coordinates (not grid points) and mapped
to the computational grid via `get_element_binary_mask`.
"""
struct KWaveArray
    elements::Vector{ArrayElement}
end

KWaveArray() = KWaveArray(ArrayElement[])

"""
    add_arc_element!(array, position, radius, diameter, focus)

Add a 2D arc element to the array.

# Arguments
- `position`: `(x, y)` center of curvature [m]
- `radius`: Radius of curvature [m]
- `diameter`: Aperture diameter [m]
- `focus`: `(fx, fy)` focus position [m]
"""
function add_arc_element!(array::KWaveArray, position::NTuple{2,<:Real},
                          radius::Real, diameter::Real, focus::NTuple{2,<:Real})
    push!(array.elements, ArcElement(Float64.(position), Float64(radius),
                                      Float64(diameter), Float64.(focus)))
    return array
end

"""
    add_bowl_element!(array, position, radius, diameter, focus)

Add a 3D bowl element to the array.
"""
function add_bowl_element!(array::KWaveArray, position::NTuple{3,<:Real},
                           radius::Real, diameter::Real, focus::NTuple{3,<:Real})
    push!(array.elements, BowlElement(Float64.(position), Float64(radius),
                                       Float64(diameter), Float64.(focus)))
    return array
end

"""
    add_disc_element!(array, position, diameter, focus)

Add a 3D flat disc element to the array.
"""
function add_disc_element!(array::KWaveArray, position::NTuple{3,<:Real},
                           diameter::Real, focus::NTuple{3,<:Real})
    push!(array.elements, DiscElement(Float64.(position), Float64(diameter),
                                       Float64.(focus)))
    return array
end

"""
    add_rect_element!(array, position, width, height, focus)

Add a 3D rectangular element to the array.
"""
function add_rect_element!(array::KWaveArray, position::NTuple{3,<:Real},
                           width::Real, height::Real, focus::NTuple{3,<:Real})
    push!(array.elements, RectElement(Float64.(position), Float64(width),
                                       Float64(height), Float64.(focus)))
    return array
end

Base.length(a::KWaveArray) = length(a.elements)

# ============================================================================
# Grid mapping
# ============================================================================

"""
    get_element_binary_mask(array, kgrid, element_index)

Get the binary mask for a single element mapped to the computational grid.

# Arguments
- `array`: KWaveArray
- `kgrid`: KWaveGrid
- `element_index`: Index of element to map

# Returns
Binary mask (BitArray) of grid size.
"""
function get_element_binary_mask(array::KWaveArray, kgrid::KWaveGrid2D, element_index::Int)
    elem = array.elements[element_index]
    return _element_to_mask(elem, kgrid)
end

function get_element_binary_mask(array::KWaveArray, kgrid::KWaveGrid3D, element_index::Int)
    elem = array.elements[element_index]
    return _element_to_mask(elem, kgrid)
end

function _element_to_mask(elem::ArcElement, kgrid::KWaveGrid2D)
    # Map physical position to grid indices
    cx_idx = argmin(abs.(kgrid.x_vec .- elem.position[1]))
    cy_idx = argmin(abs.(kgrid.y_vec .- elem.position[2]))
    fx_idx = argmin(abs.(kgrid.x_vec .- elem.focus[1]))
    fy_idx = argmin(abs.(kgrid.y_vec .- elem.focus[2]))

    r_grid = elem.radius / kgrid.dx
    d_grid = elem.diameter / kgrid.dx

    return make_arc(kgrid.Nx, kgrid.Ny, cx_idx, cy_idx, r_grid, d_grid, (fx_idx, fy_idx))
end

function _element_to_mask(elem::BowlElement, kgrid::KWaveGrid3D)
    cx = argmin(abs.(kgrid.x_vec .- elem.position[1]))
    cy = argmin(abs.(kgrid.y_vec .- elem.position[2]))
    cz = argmin(abs.(kgrid.z_vec .- elem.position[3]))
    fx = argmin(abs.(kgrid.x_vec .- elem.focus[1]))
    fy = argmin(abs.(kgrid.y_vec .- elem.focus[2]))
    fz = argmin(abs.(kgrid.z_vec .- elem.focus[3]))

    r_grid = elem.radius / kgrid.dx
    d_grid = elem.diameter / kgrid.dx

    return make_bowl(kgrid.Nx, kgrid.Ny, kgrid.Nz, (cx, cy, cz), r_grid, d_grid, (fx, fy, fz))
end

function _element_to_mask(elem::DiscElement, kgrid::KWaveGrid3D)
    cx = argmin(abs.(kgrid.x_vec .- elem.position[1]))
    cy = argmin(abs.(kgrid.y_vec .- elem.position[2]))
    cz = argmin(abs.(kgrid.z_vec .- elem.position[3]))
    r_grid = elem.diameter / (2 * kgrid.dx)
    return make_ball(kgrid.Nx, kgrid.Ny, kgrid.Nz, cx, cy, cz, r_grid)
end

function _element_to_mask(elem::RectElement, kgrid::KWaveGrid3D)
    cx = argmin(abs.(kgrid.x_vec .- elem.position[1]))
    cy = argmin(abs.(kgrid.y_vec .- elem.position[2]))
    cz = argmin(abs.(kgrid.z_vec .- elem.position[3]))

    hw = round(Int, elem.width / (2 * kgrid.dx))
    hh = round(Int, elem.height / (2 * kgrid.dy))

    mask = falses(kgrid.Nx, kgrid.Ny, kgrid.Nz)
    for j in max(1, cy - hh):min(kgrid.Ny, cy + hh)
        for i in max(1, cx - hw):min(kgrid.Nx, cx + hw)
            if 1 <= cz <= kgrid.Nz
                mask[i, j, cz] = true
            end
        end
    end
    return mask
end

"""
    get_array_binary_mask(array, kgrid)

Get the combined binary mask for all elements in the array.

# Returns
Binary mask (BitArray) of grid size with `true` at all element locations.
"""
function get_array_binary_mask(array::KWaveArray, kgrid::AbstractKWaveGrid)
    mask = falses(grid_size(kgrid)...)
    for i in 1:length(array)
        mask .|= get_element_binary_mask(array, kgrid, i)
    end
    return mask
end

"""
    get_distributed_source_signal(array, kgrid, source_signal)

Distribute a per-element source signal across grid points for each element.

For elements that map to multiple grid points, the source signal is divided
equally among those points.

# Arguments
- `array`: KWaveArray
- `kgrid`: KWaveGrid
- `source_signal`: Matrix of size (num_elements × Nt) source signals

# Returns
`(combined_mask, distributed_signal)` — combined binary mask and matrix of
grid-point signals (num_grid_points × Nt).
"""
function get_distributed_source_signal(array::KWaveArray, kgrid::AbstractKWaveGrid,
                                       source_signal::AbstractMatrix)
    n_elements = length(array)
    Nt = size(source_signal, 2)

    # Get combined mask and per-element masks
    combined_mask = falses(grid_size(kgrid)...)
    element_masks = Vector{BitArray}(undef, n_elements)

    for i in 1:n_elements
        element_masks[i] = get_element_binary_mask(array, kgrid, i)
        combined_mask .|= element_masks[i]
    end

    # Get indices of combined mask
    mask_indices = findall(combined_mask)
    n_points = length(mask_indices)
    distributed = zeros(Float64, n_points, Nt)

    # Map from global index → position in mask_indices
    idx_map = Dict{CartesianIndex, Int}()
    for (j, idx) in enumerate(mask_indices)
        idx_map[idx] = j
    end

    # Distribute each element's signal to its grid points
    for i in 1:n_elements
        elem_indices = findall(element_masks[i])
        n_elem_pts = length(elem_indices)
        if n_elem_pts == 0
            continue
        end
        for idx in elem_indices
            if haskey(idx_map, idx)
                j = idx_map[idx]
                distributed[j, :] .+= source_signal[i, :] ./ n_elem_pts
            end
        end
    end

    return combined_mask, distributed
end

"""
    combine_sensor_data(array, kgrid, sensor_data)

Combine grid-point sensor data back into per-element data by averaging
over each element's grid points.

# Arguments
- `array`: KWaveArray
- `kgrid`: KWaveGrid
- `sensor_data`: Matrix of recorded sensor data (num_sensor_points × Nt)

# Returns
Matrix of per-element sensor data (num_elements × Nt).
"""
function combine_sensor_data(array::KWaveArray, kgrid::AbstractKWaveGrid,
                             sensor_data::AbstractMatrix)
    n_elements = length(array)
    Nt = size(sensor_data, 2)
    combined = zeros(Float64, n_elements, Nt)

    # Get the sensor mask (assumed to be the combined array mask)
    combined_mask = get_array_binary_mask(array, kgrid)
    mask_indices = findall(combined_mask)

    idx_map = Dict{CartesianIndex, Int}()
    for (j, idx) in enumerate(mask_indices)
        idx_map[idx] = j
    end

    for i in 1:n_elements
        elem_mask = get_element_binary_mask(array, kgrid, i)
        elem_indices = findall(elem_mask)
        n_pts = 0
        for idx in elem_indices
            if haskey(idx_map, idx)
                j = idx_map[idx]
                if j <= size(sensor_data, 1)
                    combined[i, :] .+= sensor_data[j, :]
                    n_pts += 1
                end
            end
        end
        if n_pts > 0
            combined[i, :] ./= n_pts
        end
    end

    return combined
end
