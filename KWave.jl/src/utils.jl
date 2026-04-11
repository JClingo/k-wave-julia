# ============================================================================
# KWave.jl — Grid and matrix utilities
# ============================================================================

"""
    cart2grid(kgrid, cart_data)

Map Cartesian points onto the nearest grid points.

# Arguments
- `kgrid`: AbstractKWaveGrid
- `cart_data`: Matrix of Cartesian coordinates. For 2D: (2, n_points) with rows [x; y].
  For 3D: (3, n_points) with rows [x; y; z].

# Returns
- `grid_mask`: BitArray with `true` at grid points nearest to the Cartesian points
- `order_index`: Indices mapping Cartesian points to their grid positions
- `reorder_index`: Inverse mapping to recover original point order
"""
function cart2grid(kgrid::KWaveGrid2D, cart_data::AbstractMatrix)
    n_points = size(cart_data, 2)

    # Map each Cartesian point to nearest grid index
    ix = [argmin(abs.(kgrid.x_vec .- cart_data[1, j])) for j in 1:n_points]
    iy = [argmin(abs.(kgrid.y_vec .- cart_data[2, j])) for j in 1:n_points]

    # Create binary mask
    grid_mask = falses(kgrid.Nx, kgrid.Ny)
    for j in 1:n_points
        grid_mask[ix[j], iy[j]] = true
    end

    # Compute linear indices for ordering
    lin_idx = [LinearIndices((kgrid.Nx, kgrid.Ny))[ix[j], iy[j]] for j in 1:n_points]
    order_index = sortperm(lin_idx)
    reorder_index = invperm(order_index)

    return grid_mask, order_index, reorder_index
end

function cart2grid(kgrid::KWaveGrid3D, cart_data::AbstractMatrix)
    n_points = size(cart_data, 2)

    ix = [argmin(abs.(kgrid.x_vec .- cart_data[1, j])) for j in 1:n_points]
    iy = [argmin(abs.(kgrid.y_vec .- cart_data[2, j])) for j in 1:n_points]
    iz = [argmin(abs.(kgrid.z_vec .- cart_data[3, j])) for j in 1:n_points]

    grid_mask = falses(kgrid.Nx, kgrid.Ny, kgrid.Nz)
    for j in 1:n_points
        grid_mask[ix[j], iy[j], iz[j]] = true
    end

    lin_idx = [LinearIndices((kgrid.Nx, kgrid.Ny, kgrid.Nz))[ix[j], iy[j], iz[j]] for j in 1:n_points]
    order_index = sortperm(lin_idx)
    reorder_index = invperm(order_index)

    return grid_mask, order_index, reorder_index
end

"""
    grid2cart(kgrid, mask)

Extract Cartesian coordinates from grid points where `mask` is `true`.

# Arguments
- `kgrid`: AbstractKWaveGrid
- `mask`: BitArray of grid size

# Returns
Matrix of Cartesian coordinates (ndims × n_points).
"""
function grid2cart(kgrid::KWaveGrid2D, mask::AbstractArray{Bool})
    indices = findall(mask)
    n = length(indices)
    coords = zeros(Float64, 2, n)
    for (j, idx) in enumerate(indices)
        coords[1, j] = kgrid.x_vec[idx[1]]
        coords[2, j] = kgrid.y_vec[idx[2]]
    end
    return coords
end

function grid2cart(kgrid::KWaveGrid3D, mask::AbstractArray{Bool})
    indices = findall(mask)
    n = length(indices)
    coords = zeros(Float64, 3, n)
    for (j, idx) in enumerate(indices)
        coords[1, j] = kgrid.x_vec[idx[1]]
        coords[2, j] = kgrid.y_vec[idx[2]]
        coords[3, j] = kgrid.z_vec[idx[3]]
    end
    return coords
end

"""
    expand_matrix(matrix, expand_size; edge_val=0)

Expand a matrix by adding extra elements at each edge.

# Arguments
- `matrix`: Input array (1D, 2D, or 3D)
- `expand_size`: Number of elements to add on each side. Can be:
  - `Int`: Same expansion on all sides of all dimensions
  - `Tuple{Int,...}`: Per-dimension expansion (same on both sides)
  - `Tuple{Tuple{Int,Int},...}`: Per-dimension, per-side expansion
- `edge_val`: Value to fill expanded regions (default: 0)

# Returns
Expanded array.
"""
function expand_matrix(matrix::AbstractArray{T}, expand_size; edge_val::T=zero(T)) where T
    nd = ndims(matrix)
    orig_size = size(matrix)

    # Normalize expand_size to per-dimension, per-side format
    if expand_size isa Integer
        pads = ntuple(_ -> (Int(expand_size), Int(expand_size)), nd)
    elseif expand_size isa Tuple
        if first(expand_size) isa Tuple
            pads = expand_size
        else
            pads = ntuple(i -> (Int(expand_size[i]), Int(expand_size[i])), nd)
        end
    else
        error("expand_size must be Int, Tuple of Int, or Tuple of Tuple{Int,Int}")
    end

    # Compute new size
    new_size = ntuple(i -> orig_size[i] + pads[i][1] + pads[i][2], nd)

    # Allocate and fill
    result = fill(edge_val, new_size...)

    # Copy original data into center
    ranges = ntuple(i -> (pads[i][1]+1):(pads[i][1]+orig_size[i]), nd)
    result[ranges...] = matrix

    return result
end

"""
    resize_array(matrix, new_size)

Resize an array using linear interpolation.

# Arguments
- `matrix`: Input array (1D, 2D, or 3D)
- `new_size`: Tuple of target dimensions

# Returns
Resized array.
"""
function resize_array(matrix::AbstractVector, new_size::Tuple{Int})
    N = length(matrix)
    Nnew = new_size[1]
    if N == Nnew
        return copy(matrix)
    end
    # Create interpolation
    itp = interpolate(matrix, BSpline(Linear()))
    sitp = scale(itp, range(1, N, length=N))
    # Sample at new positions
    new_coords = range(1, N, length=Nnew)
    return [sitp(x) for x in new_coords]
end

function resize_array(matrix::AbstractMatrix, new_size::Tuple{Int, Int})
    Nx, Ny = size(matrix)
    Nx_new, Ny_new = new_size
    if (Nx, Ny) == new_size
        return copy(matrix)
    end
    itp = interpolate(matrix, BSpline(Linear()))
    sitp = scale(itp, range(1, Nx, length=Nx), range(1, Ny, length=Ny))
    xs = range(1, Nx, length=Nx_new)
    ys = range(1, Ny, length=Ny_new)
    return [sitp(x, y) for x in xs, y in ys]
end

function resize_array(matrix::AbstractArray{T, 3}, new_size::Tuple{Int, Int, Int}) where T
    Nx, Ny, Nz = size(matrix)
    Nx_new, Ny_new, Nz_new = new_size
    if (Nx, Ny, Nz) == new_size
        return copy(matrix)
    end
    itp = interpolate(matrix, BSpline(Linear()))
    sitp = scale(itp, range(1, Nx, length=Nx), range(1, Ny, length=Ny), range(1, Nz, length=Nz))
    xs = range(1, Nx, length=Nx_new)
    ys = range(1, Ny, length=Ny_new)
    zs = range(1, Nz, length=Nz_new)
    return [sitp(x, y, z) for x in xs, y in ys, z in zs]
end
