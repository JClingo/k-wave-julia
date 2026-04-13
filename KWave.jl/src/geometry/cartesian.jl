# ============================================================================
# KWave.jl — Cartesian geometry functions (off-grid point distributions)
# ============================================================================

"""
    make_cart_circle(radius, num_points, center; arc_angle=2π)

Create a set of Cartesian points evenly distributed over a circle.

# Arguments
- `radius`: Circle radius [m]
- `num_points`: Number of points on the circle
- `center`: `(cx, cy)` center position [m]
- `arc_angle`: Arc angle in radians (default: 2π for full circle)

# Returns
Matrix of size (2, num_points) with [x; y] coordinates.
"""
function make_cart_circle(radius::Real, num_points::Int,
                          center::NTuple{2,<:Real}=(0.0, 0.0);
                          arc_angle::Real=2π)
    cx, cy = Float64.(center)
    coords = zeros(Float64, 2, num_points)
    for i in 1:num_points
        θ = (i - 1) * arc_angle / num_points
        coords[1, i] = cx + radius * cos(θ)
        coords[2, i] = cy + radius * sin(θ)
    end
    return coords
end

"""
    make_cart_sphere(radius, num_points, center; ax_offset=0.0)

Create a set of Cartesian points approximately evenly distributed over a sphere.

Uses the Fibonacci sphere algorithm for approximately uniform distribution.

# Arguments
- `radius`: Sphere radius [m]
- `num_points`: Number of points on the sphere
- `center`: `(cx, cy, cz)` center position [m]
- `ax_offset`: Angular offset around axis [rad]

# Returns
Matrix of size (3, num_points) with [x; y; z] coordinates.
"""
function make_cart_sphere(radius::Real, num_points::Int,
                          center::NTuple{3,<:Real}=(0.0, 0.0, 0.0);
                          ax_offset::Float64=0.0)
    cx, cy, cz = Float64.(center)
    coords = zeros(Float64, 3, num_points)
    golden_ratio = (1 + sqrt(5)) / 2

    for i in 1:num_points
        θ = acos(1 - 2 * (i - 0.5) / num_points)
        φ = 2π * (i - 1) / golden_ratio + ax_offset
        coords[1, i] = cx + radius * sin(θ) * cos(φ)
        coords[2, i] = cy + radius * sin(θ) * sin(φ)
        coords[3, i] = cz + radius * cos(θ)
    end
    return coords
end

"""
    make_cart_arc(center, radius, diameter, focus, num_points)

Create a set of Cartesian points evenly distributed over a 2D arc.

# Arguments
- `center`: `(cx, cy)` center of arc curvature [m]
- `radius`: Radius of curvature [m]
- `diameter`: Aperture diameter [m]
- `focus`: `(fx, fy)` focus position [m]
- `num_points`: Number of points on the arc

# Returns
Matrix of size (2, num_points) with [x; y] coordinates.
"""
function make_cart_arc(center::NTuple{2,<:Real}, radius::Real, diameter::Real,
                       focus::NTuple{2,<:Real}, num_points::Int)
    cx, cy = Float64.(center)
    fx, fy = Float64.(focus)

    # Direction from center to focus
    dir_x = fx - cx
    dir_y = fy - cy
    center_angle = atan(dir_y, dir_x)

    # Half-angle subtended by aperture
    half_angle = asin(clamp(diameter / (2 * radius), -1.0, 1.0))

    coords = zeros(Float64, 2, num_points)
    for i in 1:num_points
        θ = center_angle - half_angle + (i - 1) * 2 * half_angle / (num_points - 1)
        coords[1, i] = cx + radius * cos(θ)
        coords[2, i] = cy + radius * sin(θ)
    end
    return coords
end

"""
    make_cart_bowl(center, radius, diameter, focus, num_points)

Create a set of Cartesian points evenly distributed over a 3D bowl surface.

# Arguments
- `center`: `(cx, cy, cz)` center of bowl curvature [m]
- `radius`: Radius of curvature [m]
- `diameter`: Aperture diameter [m]
- `focus`: `(fx, fy, fz)` focus position [m]
- `num_points`: Number of points on the bowl

# Returns
Matrix of size (3, num_points) with [x; y; z] coordinates.
"""
function make_cart_bowl(center::NTuple{3,<:Real}, radius::Real, diameter::Real,
                        focus::NTuple{3,<:Real}, num_points::Int)
    cx, cy, cz = Float64.(center)
    fx, fy, fz = Float64.(focus)

    # Direction from center to focus
    dir = Float64.([fx - cx, fy - cy, fz - cz])
    norm_dir = norm(dir)
    if norm_dir ≈ 0
        return zeros(Float64, 3, 0)
    end
    dir ./= norm_dir

    # Half-angle subtended by aperture
    half_angle = asin(clamp(diameter / (2 * radius), -1.0, 1.0))

    # Create points in local coordinate system (z-axis = dir), then rotate
    # Build orthonormal basis
    e3 = dir
    e1 = _perpendicular_vec(e3)
    e2 = cross(e3, e1)

    coords = zeros(Float64, 3, num_points)
    golden_ratio = (1 + sqrt(5)) / 2

    # Distribute points over spherical cap using Fibonacci-like spiral
    cos_half = cos(half_angle)
    idx = 0
    for i in 1:num_points
        # Map to cap: cos(theta) goes from cos(half_angle) to 1
        cos_θ = cos_half + (1 - cos_half) * (i - 0.5) / num_points
        sin_θ = sqrt(max(0.0, 1 - cos_θ^2))
        φ = 2π * (i - 1) / golden_ratio

        # Local coordinates
        local_x = sin_θ * cos(φ)
        local_y = sin_θ * sin(φ)
        local_z = cos_θ

        # Transform to global coordinates
        idx += 1
        coords[1, idx] = cx + radius * (e1[1] * local_x + e2[1] * local_y + e3[1] * local_z)
        coords[2, idx] = cy + radius * (e1[2] * local_x + e2[2] * local_y + e3[2] * local_z)
        coords[3, idx] = cz + radius * (e1[3] * local_x + e2[3] * local_y + e3[3] * local_z)
    end

    return coords
end

"""
    make_cart_disc(center, radius, focus, num_points; axial_offset=0.0)

Create a set of Cartesian points evenly distributed over a flat disc.

# Arguments
- `center`: `(cx, cy, cz)` center of disc [m]
- `radius`: Disc radius [m]
- `focus`: `(fx, fy, fz)` direction vector (disc normal points toward focus)
- `num_points`: Number of points on the disc
- `axial_offset`: Offset along the normal direction [m]

# Returns
Matrix of size (3, num_points) with [x; y; z] coordinates.
"""
function make_cart_disc(center::NTuple{3,<:Real}, radius::Real,
                        focus::NTuple{3,<:Real}, num_points::Int;
                        axial_offset::Float64=0.0)
    cx, cy, cz = Float64.(center)
    fx, fy, fz = Float64.(focus)

    # Normal direction
    dir = Float64.([fx - cx, fy - cy, fz - cz])
    norm_dir = norm(dir)
    if norm_dir ≈ 0
        dir = Float64.([0.0, 0.0, 1.0])
    else
        dir ./= norm_dir
    end

    # Build orthonormal basis
    e3 = dir
    e1 = _perpendicular_vec(e3)
    e2 = cross(e3, e1)

    # Offset center along normal
    cx_off = cx + axial_offset * dir[1]
    cy_off = cy + axial_offset * dir[2]
    cz_off = cz + axial_offset * dir[3]

    # Sunflower pattern for uniform disc distribution
    coords = zeros(Float64, 3, num_points)
    golden_angle = π * (3 - sqrt(5))

    for i in 1:num_points
        r = radius * sqrt((i - 0.5) / num_points)
        θ = i * golden_angle

        local_x = r * cos(θ)
        local_y = r * sin(θ)

        coords[1, i] = cx_off + e1[1] * local_x + e2[1] * local_y
        coords[2, i] = cy_off + e1[2] * local_x + e2[2] * local_y
        coords[3, i] = cz_off + e1[3] * local_x + e2[3] * local_y
    end

    return coords
end

"""
    make_cart_rect(center, Lx, Ly, focus, num_points)

Create a set of Cartesian points evenly distributed over a rectangle.

# Arguments
- `center`: `(cx, cy, cz)` center of rectangle [m]
- `Lx`: Rectangle width [m]
- `Ly`: Rectangle height [m]
- `focus`: `(fx, fy, fz)` direction vector (normal direction)
- `num_points`: Approximate number of points

# Returns
Matrix of size (3, N) with [x; y; z] coordinates.
"""
function make_cart_rect(center::NTuple{3,<:Real}, Lx::Real, Ly::Real,
                        focus::NTuple{3,<:Real}, num_points::Int)
    cx, cy, cz = Float64.(center)
    fx, fy, fz = Float64.(focus)

    # Normal direction
    dir = Float64.([fx - cx, fy - cy, fz - cz])
    norm_dir = norm(dir)
    if norm_dir ≈ 0
        dir = Float64.([0.0, 0.0, 1.0])
    else
        dir ./= norm_dir
    end

    # Build orthonormal basis
    e3 = dir
    e1 = _perpendicular_vec(e3)
    e2 = cross(e3, e1)

    # Compute grid dimensions for approximately num_points
    aspect = Lx / Ly
    ny_approx = sqrt(num_points / aspect)
    ny = max(1, round(Int, ny_approx))
    nx = max(1, round(Int, num_points / ny))
    actual_points = nx * ny

    coords = zeros(Float64, 3, actual_points)
    idx = 0
    for iy in 1:ny, ix in 1:nx
        idx += 1
        local_x = Lx * ((ix - 0.5) / nx - 0.5)
        local_y = Ly * ((iy - 0.5) / ny - 0.5)
        coords[1, idx] = cx + e1[1] * local_x + e2[1] * local_y
        coords[2, idx] = cy + e1[2] * local_x + e2[2] * local_y
        coords[3, idx] = cz + e1[3] * local_x + e2[3] * local_y
    end

    return coords
end

# ============================================================================
# Helper: find a vector perpendicular to a given unit vector
# ============================================================================

function _perpendicular_vec(v::Vector{Float64})
    # Pick the axis least aligned with v
    if abs(v[1]) <= abs(v[2]) && abs(v[1]) <= abs(v[3])
        candidate = Float64[1.0, 0.0, 0.0]
    elseif abs(v[2]) <= abs(v[3])
        candidate = Float64[0.0, 1.0, 0.0]
    else
        candidate = Float64[0.0, 0.0, 1.0]
    end
    perp = cross(v, candidate)
    perp ./= norm(perp)
    return perp
end
