# ============================================================================
# KWave.jl — Geometry / shape creation functions
# ============================================================================

# ============================================================================
# 2D shapes
# ============================================================================

"""
    make_disc(Nx, Ny, cx, cy, radius; plot_disc=false)

Create a 2D binary mask of a filled disc (circle).

# Arguments
- `Nx`, `Ny`: Grid dimensions
- `cx`, `cy`: Center point (1-based grid indices)
- `radius`: Disc radius in grid points

# Returns
`BitMatrix` of size (Nx, Ny) with `true` inside the disc.
"""
function make_disc(Nx::Int, Ny::Int, cx::Int, cy::Int, radius::Real;
                   plot_disc::Bool=false)
    disc = falses(Nx, Ny)
    r2 = radius^2
    for j in 1:Ny
        for i in 1:Nx
            if (i - cx)^2 + (j - cy)^2 <= r2
                disc[i, j] = true
            end
        end
    end
    return disc
end

function make_disc(Nx::Int, Ny::Int, cx::Real, cy::Real, radius::Real; kwargs...)
    return make_disc(Nx, Ny, round(Int, cx), round(Int, cy), radius; kwargs...)
end

"""
    make_circle(Nx, Ny, cx, cy, radius; arc_angle=2π, plot_circle=false)

Create a 2D binary mask of a circle perimeter (single-pixel width).

Uses the midpoint circle algorithm for rasterization.

# Arguments
- `Nx`, `Ny`: Grid dimensions
- `cx`, `cy`: Center point (1-based grid indices)
- `radius`: Circle radius in grid points
- `arc_angle`: Arc angle in radians (default: 2π for full circle)

# Returns
`BitMatrix` of size (Nx, Ny) with `true` on the circle perimeter.
"""
function make_circle(Nx::Int, Ny::Int, cx::Int, cy::Int, radius::Int;
                     arc_angle::Float64=2π, plot_circle::Bool=false)
    circle = falses(Nx, Ny)

    if arc_angle >= 2π
        _midpoint_circle!(circle, Nx, Ny, cx, cy, radius)
    else
        n_points = max(100, ceil(Int, 4 * radius * arc_angle / (2π)))
        for i in 0:n_points-1
            θ = i * arc_angle / n_points
            px = round(Int, cx + radius * cos(θ))
            py = round(Int, cy + radius * sin(θ))
            if 1 <= px <= Nx && 1 <= py <= Ny
                circle[px, py] = true
            end
        end
    end

    return circle
end

function make_circle(Nx::Int, Ny::Int, cx::Real, cy::Real, radius::Real; kwargs...)
    return make_circle(Nx, Ny, round(Int, cx), round(Int, cy), round(Int, radius); kwargs...)
end

"""
Midpoint circle algorithm for rasterizing a full circle perimeter.
"""
function _midpoint_circle!(mask::BitMatrix, Nx::Int, Ny::Int, cx::Int, cy::Int, r::Int)
    x = r
    y = 0
    d = 1 - r

    while x >= y
        _set_point!(mask, Nx, Ny, cx + x, cy + y)
        _set_point!(mask, Nx, Ny, cx - x, cy + y)
        _set_point!(mask, Nx, Ny, cx + x, cy - y)
        _set_point!(mask, Nx, Ny, cx - x, cy - y)
        _set_point!(mask, Nx, Ny, cx + y, cy + x)
        _set_point!(mask, Nx, Ny, cx - y, cy + x)
        _set_point!(mask, Nx, Ny, cx + y, cy - x)
        _set_point!(mask, Nx, Ny, cx - y, cy - x)

        y += 1
        if d <= 0
            d += 2y + 1
        else
            x -= 1
            d += 2(y - x) + 1
        end
    end
end

function _set_point!(mask::BitMatrix, Nx::Int, Ny::Int, px::Int, py::Int)
    if 1 <= px <= Nx && 1 <= py <= Ny
        mask[px, py] = true
    end
end

"""
    make_arc(Nx, Ny, cx, cy, radius, diameter, focus; plot_arc=false)

Create a 2D binary mask of a circular arc (single-pixel width).

# Arguments
- `Nx`, `Ny`: Grid dimensions
- `cx`, `cy`: Center point of the arc's curvature (1-based grid indices)
- `radius`: Radius of curvature in grid points
- `diameter`: Arc diameter (aperture) in grid points
- `focus`: Focus position `(fx, fy)` — determines arc orientation

# Returns
`BitMatrix` of size (Nx, Ny) with `true` on the arc.
"""
function make_arc(Nx::Int, Ny::Int, cx::Int, cy::Int, radius::Real,
                  diameter::Real, focus::Tuple{Int,Int};
                  plot_arc::Bool=false)
    arc = falses(Nx, Ny)

    fx, fy = focus
    # Direction from center to focus
    dx = fx - cx
    dy = fy - cy
    norm_d = sqrt(dx^2 + dy^2)
    if norm_d ≈ 0
        return arc
    end
    # Center angle (direction of arc)
    center_angle = atan(dy, dx)

    # Half-angle subtended by the aperture
    half_angle = asin(clamp(diameter / (2 * radius), -1.0, 1.0))

    # Number of points to rasterize
    n_points = max(200, ceil(Int, 4 * radius * 2 * half_angle / (2π)))

    for i in 0:n_points
        θ = center_angle - half_angle + i * 2 * half_angle / n_points
        px = round(Int, cx + radius * cos(θ))
        py = round(Int, cy + radius * sin(θ))
        if 1 <= px <= Nx && 1 <= py <= Ny
            arc[px, py] = true
        end
    end

    return arc
end

function make_arc(Nx::Int, Ny::Int, cx::Real, cy::Real, radius::Real,
                  diameter::Real, focus::Tuple{<:Real,<:Real}; kwargs...)
    return make_arc(Nx, Ny, round(Int, cx), round(Int, cy), radius, diameter,
                    (round(Int, focus[1]), round(Int, focus[2])); kwargs...)
end

"""
    make_line(Nx, Ny, startpoint, endpoint; plot_line=false)

Create a 2D binary mask of a line segment (single-pixel width).

Uses Bresenham's line algorithm.

# Arguments
- `Nx`, `Ny`: Grid dimensions
- `startpoint`: `(x1, y1)` start point (1-based grid indices)
- `endpoint`: `(x2, y2)` end point (1-based grid indices)

# Returns
`BitMatrix` of size (Nx, Ny) with `true` on the line.
"""
function make_line(Nx::Int, Ny::Int, startpoint::Tuple{Int,Int}, endpoint::Tuple{Int,Int};
                   plot_line::Bool=false)
    line = falses(Nx, Ny)

    x1, y1 = startpoint
    x2, y2 = endpoint

    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = x1 < x2 ? 1 : -1
    sy = y1 < y2 ? 1 : -1
    err = dx - dy

    x, y = x1, y1
    while true
        if 1 <= x <= Nx && 1 <= y <= Ny
            line[x, y] = true
        end
        if x == x2 && y == y2
            break
        end
        e2 = 2 * err
        if e2 > -dy
            err -= dy
            x += sx
        end
        if e2 < dx
            err += dx
            y += sy
        end
    end

    return line
end

"""
    make_multi_arc(Nx, Ny, arcs; plot_arcs=false)

Create a 2D binary mask of multiple circular arcs.

# Arguments
- `Nx`, `Ny`: Grid dimensions
- `arcs`: Vector of tuples `(cx, cy, radius, diameter, focus)` for each arc

# Returns
`BitMatrix` of size (Nx, Ny) with `true` on all arcs.
"""
function make_multi_arc(Nx::Int, Ny::Int,
                        arcs::Vector{<:Tuple}; plot_arcs::Bool=false)
    mask = falses(Nx, Ny)
    for arc_params in arcs
        cx, cy, radius, diameter, focus = arc_params
        arc_mask = make_arc(Nx, Ny, cx, cy, radius, diameter, focus)
        mask .|= arc_mask
    end
    return mask
end

# ============================================================================
# 3D shapes
# ============================================================================

"""
    make_ball(Nx, Ny, Nz, cx, cy, cz, radius; plot_ball=false)

Create a 3D binary mask of a filled sphere (ball).

# Arguments
- `Nx`, `Ny`, `Nz`: Grid dimensions
- `cx`, `cy`, `cz`: Center point (1-based grid indices)
- `radius`: Sphere radius in grid points

# Returns
`BitArray{3}` of size (Nx, Ny, Nz) with `true` inside the sphere.
"""
function make_ball(Nx::Int, Ny::Int, Nz::Int, cx::Int, cy::Int, cz::Int, radius::Real;
                   plot_ball::Bool=false)
    ball = falses(Nx, Ny, Nz)
    r2 = radius^2
    for k in 1:Nz
        for j in 1:Ny
            for i in 1:Nx
                if (i - cx)^2 + (j - cy)^2 + (k - cz)^2 <= r2
                    ball[i, j, k] = true
                end
            end
        end
    end
    return ball
end

function make_ball(Nx::Int, Ny::Int, Nz::Int, cx::Real, cy::Real, cz::Real, radius::Real; kwargs...)
    return make_ball(Nx, Ny, Nz, round(Int, cx), round(Int, cy), round(Int, cz), radius; kwargs...)
end

"""
    make_sphere(Nx, Ny, Nz, cx, cy, cz, radius; plot_sphere=false)

Create a 3D binary mask of a spherical shell (surface only, single-voxel thickness).

# Arguments
- `Nx`, `Ny`, `Nz`: Grid dimensions
- `cx`, `cy`, `cz`: Center point (1-based grid indices)
- `radius`: Sphere radius in grid points

# Returns
`BitArray{3}` of size (Nx, Ny, Nz) with `true` on the sphere surface.
"""
function make_sphere(Nx::Int, Ny::Int, Nz::Int, cx::Int, cy::Int, cz::Int, radius::Real;
                     plot_sphere::Bool=false)
    sphere = falses(Nx, Ny, Nz)
    r2_inner = (radius - 0.5)^2
    r2_outer = (radius + 0.5)^2
    for k in 1:Nz
        for j in 1:Ny
            for i in 1:Nx
                d2 = (i - cx)^2 + (j - cy)^2 + (k - cz)^2
                if r2_inner <= d2 <= r2_outer
                    sphere[i, j, k] = true
                end
            end
        end
    end
    return sphere
end

function make_sphere(Nx::Int, Ny::Int, Nz::Int, cx::Real, cy::Real, cz::Real, radius::Real; kwargs...)
    return make_sphere(Nx, Ny, Nz, round(Int, cx), round(Int, cy), round(Int, cz), radius; kwargs...)
end

"""
    make_bowl(Nx, Ny, Nz, center, radius, diameter, focus; plot_bowl=false)

Create a 3D binary mask of a bowl (spherical section) surface.

# Arguments
- `Nx`, `Ny`, `Nz`: Grid dimensions
- `center`: `(cx, cy, cz)` center of the bowl's curvature
- `radius`: Radius of curvature in grid points
- `diameter`: Bowl aperture diameter in grid points
- `focus`: `(fx, fy, fz)` focus position — determines bowl orientation

# Returns
`BitArray{3}` with `true` on the bowl surface.
"""
function make_bowl(Nx::Int, Ny::Int, Nz::Int,
                   center::NTuple{3,Int}, radius::Real, diameter::Real,
                   focus::NTuple{3,Int};
                   plot_bowl::Bool=false)
    cx, cy, cz = center
    fx, fy, fz = focus

    bowl = falses(Nx, Ny, Nz)

    # Direction from center to focus
    dir = Float64.([fx - cx, fy - cy, fz - cz])
    norm_dir = norm(dir)
    if norm_dir ≈ 0
        return bowl
    end
    dir ./= norm_dir

    # Half-angle subtended by aperture
    half_angle = asin(clamp(diameter / (2 * radius), -1.0, 1.0))

    r2_inner = (radius - 0.5)^2
    r2_outer = (radius + 0.5)^2

    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        d2 = (i - cx)^2 + (j - cy)^2 + (k - cz)^2
        if r2_inner <= d2 <= r2_outer
            # Check if point is within the cone defined by half_angle
            r_vec = Float64.([i - cx, j - cy, k - cz])
            r_mag = sqrt(d2)
            cos_angle = dot(r_vec, dir) / r_mag
            if cos_angle >= cos(half_angle)
                bowl[i, j, k] = true
            end
        end
    end

    return bowl
end

function make_bowl(Nx::Int, Ny::Int, Nz::Int,
                   center::NTuple{3,<:Real}, radius::Real, diameter::Real,
                   focus::NTuple{3,<:Real}; kwargs...)
    return make_bowl(Nx, Ny, Nz,
                     (round(Int, center[1]), round(Int, center[2]), round(Int, center[3])),
                     radius, diameter,
                     (round(Int, focus[1]), round(Int, focus[2]), round(Int, focus[3]));
                     kwargs...)
end

"""
    make_multi_bowl(Nx, Ny, Nz, bowls; plot_bowls=false)

Create a 3D binary mask of multiple bowl surfaces.

# Arguments
- `Nx`, `Ny`, `Nz`: Grid dimensions
- `bowls`: Vector of tuples `(center, radius, diameter, focus)` for each bowl

# Returns
`BitArray{3}` with `true` on all bowl surfaces.
"""
function make_multi_bowl(Nx::Int, Ny::Int, Nz::Int,
                         bowls::Vector{<:Tuple}; plot_bowls::Bool=false)
    mask = falses(Nx, Ny, Nz)
    for bowl_params in bowls
        center, radius, diameter, focus = bowl_params
        bowl_mask = make_bowl(Nx, Ny, Nz, center, radius, diameter, focus)
        mask .|= bowl_mask
    end
    return mask
end

"""
    make_spherical_section(Nx, Ny, Nz, center, radius, height; plot_section=false)

Create a 3D binary mask of a spherical cap (section).

The cap extends from the sphere surface inward by `height` grid points
along the z-axis from the top of the sphere.

# Arguments
- `Nx`, `Ny`, `Nz`: Grid dimensions
- `center`: `(cx, cy, cz)` center of sphere
- `radius`: Sphere radius in grid points
- `height`: Height of the spherical cap in grid points

# Returns
`BitArray{3}` with `true` within the spherical cap.
"""
function make_spherical_section(Nx::Int, Ny::Int, Nz::Int,
                                center::NTuple{3,Int}, radius::Real, height::Real;
                                plot_section::Bool=false)
    cx, cy, cz = center
    section = falses(Nx, Ny, Nz)
    r2 = radius^2

    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        d2 = (i - cx)^2 + (j - cy)^2 + (k - cz)^2
        if d2 <= r2
            # Check if point is within cap height from top of sphere
            z_offset = k - cz
            if z_offset >= radius - height
                section[i, j, k] = true
            end
        end
    end

    return section
end

function make_spherical_section(Nx::Int, Ny::Int, Nz::Int,
                                center::NTuple{3,<:Real}, radius::Real, height::Real; kwargs...)
    return make_spherical_section(Nx, Ny, Nz,
                                  (round(Int, center[1]), round(Int, center[2]), round(Int, center[3])),
                                  radius, height; kwargs...)
end

"""
    make_line(Nx, Ny, Nz, startpoint, endpoint; plot_line=false)

Create a 3D binary mask of a line segment (single-voxel width).

Uses 3D Bresenham's line algorithm.

# Arguments
- `Nx`, `Ny`, `Nz`: Grid dimensions
- `startpoint`: `(x1, y1, z1)` start point
- `endpoint`: `(x2, y2, z2)` end point

# Returns
`BitArray{3}` with `true` on the line.
"""
function make_line(Nx::Int, Ny::Int, Nz::Int,
                   startpoint::NTuple{3,Int}, endpoint::NTuple{3,Int};
                   plot_line::Bool=false)
    line = falses(Nx, Ny, Nz)

    x1, y1, z1 = startpoint
    x2, y2, z2 = endpoint

    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    dz = abs(z2 - z1)
    sx = x1 < x2 ? 1 : -1
    sy = y1 < y2 ? 1 : -1
    sz = z1 < z2 ? 1 : -1

    # Driving axis is the one with the largest delta
    dm = max(dx, dy, dz)

    x, y, z = x1, y1, z1
    if dm == 0
        if 1 <= x <= Nx && 1 <= y <= Ny && 1 <= z <= Nz
            line[x, y, z] = true
        end
        return line
    end

    # Use parametric stepping along the driving axis
    ex = dm ÷ 2
    ey = dm ÷ 2
    ez = dm ÷ 2

    for _ in 0:dm
        if 1 <= x <= Nx && 1 <= y <= Ny && 1 <= z <= Nz
            line[x, y, z] = true
        end

        ex -= dx
        ey -= dy
        ez -= dz

        if ex < 0
            ex += dm
            x += sx
        end
        if ey < 0
            ey += dm
            y += sy
        end
        if ez < 0
            ez += dm
            z += sz
        end
    end

    return line
end
