# ============================================================================
# KWave.jl — Geometry / shape creation functions
# ============================================================================

"""
    make_disc(Nx, Ny, cx, cy, radius; plot_disc=false)

Create a 2D binary mask of a filled disc (circle).

# Arguments
- `Nx`, `Ny`: Grid dimensions
- `cx`, `cy`: Center point (1-based grid indices)
- `radius`: Disc radius in grid points
- `plot_disc`: If true, display the disc (requires Makie extension)

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

# Convenience method accepting Float64 center
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
- `plot_circle`: If true, display the circle

# Returns
`BitMatrix` of size (Nx, Ny) with `true` on the circle perimeter.
"""
function make_circle(Nx::Int, Ny::Int, cx::Int, cy::Int, radius::Int;
                     arc_angle::Float64=2π, plot_circle::Bool=false)
    circle = falses(Nx, Ny)

    if arc_angle >= 2π
        # Full circle: use midpoint circle algorithm
        _midpoint_circle!(circle, Nx, Ny, cx, cy, radius)
    else
        # Partial arc: parametric approach
        # Estimate number of points needed
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
        # Set all 8 octant symmetric points
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
