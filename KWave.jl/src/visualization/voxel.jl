# ============================================================================
# KWave.jl — Voxel rendering and advanced 3D visualization
# ============================================================================
# Interface definitions for voxel_plot and 3D visualization.
# Actual rendering provided by Makie extensions.
# ============================================================================

"""
    voxel_plot(volume; kwargs...)

Render a 3D voxel plot of a volumetric field using mesh scatter or volume rendering.

For isosurface rendering, extracts surfaces at specified threshold levels.
For volume rendering, uses transfer functions to map values to color and opacity.

# Arguments
- `volume`: 3D array of field values

# Keyword Arguments
- `mode`: Rendering mode — `:volume`, `:isosurface`, or `:scatter` (default: `:volume`)
- `threshold`: Isosurface threshold (default: 0.5 * maximum)
- `colormap`: Colormap (default: k-Wave colormap)
- `opacity`: Global opacity multiplier (default: 1.0)
- `db_scale`: Display in dB (default: false)
- `db_range`: Dynamic range in dB (default: 40)
- `dx`, `dy`, `dz`: Grid spacings for correct aspect ratio (default: 1.0)

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function voxel_plot(volume::AbstractArray{<:Real, 3};
                    mode::Symbol=:volume,
                    threshold::Union{Nothing, Real}=nothing,
                    colormap=nothing,
                    opacity::Real=1.0,
                    db_scale::Bool=false,
                    db_range::Real=40,
                    dx::Real=1.0,
                    dy::Real=1.0,
                    dz::Real=1.0)
    @warn "voxel_plot requires a Makie backend (GLMakie). Load it first."
    return nothing
end

"""
    isosurface_plot(volume, threshold; kwargs...)

Extract and render an isosurface from a 3D volume.

# Arguments
- `volume`: 3D array
- `threshold`: Isosurface level

# Keyword Arguments
- `colormap`: Colormap (default: k-Wave colormap)
- `alpha`: Surface opacity (default: 0.8)
- `dx`, `dy`, `dz`: Grid spacings

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function isosurface_plot(volume::AbstractArray{<:Real, 3},
                         threshold::Real;
                         colormap=nothing,
                         alpha::Real=0.8,
                         dx::Real=1.0,
                         dy::Real=1.0,
                         dz::Real=1.0)
    @warn "isosurface_plot requires a Makie backend (GLMakie). Load it first."
    return nothing
end

"""
    max_intensity_projection(volume; dims=3, db_scale=false, db_range=40)

Compute maximum intensity projection (MIP) of a 3D volume.

# Arguments
- `volume`: 3D array

# Keyword Arguments
- `dims`: Projection dimension(s) — 1, 2, 3, or :all (default: 3)
- `db_scale`: Apply dB scaling (default: false)
- `db_range`: Dynamic range in dB (default: 40)

# Returns
2D projection(s). For `:all`, returns a named tuple `(xy=..., xz=..., yz=...)`.
"""
function max_intensity_projection(volume::AbstractArray{<:Real, 3};
                                  dims::Union{Int, Symbol}=3,
                                  db_scale::Bool=false,
                                  db_range::Real=40)
    if dims == :all
        mip_xy = dropdims(maximum(abs, volume; dims=3); dims=3)
        mip_xz = dropdims(maximum(abs, volume; dims=2); dims=2)
        mip_yz = dropdims(maximum(abs, volume; dims=1); dims=1)
        if db_scale
            mip_xy = _apply_db_scale(mip_xy, db_range)
            mip_xz = _apply_db_scale(mip_xz, db_range)
            mip_yz = _apply_db_scale(mip_yz, db_range)
        end
        return (xy=mip_xy, xz=mip_xz, yz=mip_yz)
    else
        mip = dropdims(maximum(abs, volume; dims=dims); dims=dims)
        if db_scale
            mip = _apply_db_scale(mip, db_range)
        end
        return mip
    end
end

"""Apply dB scaling to a field, clamping to db_range below maximum."""
function _apply_db_scale(field::AbstractArray, db_range::Real)
    field_max = maximum(abs, field)
    if field_max == 0
        return zeros(Float64, size(field))
    end
    db_field = 20 .* log10.(abs.(field) ./ field_max .+ eps(Float64))
    return max.(db_field, -db_range)
end
