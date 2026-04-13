# ============================================================================
# KWave.jl — Visualization plot function interfaces
# ============================================================================
# These define the interface. Actual rendering is provided by Makie extensions.

"""
    beam_plot(field; kwargs...)

Display a beam pattern plot with orthogonal plane slicing.

For 2D fields: shows a heatmap with optional dB scaling.
For 3D fields: shows three orthogonal slices through the field maximum.

# Arguments
- `field`: 2D or 3D pressure/intensity field
- `db_scale`: If true, display in dB (default: false)
- `db_range`: Dynamic range in dB when using dB scale (default: 40)
- `slice_dim`: For 3D, which dimension to slice (default: auto at max)
- `colormap`: Colormap to use (default: k-Wave colormap)

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function beam_plot(field::AbstractArray;
                   db_scale::Bool=false,
                   db_range::Real=40,
                   slice_dim::Union{Nothing,Int}=nothing,
                   colormap=nothing)
    @warn "beam_plot requires a Makie backend (GLMakie or CairoMakie). Load one first."
    return nothing
end

"""
    fly_through(volume; kwargs...)

Interactive slider-controlled slice viewer for 3D volumes.

# Arguments
- `volume`: 3D array
- `dim`: Dimension to slice through (default: 3)
- `colormap`: Colormap (default: k-Wave colormap)
- `db_scale`: Display in dB (default: false)

# Returns
Plot figure with slider control (or `nothing` without a Makie backend).
"""
function fly_through(volume::AbstractArray{<:Real, 3};
                     dim::Int=3,
                     colormap=nothing,
                     db_scale::Bool=false)
    @warn "fly_through requires a Makie backend (GLMakie). Load it first."
    return nothing
end

"""
    voxel_plot(field; kwargs...)

Volume rendering of a 3D pressure/intensity field.

# Arguments
- `field`: 3D array
- `mode`: Rendering mode — `:volume` (default) or `:mip`
- `db_scale`: Display in dB (default: false)
- `db_range`: Dynamic range in dB (default: 30)
- `dx`, `dy`, `dz`: Voxel spacings for axis labels

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function voxel_plot(field::AbstractArray{<:Real, 3};
                    mode::Symbol=:volume,
                    db_scale::Bool=false,
                    db_range::Real=30,
                    dx::Real=1.0, dy::Real=1.0, dz::Real=1.0)
    @warn "voxel_plot requires a Makie backend (GLMakie). Load it first."
    return nothing
end

"""
    isosurface_plot(field, threshold; kwargs...)

Isosurface rendering of a 3D field at a given threshold value.

# Arguments
- `field`: 3D array
- `threshold`: Isosurface value
- `alpha`: Surface transparency (default: 0.8)
- `dx`, `dy`, `dz`: Voxel spacings for axis labels

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function isosurface_plot(field::AbstractArray{<:Real, 3}, threshold::Real;
                         alpha::Real=0.8,
                         dx::Real=1.0, dy::Real=1.0, dz::Real=1.0)
    @warn "isosurface_plot requires a Makie backend (GLMakie). Load it first."
    return nothing
end

"""
    max_intensity_projection(field; kwargs...)

Maximum intensity projection of a 3D field along all three axes.

# Arguments
- `field`: 3D array
- `colormap`: Colormap (default: k-Wave colormap)

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function max_intensity_projection(field::AbstractArray{<:Real, 3};
                                  colormap=nothing)
    @warn "max_intensity_projection requires a Makie backend (GLMakie or CairoMakie). Load one first."
    return nothing
end

"""
    overlay_plot(background, overlay; kwargs...)

Display an overlay of two fields with transparency.

# Arguments
- `background`: Background field (typically tissue/medium map)
- `overlay`: Overlay field (typically pressure/intensity)
- `alpha`: Overlay transparency (default: 0.5)
- `background_cmap`: Colormap for background (default: :grays)
- `overlay_cmap`: Colormap for overlay (default: k-Wave colormap)

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function overlay_plot(background::AbstractMatrix, overlay::AbstractMatrix;
                      alpha::Real=0.5,
                      background_cmap=nothing,
                      overlay_cmap=nothing)
    @warn "overlay_plot requires a Makie backend (GLMakie or CairoMakie). Load one first."
    return nothing
end

"""
    stacked_plot(signals; kwargs...)

Display multiple time-domain signals stacked vertically with offsets.

# Arguments
- `signals`: Matrix of signals (num_signals × num_time_steps) or vector of vectors
- `dt`: Time step [s] (default: 1.0)
- `labels`: Signal labels (default: nothing)
- `spacing`: Vertical spacing between signals (default: :auto)

# Returns
Plot figure (or `nothing` without a Makie backend).
"""
function stacked_plot(signals::AbstractMatrix;
                      dt::Real=1.0,
                      labels::Union{Nothing, Vector{String}}=nothing,
                      spacing::Union{Symbol, Real}=:auto)
    @warn "stacked_plot requires a Makie backend (GLMakie or CairoMakie). Load one first."
    return nothing
end
