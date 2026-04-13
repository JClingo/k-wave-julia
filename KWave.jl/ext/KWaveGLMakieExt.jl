# ============================================================================
# KWave.jl GLMakie Extension — Real-time simulation visualization
# ============================================================================

module KWaveGLMakieExt

using KWave
using GLMakie

"""
    GLMakieDisplay <: KWave.SimulationDisplay

Real-time simulation display using GLMakie.
"""
mutable struct GLMakieDisplay <: KWave.SimulationDisplay
    fig::Figure
    ax::Axis
    obs_field::Observable{Matrix{Float64}}
    obs_title::Observable{String}
    hm::Any
    plot_scale::Union{Symbol, Tuple}
end

function KWave.create_sim_display(grid::KWave.KWaveGrid2D;
                                   plot_layout::Symbol=:default,
                                   plot_scale::Union{Symbol, Tuple}=:auto)
    Nx, Ny = grid.Nx, grid.Ny

    fig = Figure(size=(700, 600))

    obs_field = Observable(zeros(Nx, Ny))
    obs_title = Observable("k-Wave 2D: t = 0 / 0")

    ax = Axis(fig[1, 1]; title=obs_title,
              xlabel="x", ylabel="y",
              aspect=DataAspect())

    cmap_tuples = KWave.get_color_map()
    cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]

    if plot_scale isa Tuple
        clims = plot_scale
    else
        clims = (-1.0, 1.0)
    end

    hm = heatmap!(ax, obs_field; colormap=cmap, colorrange=clims)
    Colorbar(fig[1, 2], hm)

    display(fig)

    return GLMakieDisplay(fig, ax, obs_field, obs_title, hm, plot_scale)
end

function KWave.update_sim_display!(disp::GLMakieDisplay, p_field::Matrix, t_index::Int, Nt::Int)
    disp.obs_field[] = copy(p_field)
    disp.obs_title[] = "k-Wave 2D: t = $t_index / $Nt"

    if disp.plot_scale == :auto || disp.plot_scale == :symmetric
        p_max = maximum(abs, p_field)
        if p_max > 0
            disp.hm.colorrange[] = (-p_max, p_max)
        end
    end

    yield()
    return nothing
end

function KWave.close_sim_display!(disp::GLMakieDisplay)
    return nothing
end

"""
    GLMakieRecorder <: KWave.MovieRecorder

Movie recorder using GLMakie's record functionality.
"""
mutable struct GLMakieRecorder <: KWave.MovieRecorder
    fig::Figure
    ax::Axis
    obs_field::Observable{Matrix{Float64}}
    obs_title::Observable{String}
    hm::Any
    path::String
    frames::Vector{Matrix{Float64}}
    fps::Int
end

function KWave.create_movie_recorder(path::String, grid::KWave.KWaveGrid2D; fps::Int=15)
    Nx, Ny = grid.Nx, grid.Ny

    fig = Figure(size=(700, 600))
    obs_field = Observable(zeros(Nx, Ny))
    obs_title = Observable("k-Wave 2D")

    ax = Axis(fig[1, 1]; title=obs_title, xlabel="x", ylabel="y", aspect=DataAspect())

    cmap_tuples = KWave.get_color_map()
    cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]

    hm = heatmap!(ax, obs_field; colormap=cmap, colorrange=(-1.0, 1.0))
    Colorbar(fig[1, 2], hm)

    return GLMakieRecorder(fig, ax, obs_field, obs_title, hm, path, Matrix{Float64}[], fps)
end

function KWave.record_frame!(rec::GLMakieRecorder, p_field::Matrix, t_index::Int)
    push!(rec.frames, copy(p_field))
    return nothing
end

function KWave.finalize_movie!(rec::GLMakieRecorder)
    if isempty(rec.frames)
        return nothing
    end

    p_max_global = maximum(maximum(abs, f) for f in rec.frames)
    if p_max_global > 0
        rec.hm.colorrange[] = (-p_max_global, p_max_global)
    end

    record(rec.fig, rec.path, enumerate(rec.frames); framerate=rec.fps) do (i, frame)
        rec.obs_field[] = frame
        rec.obs_title[] = "k-Wave 2D: frame $i / $(length(rec.frames))"
    end

    return nothing
end

# ============================================================================
# beam_plot
# ============================================================================

function KWave.beam_plot(field::AbstractMatrix;
                          db_scale::Bool=false,
                          db_range::Real=40,
                          slice_dim::Union{Nothing,Int}=nothing,
                          colormap=nothing)
    fig = Figure(size=(700, 600))

    cmap = if colormap === nothing
        cmap_tuples = KWave.get_color_map()
        [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    else
        colormap
    end

    plot_data = if db_scale
        ref = maximum(abs, field)
        ref = ref > 0 ? ref : 1.0
        clamp.(20 .* log10.(abs.(field) ./ ref .+ eps()), -db_range, 0.0)
    else
        field
    end

    ax = Axis(fig[1, 1]; title="Beam Pattern", xlabel="x", ylabel="y", aspect=DataAspect())
    clims = db_scale ? (-db_range, 0.0) : extrema(plot_data)
    hm = heatmap!(ax, plot_data; colormap=cmap, colorrange=clims)
    Colorbar(fig[1, 2], hm; label=db_scale ? "dB" : "Pressure")

    return fig
end

function KWave.beam_plot(field::AbstractArray{<:Real, 3};
                          db_scale::Bool=false,
                          db_range::Real=40,
                          slice_dim::Union{Nothing,Int}=nothing,
                          colormap=nothing)
    cmap = if colormap === nothing
        cmap_tuples = KWave.get_color_map()
        [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    else
        colormap
    end

    max_idx = argmax(abs.(field))
    ix, iy, iz = Tuple(max_idx)

    fig = Figure(size=(1200, 400))
    slices = [
        (field[ix, :, :], "YZ (x=$ix)"),
        (field[:, iy, :], "XZ (y=$iy)"),
        (field[:, :, iz], "XY (z=$iz)"),
    ]

    for (i, (slice_data, title)) in enumerate(slices)
        plot_data = if db_scale
            ref = maximum(abs, field)
            ref = ref > 0 ? ref : 1.0
            clamp.(20 .* log10.(abs.(slice_data) ./ ref .+ eps()), -db_range, 0.0)
        else
            slice_data
        end
        ax = Axis(fig[1, i]; title=title, aspect=DataAspect())
        clims = db_scale ? (-db_range, 0.0) : extrema(plot_data)
        heatmap!(ax, plot_data; colormap=cmap, colorrange=clims)
    end

    return fig
end

# ============================================================================
# overlay_plot
# ============================================================================

function KWave.overlay_plot(background::AbstractMatrix, overlay::AbstractMatrix;
                             alpha::Real=0.5,
                             background_cmap=nothing,
                             overlay_cmap=nothing)
    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1]; aspect=DataAspect())

    bg_cmap = background_cmap === nothing ? :grays : background_cmap
    heatmap!(ax, background; colormap=bg_cmap)

    ov_cmap = if overlay_cmap === nothing
        cmap_tuples = KWave.get_color_map()
        [RGBAf(r, g, b, Float32(alpha)) for (r, g, b) in cmap_tuples]
    else
        overlay_cmap
    end
    heatmap!(ax, overlay; colormap=ov_cmap)

    return fig
end

# ============================================================================
# stacked_plot
# ============================================================================

function KWave.stacked_plot(signals::AbstractMatrix;
                             dt::Real=1.0,
                             labels::Union{Nothing, Vector{String}}=nothing,
                             spacing::Union{Symbol, Real}=:auto)
    n_signals, Nt = size(signals)
    t = (0:Nt-1) .* dt

    offset = if spacing == :auto
        maximum(abs, signals) * 1.5
    else
        Float64(spacing)
    end

    fig = Figure(size=(800, max(300, 50 * n_signals)))
    ax = Axis(fig[1, 1]; xlabel="Time", ylabel="Signal")

    for i in 1:n_signals
        y = signals[i, :] .+ (n_signals - i) * offset
        label = labels !== nothing ? labels[i] : "ch $i"
        lines!(ax, t, y; label=label)
    end

    if n_signals <= 20
        axislegend(ax)
    end

    return fig
end

# ============================================================================
# fly_through (interactive slider — GLMakie only)
# ============================================================================

function KWave.fly_through(volume::AbstractArray{<:Real, 3};
                            dim::Int=3,
                            colormap=nothing,
                            db_scale::Bool=false)
    cmap = if colormap === nothing
        cmap_tuples = KWave.get_color_map()
        [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    else
        colormap
    end

    n = size(volume, dim)
    fig = Figure(size=(700, 650))
    ax = Axis(fig[1, 1]; title="Slice View", aspect=DataAspect())
    sl = Slider(fig[2, 1]; range=1:n, startvalue=n ÷ 2)

    slice_obs = @lift begin
        i = $(sl.value)
        raw = if dim == 1
            volume[i, :, :]
        elseif dim == 2
            volume[:, i, :]
        else
            volume[:, :, i]
        end
        if db_scale
            ref = maximum(abs, volume)
            ref = ref > 0 ? ref : 1.0
            clamp.(20 .* log10.(abs.(raw) ./ ref .+ eps()), -60.0, 0.0)
        else
            raw
        end
    end

    heatmap!(ax, slice_obs; colormap=cmap)
    Label(fig[2, 2], @lift("Slice $( $(sl.value) ) / $n"); tellwidth=false)

    return fig
end

# ============================================================================
# voxel_plot
# ============================================================================

function KWave.voxel_plot(field::AbstractArray{<:Real, 3};
                           mode::Symbol=:volume,
                           db_scale::Bool=false,
                           db_range::Real=30,
                           dx::Real=1.0, dy::Real=1.0, dz::Real=1.0)
    plot_data = if db_scale
        ref = maximum(abs, field)
        ref = ref > 0 ? ref : 1.0
        norm = clamp.((20 .* log10.(abs.(field) ./ ref .+ eps()) .+ db_range) ./ db_range, 0.0, 1.0)
        Float32.(norm)
    else
        fmax = maximum(abs, field)
        fmax > 0 ? Float32.(abs.(field) ./ fmax) : Float32.(abs.(field))
    end

    Nx, Ny, Nz = size(field)
    fig = Figure(size=(700, 600))
    ax = Axis3(fig[1, 1]; title="Volume Render",
               xlabel="x ($(round(Nx*dx*1e3, digits=1)) mm)",
               ylabel="y ($(round(Ny*dy*1e3, digits=1)) mm)",
               zlabel="z ($(round(Nz*dz*1e3, digits=1)) mm)")

    cmap_tuples = KWave.get_color_map()
    cmap = [RGBAf(r, g, b, t) for (t, (r, g, b)) in zip(range(0, 1, length(cmap_tuples)), cmap_tuples)]

    volume!(ax, plot_data; algorithm=:absorptionrgba, colormap=cmap)

    return fig
end

# ============================================================================
# isosurface_plot
# ============================================================================

function KWave.isosurface_plot(field::AbstractArray{<:Real, 3}, threshold::Real;
                                alpha::Real=0.8,
                                dx::Real=1.0, dy::Real=1.0, dz::Real=1.0)
    Nx, Ny, Nz = size(field)
    fig = Figure(size=(700, 600))
    ax = Axis3(fig[1, 1]; title="Isosurface at $(round(threshold, sigdigits=3))",
               xlabel="x", ylabel="y", zlabel="z")

    cmap_tuples = KWave.get_color_map()
    cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    isorange = Float32(max(abs(threshold) * 0.05, 1e-6))

    volume!(ax, Float32.(field);
            algorithm=:iso,
            isovalue=Float32(threshold),
            isorange=isorange,
            colormap=cmap,
            alpha=Float32(alpha))

    return fig
end

# ============================================================================
# max_intensity_projection
# ============================================================================

function KWave._plot_max_intensity_projection(field::AbstractArray{<:Real, 3};
                                              colormap=nothing)
    cmap = if colormap === nothing
        cmap_tuples = KWave.get_color_map()
        [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    else
        colormap
    end

    mip_xy = dropdims(maximum(field; dims=3); dims=3)
    mip_xz = dropdims(maximum(field; dims=2); dims=2)
    mip_yz = dropdims(maximum(field; dims=1); dims=1)

    fig = Figure(size=(1200, 400))
    for (i, (proj, title)) in enumerate([(mip_xy, "XY (max z)"), (mip_xz, "XZ (max y)"), (mip_yz, "YZ (max x)")])
        ax = Axis(fig[1, i]; title=title, aspect=DataAspect())
        clims = extrema(proj)
        clims = clims[1] == clims[2] ? (clims[1] - 1, clims[2] + 1) : clims
        hm = heatmap!(ax, proj; colormap=cmap, colorrange=clims)
        i == 3 && Colorbar(fig[1, 4], hm)
    end

    return fig
end

end # module KWaveGLMakieExt
