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

end # module KWaveGLMakieExt
