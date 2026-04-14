# ============================================================================
# KWave.jl WGLMakie Extension — Web-based dashboard prototype
# ============================================================================

module KWaveWGLMakieExt

using KWave
using WGLMakie
using Bonito

function _kwave_colormap(; alpha::Union{Nothing, Real}=nothing)
    cmap_tuples = KWave.get_color_map()
    if alpha === nothing
        return [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    end
    return [RGBAf(r, g, b, Float32(alpha)) for (r, g, b) in cmap_tuples]
end

"""
    KWaveDashboard

Web-based simulation monitoring dashboard using WGLMakie + Bonito.

Provides a browser-accessible interface for:
- Launching and monitoring simulations
- Viewing field snapshots and time series
- Exporting results
"""
mutable struct KWaveDashboard
    app::Bonito.App
    port::Int
end

"""
    create_dashboard(; port=9284)

Create and start a web-based simulation dashboard.

# Arguments
- `port`: HTTP port to serve on (default: 9284)

# Returns
`KWaveDashboard` instance. Open http://localhost:<port> in a browser.
"""
function create_dashboard(; port::Int=9284)
    app = Bonito.App() do session
        fig = Figure(size=(900, 700))

        # Status panel
        status_text = Observable("Ready")

        # Field display
        field_obs = Observable(zeros(Float32, 64, 64))
        ax = Axis(fig[1, 1]; title="Pressure Field", aspect=DataAspect())

        cmap_tuples = KWave.get_color_map()
        cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
        heatmap!(ax, field_obs; colormap=cmap, colorrange=(-1.0, 1.0))

        # Info panel
        info_label = Label(fig[2, 1], status_text; fontsize=14)

        return Bonito.DOM.div(
            Bonito.DOM.h2("KWave.jl Simulation Dashboard"),
            fig,
        )
    end

    server = Bonito.Server(app, "0.0.0.0", port)

    println("KWave Dashboard running at http://localhost:$port")

    return KWaveDashboard(app, port)
end

# ============================================================================
# WGLMakie-based simulation display
# ============================================================================

mutable struct WGLMakieDisplay <: KWave.SimulationDisplay
    fig::Figure
    ax::Axis
    obs_field::Observable{Matrix{Float64}}
    obs_title::Observable{String}
    hm::Any
    plot_scale::Union{Symbol, Tuple}
end

function KWave.create_sim_display(grid::KWave.KWaveGrid2D,
                                   ::Val{:wglmakie};
                                   plot_layout::Symbol=:default,
                                   plot_scale::Union{Symbol, Tuple}=:auto)
    Nx, Ny = grid.Nx, grid.Ny

    fig = Figure(size=(700, 600))
    obs_field = Observable(zeros(Nx, Ny))
    obs_title = Observable("k-Wave 2D: t = 0 / 0")

    ax = Axis(fig[1, 1]; title=obs_title, xlabel="x", ylabel="y", aspect=DataAspect())

    cmap_tuples = KWave.get_color_map()
    cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]

    clims = plot_scale isa Tuple ? plot_scale : (-1.0, 1.0)
    hm = heatmap!(ax, obs_field; colormap=cmap, colorrange=clims)
    Colorbar(fig[1, 2], hm)

    return WGLMakieDisplay(fig, ax, obs_field, obs_title, hm, plot_scale)
end

function KWave.update_sim_display!(disp::WGLMakieDisplay, p_field::Matrix, t_index::Int, Nt::Int)
    disp.obs_field[] = copy(p_field)
    disp.obs_title[] = "k-Wave 2D: t = $t_index / $Nt"

    if disp.plot_scale == :auto || disp.plot_scale == :symmetric
        p_max = maximum(abs, p_field)
        if p_max > 0
            disp.hm.colorrange[] = (-p_max, p_max)
        end
    end
    return nothing
end

function KWave.close_sim_display!(disp::WGLMakieDisplay)
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
    cmap = colormap === nothing ? _kwave_colormap() : colormap

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
    cmap = colormap === nothing ? _kwave_colormap() : colormap

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

    ov_cmap = overlay_cmap === nothing ? _kwave_colormap(alpha=alpha) : overlay_cmap
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

end # module KWaveWGLMakieExt
