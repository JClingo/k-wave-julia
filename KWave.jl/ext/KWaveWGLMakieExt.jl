# ============================================================================
# KWave.jl WGLMakie Extension — Web-based dashboard prototype
# ============================================================================

module KWaveWGLMakieExt

using KWave
using WGLMakie
using Bonito

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

end # module KWaveWGLMakieExt
