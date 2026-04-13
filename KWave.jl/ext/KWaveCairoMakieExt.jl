# ============================================================================
# KWave.jl CairoMakie Extension — Publication-quality static plots
# ============================================================================

module KWaveCairoMakieExt

using KWave
using CairoMakie

# ============================================================================
# beam_plot
# ============================================================================

function KWave.beam_plot(field::AbstractMatrix;
                          db_scale::Bool=false,
                          db_range::Real=40,
                          slice_dim::Union{Nothing,Int}=nothing,
                          colormap=nothing)
    fig = Figure(size=(700, 600))

    if colormap === nothing
        cmap_tuples = KWave.get_color_map()
        cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    else
        cmap = colormap
    end

    plot_data = if db_scale
        ref = maximum(abs, field)
        ref = ref > 0 ? ref : 1.0
        db_vals = 20 .* log10.(abs.(field) ./ ref .+ eps())
        clamp.(db_vals, -db_range, 0.0)
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
    if colormap === nothing
        cmap_tuples = KWave.get_color_map()
        cmap = [RGBf(r, g, b) for (r, g, b) in cmap_tuples]
    else
        cmap = colormap
    end

    # Find maximum location for default slicing
    max_idx = argmax(abs.(field))
    ix, iy, iz = Tuple(max_idx)

    fig = Figure(size=(1200, 400))

    slices = [
        (field[ix, :, :], "YZ (x=$ix)", 1),
        (field[:, iy, :], "XZ (y=$iy)", 2),
        (field[:, :, iz], "XY (z=$iz)", 3),
    ]

    for (i, (slice_data, title, _)) in enumerate(slices)
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

end # module KWaveCairoMakieExt
