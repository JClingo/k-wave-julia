# ============================================================================
# KWave.jl — k-Wave colormap
# ============================================================================

"""
    get_color_map(; num_colors=256)

Return the k-Wave perceptually uniform colormap.

This is a diverging colormap going from cool blue through white to warm red,
designed for visualizing pressure fields (negative/positive values).

The colormap data is based on the MATLAB k-Wave `getColorMap.m` function.

# Arguments
- `num_colors`: Number of colormap entries (default: 256)

# Returns
Vector of (R, G, B) tuples, each component in [0, 1].
"""
function get_color_map(; num_colors::Int=256)
    # Anchor colors for the k-Wave colormap (from MATLAB getColorMap.m)
    # Format: (R, G, B) at normalized positions along the colormap
    # This is a blue-white-red diverging scheme
    anchors = [
        (0.0, (0.0, 0.0, 0.5)),    # dark blue
        (0.15, (0.0, 0.0, 1.0)),   # blue
        (0.35, (0.0, 0.8, 1.0)),   # cyan
        (0.5, (1.0, 1.0, 1.0)),    # white (center)
        (0.65, (1.0, 0.8, 0.0)),   # yellow
        (0.85, (1.0, 0.0, 0.0)),   # red
        (1.0, (0.5, 0.0, 0.0)),    # dark red
    ]

    cmap = Vector{Tuple{Float64, Float64, Float64}}(undef, num_colors)

    for i in 1:num_colors
        t = (i - 1) / (num_colors - 1)

        # Find the two anchor points that bracket t
        idx = 1
        for j in 1:length(anchors)-1
            if t >= anchors[j][1]
                idx = j
            end
        end

        t0, (r0, g0, b0) = anchors[idx]
        t1, (r1, g1, b1) = anchors[min(idx + 1, length(anchors))]

        # Linear interpolation between anchors
        if t1 > t0
            s = (t - t0) / (t1 - t0)
        else
            s = 0.0
        end
        s = clamp(s, 0.0, 1.0)

        r = r0 + s * (r1 - r0)
        g = g0 + s * (g1 - g0)
        b = b0 + s * (b1 - b0)

        cmap[i] = (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
    end

    return cmap
end
