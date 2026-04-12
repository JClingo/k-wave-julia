# ============================================================================
# KWave.jl — Simulation display interface
# ============================================================================
# This file defines the interface for real-time simulation monitoring.
# The actual rendering is provided by package extensions (GLMakie, CairoMakie).
# Without a visualization backend loaded, these are no-ops.

"""
    SimulationDisplay

Abstract type for simulation display backends.
"""
abstract type SimulationDisplay end

"""
    NullDisplay <: SimulationDisplay

No-op display when no visualization backend is loaded.
"""
struct NullDisplay <: SimulationDisplay end

"""
    create_sim_display(grid, plot_layout, plot_scale)

Create a simulation display for real-time monitoring.
Returns a `NullDisplay` unless a visualization backend (GLMakie) is loaded.

Override this function in a package extension to provide actual rendering.
"""
function create_sim_display(grid::AbstractKWaveGrid;
                            plot_layout::Symbol=:default,
                            plot_scale::Union{Symbol, Tuple}=:auto)
    return NullDisplay()
end

"""
    update_sim_display!(display, p_field, t_index, Nt)

Update the simulation display with the current pressure field.
No-op for `NullDisplay`.
"""
function update_sim_display!(display::NullDisplay, p_field, t_index::Int, Nt::Int)
    return nothing
end

"""
    close_sim_display!(display)

Close the simulation display window.
"""
function close_sim_display!(display::NullDisplay)
    return nothing
end

"""
    MovieRecorder

Abstract type for movie recording backends.
"""
abstract type MovieRecorder end

"""
    NullRecorder <: MovieRecorder

No-op recorder when no visualization backend is loaded.
"""
struct NullRecorder <: MovieRecorder end

"""
    create_movie_recorder(path, grid; fps=15)

Create a movie recorder for saving simulation frames.
Returns a `NullRecorder` unless a visualization backend is loaded.
"""
function create_movie_recorder(path::String, grid::AbstractKWaveGrid; fps::Int=15)
    return NullRecorder()
end

"""
    record_frame!(recorder, p_field, t_index)

Record a single frame of the simulation.
"""
function record_frame!(recorder::NullRecorder, p_field, t_index::Int)
    return nothing
end

"""
    finalize_movie!(recorder)

Finalize and save the movie file.
"""
function finalize_movie!(recorder::NullRecorder)
    return nothing
end
