# ============================================================================
# KWave.jl — Shared types, enums, and type aliases
# ============================================================================

"""
    SourceMode

Mode for time-varying sources.

- `Dirichlet`: Replace field value at source points
- `Additive`: Add to field value with correction for staggered grid
- `AdditiveNoCorrection`: Add to field value without correction
"""
@enum SourceMode Dirichlet Additive AdditiveNoCorrection

"""
Type alias for fields that can be either a scalar or an array (e.g., medium properties).
"""
const ScalarOrArray{T} = Union{T, AbstractArray{T}}
