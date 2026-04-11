# ============================================================================
# KWave.jl — Unit conversion functions
# ============================================================================

"""
    db2neper(alpha_db)

Convert attenuation from decibels to nepers.

    α_Np = α_dB / (20 * log10(e))

# Arguments
- `alpha_db`: Attenuation in decibels [dB]

# Returns
Attenuation in nepers [Np].
"""
db2neper(alpha_db) = alpha_db ./ (20 * log10(ℯ))

"""
    neper2db(alpha_neper)

Convert attenuation from nepers to decibels.

    α_dB = α_Np * 20 * log10(e)

# Arguments
- `alpha_neper`: Attenuation in nepers [Np]

# Returns
Attenuation in decibels [dB].
"""
neper2db(alpha_neper) = alpha_neper .* (20 * log10(ℯ))
