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
    db2neper(alpha_coeff, alpha_power)

Convert power-law absorption coefficient from dB/(MHz^y cm) to Nepers/((rad/s)^y m).

    α_Np = α_dB * (1e-2 / (20 * log10(e))) * (2π * 1e6)^(-y)

# Arguments
- `alpha_coeff`: Absorption coefficient [dB/(MHz^y cm)]
- `alpha_power`: Power law exponent y

# Returns
Absorption coefficient in [Nepers/((rad/s)^y m)].
"""
db2neper(alpha_coeff, alpha_power) = alpha_coeff * (1e-2 / (20 * log10(ℯ))) * (2π * 1e6)^(-alpha_power)

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
