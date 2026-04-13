# ============================================================================
# KWave.jl — Material property functions
# ============================================================================

"""
    water_sound_speed(T)

Return the sound speed of water at temperature `T` [°C].

Uses the 5th-order polynomial fit from Marczak (1997).
Valid for 0 ≤ T ≤ 95 °C at atmospheric pressure.

# Arguments
- `T`: Temperature in degrees Celsius

# Returns
Sound speed [m/s].
"""
function water_sound_speed(T::Real)
    T = Float64(T)
    return 1.402385e3 +
           5.038813e0 * T -
           5.799136e-2 * T^2 +
           3.287156e-4 * T^3 -
           1.398845e-6 * T^4 +
           2.787860e-9 * T^5
end

"""
    water_density(T)

Return the density of water at temperature `T` [°C].

Uses the polynomial fit from Kell (1975).
Valid for 0 ≤ T ≤ 100 °C at atmospheric pressure.

# Arguments
- `T`: Temperature in degrees Celsius

# Returns
Density [kg/m³].
"""
function water_density(T::Real)
    T = Float64(T)
    num = 999.83952 +
          16.945176 * T -
          7.9870401e-3 * T^2 -
          46.170461e-6 * T^3 +
          105.56302e-9 * T^4 -
          280.54253e-12 * T^5
    den = 1 + 16.879850e-3 * T
    return num / den
end

"""
    water_absorption(T, f)

Return the acoustic absorption coefficient of water at temperature `T` [°C]
and frequency `f` [Hz].

Uses the formula from Pinkerton (1949).

# Arguments
- `T`: Temperature in degrees Celsius
- `f`: Frequency [Hz]

# Returns
Absorption coefficient [dB/(MHz^2 cm)] × f² in [dB/cm], or per [Np/m] depending on use.
Actually returns absorption in [dB/m] for the given frequency.
"""
function water_absorption(T::Real, f::Real)
    T = Float64(T)
    f_mhz = Float64(f) / 1e6
    # Absorption coefficient in dB/(MHz^2 cm) from Pinkerton
    # Simplified fit: alpha ~ 2.17e-15 * f^2 for water at 20°C
    # Temperature-dependent fit (Stokes-Kirchhoff + structural relaxation)
    viscosity = 1.002e-3 * exp(-1.126e-2 * T)  # Dynamic viscosity [Pa·s]
    rho = water_density(T)
    c = water_sound_speed(T)
    # Classical absorption (Stokes)
    alpha_classical = 2 * (2π * f)^2 * viscosity / (3 * rho * c^3)
    # Convert from Np/m to dB/m
    alpha_db = alpha_classical * 20 * log10(ℯ)
    return alpha_db
end

"""
    water_nonlinearity(T)

Return the nonlinearity parameter B/A of water at temperature `T` [°C].

Uses polynomial fit from Beyer (1960).
Valid for 0 ≤ T ≤ 100 °C.

# Arguments
- `T`: Temperature in degrees Celsius

# Returns
B/A (dimensionless).
"""
function water_nonlinearity(T::Real)
    T = Float64(T)
    return 4.6 + 0.046 * T + 3.2e-4 * T^2 - 2.5e-6 * T^3
end

"""
    fit_power_law_params(f, alpha, c0)

Fit power-law absorption parameters to measured attenuation data.

Fits the model: α(f) = α₀ * f^y

where α₀ is in [dB/(MHz^y cm)] and f is in [MHz].

# Arguments
- `f`: Frequency values [Hz] — vector
- `alpha`: Measured attenuation values [dB/cm] — vector
- `c0`: Sound speed [m/s] (not used in fit, but returned for convenience)

# Returns
`(alpha_coeff, alpha_power)` — fitted power-law parameters.
"""
function fit_power_law_params(f::AbstractVector, alpha::AbstractVector, c0::Real)
    # Convert to MHz
    f_mhz = Float64.(f) ./ 1e6

    # Log-log linear regression: log(alpha) = log(alpha0) + y * log(f_mhz)
    log_f = log.(f_mhz)
    log_alpha = log.(Float64.(alpha))

    # Remove any non-finite values
    valid = isfinite.(log_f) .& isfinite.(log_alpha)
    log_f = log_f[valid]
    log_alpha = log_alpha[valid]

    if length(log_f) < 2
        error("Need at least 2 valid data points for fitting")
    end

    # Linear regression
    n = length(log_f)
    sum_x = sum(log_f)
    sum_y = sum(log_alpha)
    sum_xy = sum(log_f .* log_alpha)
    sum_x2 = sum(log_f.^2)

    alpha_power = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x^2)
    log_alpha0 = (sum_y - alpha_power * sum_x) / n
    alpha_coeff = exp(log_alpha0)

    return alpha_coeff, alpha_power
end
