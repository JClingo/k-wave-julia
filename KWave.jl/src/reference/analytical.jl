# ============================================================================
# KWave.jl — Analytical reference solutions
# ============================================================================

"""
    focused_annulus_oneil(radius_outer, radius_inner, amplitude, phase, freq, c0, density, axial_positions)

Compute the O'Neil solution for the on-axis pressure from a focused annular transducer.

Based on O'Neil, H.T., "Theory of Focusing Radiators," JASA (1949).

# Arguments
- `radius_outer`: Outer radius of the annulus [m]
- `radius_inner`: Inner radius of the annulus [m]
- `amplitude`: Source pressure amplitude [Pa]
- `phase`: Source phase [rad]
- `freq`: Source frequency [Hz]
- `c0`: Sound speed [m/s]
- `density`: Medium density [kg/m³]
- `axial_positions`: Vector of axial positions [m]

# Returns
Complex pressure along the axis.
"""
function focused_annulus_oneil(
    radius_outer::Real,
    radius_inner::Real,
    amplitude::Real,
    phase::Real,
    freq::Real,
    c0::Real,
    density::Real,
    axial_positions::AbstractVector;
    focus_distance::Real=Inf,
)
    k = 2π * freq / c0
    omega = 2π * freq

    p_axis = zeros(ComplexF64, length(axial_positions))

    for (i, z) in enumerate(axial_positions)
        if abs(z) < eps(Float64)
            p_axis[i] = 0.0
            continue
        end

        # Distance from outer edge to field point on axis
        r_outer = sqrt(z^2 + radius_outer^2)
        # Distance from inner edge to field point on axis
        r_inner = sqrt(z^2 + radius_inner^2)

        if focus_distance == Inf
            # Unfocused: flat piston annulus
            # p(z) = ρ*c*u0 * (exp(ikz) - exp(ik*r_outer) + exp(ik*r_inner) - exp(ikz))
            # Simplified using Rayleigh integral
            p_axis[i] = density * c0 * amplitude * exp(im * phase) * (
                exp(im * k * z) - exp(im * k * r_outer) +
                exp(im * k * r_inner) - exp(im * k * z)
            )
            # For flat annulus: p(z) = ρcu₀[exp(ikz) - exp(ik·r_outer)]  + ρcu₀[exp(ik·r_inner) - exp(ikz)]
            # Actually the standard O'Neil result for a flat annulus:
            phase_outer = exp(im * k * r_outer)
            phase_inner = exp(im * k * r_inner)
            phase_z = exp(im * k * z)

            p_axis[i] = density * c0 * amplitude * exp(im * phase) * (
                phase_z - phase_outer + phase_inner - phase_z
            )
        else
            # Focused annulus — geometric focusing
            # Delay-corrected phases for curved element
            delay_outer = (sqrt(focus_distance^2 + radius_outer^2) - focus_distance) / c0
            delay_inner = (sqrt(focus_distance^2 + radius_inner^2) - focus_distance) / c0

            r_outer_focus = sqrt(z^2 + radius_outer^2)
            r_inner_focus = sqrt(z^2 + radius_inner^2)

            p_outer = exp(im * k * (r_outer_focus - c0 * delay_outer))
            p_inner = exp(im * k * (r_inner_focus - c0 * delay_inner))
            p_z_outer = exp(im * k * (z - c0 * delay_outer))
            p_z_inner = exp(im * k * (z - c0 * delay_inner))

            p_axis[i] = density * c0 * amplitude * exp(im * phase) * (
                p_z_outer - p_outer + p_inner - p_z_inner
            )
        end
    end

    return p_axis
end

"""
    focused_bowl_oneil(radius, diameter, amplitude, phase, freq, c0, density, axial_positions)

Compute the O'Neil solution for the on-axis pressure from a focused spherical bowl transducer.

# Arguments
- `radius`: Radius of curvature [m]
- `diameter`: Bowl aperture diameter [m]
- `amplitude`: Source velocity amplitude [m/s]
- `phase`: Source phase [rad]
- `freq`: Source frequency [Hz]
- `c0`: Sound speed [m/s]
- `density`: Medium density [kg/m³]
- `axial_positions`: Vector of axial positions [m]

# Returns
Complex pressure along the axis.
"""
function focused_bowl_oneil(
    radius::Real,
    diameter::Real,
    amplitude::Real,
    phase::Real,
    freq::Real,
    c0::Real,
    density::Real,
    axial_positions::AbstractVector,
)
    k = 2π * freq / c0
    a = diameter / 2  # aperture radius

    # Bowl geometry
    # h = R - sqrt(R² - a²) is the bowl height
    R = Float64(radius)
    h = R - sqrt(R^2 - a^2)

    p_axis = zeros(ComplexF64, length(axial_positions))

    for (i, z) in enumerate(axial_positions)
        if abs(z) < eps(Float64)
            p_axis[i] = 0.0
            continue
        end

        # Distance from bowl rim to axial point
        r1 = sqrt((z - h)^2 + a^2)
        # Distance from bowl apex to axial point
        r2 = abs(z)

        # O'Neil formula for focused bowl
        # p(z) = ρcu₀R/(z) * [exp(ikr2) - exp(ikr1)]
        # with geometric focusing correction
        p_axis[i] = im * density * c0 * amplitude * exp(im * phase) * R / z * (
            exp(im * k * r2) - exp(im * k * r1)
        )
    end

    return p_axis
end

"""
    mendousse(x, t, source_amplitude, freq, c0, density, BonA, alpha_0;
              num_harmonics=50)

Compute the Mendousse solution for nonlinear plane wave propagation.

Provides an analytical solution for the propagation of a finite-amplitude
sinusoidal plane wave in a thermoviscous fluid.

Based on Mendousse, "Nonlinear Dissipative Distortion of Progressive Sound
Waves at Moderate Amplitudes," JASA (1953).

# Arguments
- `x`: Propagation distance(s) [m]
- `t`: Time points [s]
- `source_amplitude`: Source velocity amplitude [m/s]
- `freq`: Source frequency [Hz]
- `c0`: Small-signal sound speed [m/s]
- `density`: Medium density [kg/m³]
- `BonA`: Nonlinearity parameter B/A
- `alpha_0`: Absorption coefficient [Np/m] at the source frequency

# Keyword Arguments
- `num_harmonics`: Number of harmonics in the Fubini-Blackstock series (default: 50)

# Returns
Pressure waveform `p(t)` at distance `x` (or matrix for vector `x`).
"""
function mendousse(
    x::Union{Real, AbstractVector},
    t::AbstractVector,
    source_amplitude::Real,
    freq::Real,
    c0::Real,
    density::Real,
    BonA::Real,
    alpha_0::Real;
    num_harmonics::Int=50,
)
    omega = 2π * freq
    k = omega / c0
    beta = 1 + BonA / 2  # coefficient of nonlinearity
    u0 = Float64(source_amplitude)
    p0 = density * c0 * u0

    # Goldberg number: ratio of nonlinear to dissipative effects
    # Γ = β * k * u0 / (α₀ * c0)
    # Shock distance: x_s = 1 / (β * k * Ma) where Ma = u0/c0
    x_shock = c0 / (beta * omega * u0 / c0)

    x_vec = x isa Real ? [Float64(x)] : Float64.(x)

    result = zeros(Float64, length(x_vec), length(t))

    for (ix, xval) in enumerate(x_vec)
        # Normalized distance (σ = x / x_shock)
        sigma = xval / x_shock

        # Fubini-Blackstock solution (weak shock regime)
        # Valid when σ < 1 (before shock) and with absorption
        for (it, tval) in enumerate(t)
            retarded_time = tval - xval / c0
            theta = omega * retarded_time

            p_val = 0.0
            for n in 1:num_harmonics
                # Bessel function coefficient for nth harmonic
                Jn = _bessel_fubini(n, sigma)

                # Absorption at nth harmonic
                atten = exp(-n^2 * alpha_0 * xval)

                p_val += (2 / (n * sigma)) * Jn * atten * sin(n * theta)
            end

            result[ix, it] = p0 * p_val
        end
    end

    return x isa Real ? result[1, :] : result
end

"""
Compute the Fubini coefficient (2/nσ) * J_n(nσ) for the Fubini solution.
Uses the series expansion of the Bessel function J_n(nσ).
"""
function _bessel_fubini(n::Int, sigma::Real)
    # J_n(nσ) via series expansion
    x = n * sigma
    if abs(x) < 1e-15
        return n == 1 ? 0.5 : 0.0
    end

    # Compute Bessel J_n(x) using the series definition
    term = (x / 2)^n / factorial(big(n))
    result = Float64(term)
    for m in 1:30
        term *= -(x / 2)^2 / (m * (n + m))
        result += Float64(term)
        if abs(Float64(term)) < 1e-15 * abs(result)
            break
        end
    end

    return result
end
