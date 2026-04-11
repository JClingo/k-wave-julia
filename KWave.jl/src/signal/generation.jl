# ============================================================================
# KWave.jl — Signal generation functions
# ============================================================================

"""
    tone_burst(sample_freq, signal_freq, num_cycles; envelope=:gaussian, signal_offset=0, plot_signal=false)

Generate a windowed tone burst signal.

# Arguments
- `sample_freq`: Sampling frequency [Hz]
- `signal_freq`: Signal frequency [Hz]
- `num_cycles`: Number of sinusoidal cycles in the burst
- `envelope`: Window envelope type — `:gaussian` (default), `:rectangular`, or `:cosine`
- `signal_offset`: Number of zero-padding samples before the burst (default: 0)
- `plot_signal`: If true, display the signal

# Returns
`Vector{Float64}` containing the tone burst signal.

# Example
```julia
signal = tone_burst(40e6, 1e6, 5)  # 5-cycle 1 MHz burst at 40 MHz sample rate
```
"""
function tone_burst(sample_freq::Real, signal_freq::Real, num_cycles::Int;
                    envelope::Symbol=:gaussian, signal_offset::Int=0,
                    plot_signal::Bool=false)
    # Period and number of samples per burst
    tone_length = num_cycles / signal_freq
    dt = 1.0 / sample_freq
    N = round(Int, tone_length / dt)

    # Time vector for the burst
    t = (0:N-1) .* dt

    # Generate carrier
    carrier = sin.(2π * signal_freq .* t)

    # Generate envelope
    if envelope == :gaussian
        # Gaussian window centered on the burst
        t_center = tone_length / 2
        sigma = tone_length / (2 * num_cycles)  # width scaled to cycles
        env = exp.(-((t .- t_center) ./ sigma).^2 / 2)
    elseif envelope == :rectangular
        env = ones(Float64, N)
    elseif envelope == :cosine
        # Cosine taper (Tukey-like)
        env = zeros(Float64, N)
        for i in 1:N
            env[i] = 0.5 * (1.0 - cos(2π * (i - 1) / (N - 1)))
        end
    else
        error("Unknown envelope type: $envelope. Use :gaussian, :rectangular, or :cosine")
    end

    # Apply envelope to carrier
    signal = carrier .* env

    # Add zero-padding offset
    if signal_offset > 0
        signal = vcat(zeros(Float64, signal_offset), signal)
    end

    return signal
end

"""
    gaussian_pulse(x; magnitude=1.0, mean=0.0, variance=1.0)

Compute a Gaussian function.

    f(x) = magnitude * exp(-((x - mean)^2) / (2 * variance))

# Arguments
- `x`: Input value(s) — scalar or array
- `magnitude`: Peak amplitude (default: 1.0)
- `mean`: Center position (default: 0.0)
- `variance`: Variance (width squared) (default: 1.0)

# Returns
Gaussian function value(s), same type as input `x`.
"""
function gaussian_pulse(x; magnitude::Real=1.0, mean::Real=0.0, variance::Real=1.0)
    return magnitude .* exp.(-((x .- mean).^2) ./ (2 .* variance))
end
