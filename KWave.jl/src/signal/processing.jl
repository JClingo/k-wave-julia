# ============================================================================
# KWave.jl — Signal processing functions
# ============================================================================

"""
    add_noise(signal, snr_db; mode=:rms)

Add Gaussian noise to a signal at a specified signal-to-noise ratio.

# Arguments
- `signal`: Input signal vector
- `snr_db`: Desired SNR in decibels
- `mode`: How to compute signal level — `:rms` (default) or `:peak`

# Returns
Signal with added noise.
"""
function add_noise(signal::AbstractVector, snr_db::Real; mode::Symbol=:rms)
    if mode == :peak
        signal_level = maximum(abs, signal)
    else
        signal_level = sqrt(mean(signal.^2))
    end
    noise_level = signal_level / 10^(snr_db / 20)
    noise = noise_level .* randn(length(signal))
    return signal .+ noise
end

"""
    create_cw_signals(kgrid, source_freq, source_strength, source_phase; num_cycles=nothing)

Generate continuous wave source signals for each source point.

# Arguments
- `kgrid`: KWaveGrid
- `source_freq`: Source frequency [Hz]
- `source_strength`: Source amplitude(s) — scalar or vector (one per source)
- `source_phase`: Source phase(s) [rad] — scalar or vector
- `num_cycles`: Number of ramp-up cycles (default: nothing = no ramp)

# Returns
Matrix of CW signals (num_sources × Nt).
"""
function create_cw_signals(kgrid::AbstractKWaveGrid, source_freq::Real,
                           source_strength, source_phase;
                           num_cycles::Union{Nothing,Int}=nothing)
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    n_sources = max(length(source_strength), length(source_phase))
    strengths = source_strength isa Real ? fill(Float64(source_strength), n_sources) : Float64.(source_strength)
    phases = source_phase isa Real ? fill(Float64(source_phase), n_sources) : Float64.(source_phase)

    signals = zeros(Float64, n_sources, Nt)
    t = (0:Nt-1) .* dt

    for s in 1:n_sources
        @. signals[s, :] = strengths[s] * sin(2π * source_freq * t + phases[s])
    end

    # Apply ramp-up window if requested
    if num_cycles !== nothing && num_cycles > 0
        ramp_length = round(Int, num_cycles / source_freq / dt)
        ramp_length = min(ramp_length, Nt)
        ramp = (1 .- cos.(π .* (0:ramp_length-1) ./ ramp_length)) ./ 2
        for s in 1:n_sources
            signals[s, 1:ramp_length] .*= ramp
        end
    end

    return signals
end

"""
    extract_amp_phase(data, Fs, source_freq; dim=1, fft_padding=1, window=:tukey)

Extract amplitude and phase from a CW signal at a specified frequency.

# Arguments
- `data`: Time-domain signal — vector or matrix (each row/column is a signal)
- `Fs`: Sample frequency [Hz]
- `source_freq`: Frequency to extract [Hz]
- `dim`: Dimension along which time varies (default: 1 for columns)
- `fft_padding`: FFT zero-padding factor (default: 1)
- `window`: Window type to apply (default: `:tukey`)

# Returns
`(amplitude, phase)` — arrays or scalars of amplitude and phase at `source_freq`.
"""
function extract_amp_phase(data::AbstractVector, Fs::Real, source_freq::Real;
                           fft_padding::Int=1, window::Symbol=:tukey)
    N = length(data)
    Nfft = N * fft_padding

    # Apply window
    win = get_win(N, window)
    windowed = data .* win

    # FFT
    S = fft(windowed, Nfft)

    # Find frequency bin closest to source_freq
    freqs = (0:Nfft-1) .* (Fs / Nfft)
    idx = argmin(abs.(freqs .- source_freq))

    # Scale amplitude (compensate for windowing and one-sided spectrum)
    win_correction = N / sum(win)
    amplitude = 2 * abs(S[idx]) / N * win_correction
    phase = angle(S[idx])

    return amplitude, phase
end

function extract_amp_phase(data::AbstractMatrix, Fs::Real, source_freq::Real;
                           dim::Int=2, kwargs...)
    if dim == 2
        n_signals = size(data, 1)
        amps = zeros(Float64, n_signals)
        phases = zeros(Float64, n_signals)
        for i in 1:n_signals
            amps[i], phases[i] = extract_amp_phase(data[i, :], Fs, source_freq; kwargs...)
        end
    else
        n_signals = size(data, 2)
        amps = zeros(Float64, n_signals)
        phases = zeros(Float64, n_signals)
        for i in 1:n_signals
            amps[i], phases[i] = extract_amp_phase(data[:, i], Fs, source_freq; kwargs...)
        end
    end
    return amps, phases
end

"""
    log_compression(signal, compression_ratio)

Apply log compression to a signal for display purposes.

    output = log10(1 + compression_ratio * |signal| / max(|signal|)) / log10(1 + compression_ratio)

# Arguments
- `signal`: Input signal/image
- `compression_ratio`: Compression ratio (higher = more compression)

# Returns
Log-compressed signal with values in [0, 1].
"""
function log_compression(signal::AbstractArray, compression_ratio::Real)
    sig_max = maximum(abs, signal)
    if sig_max ≈ 0
        return zeros(Float64, size(signal))
    end
    normalized = abs.(signal) ./ sig_max
    return log10.(1 .+ compression_ratio .* normalized) ./ log10(1 + compression_ratio)
end

"""
    envelope_detection(signal; dim=1)

Extract the envelope of a signal using the Hilbert transform.

# Arguments
- `signal`: Input signal vector or matrix
- `dim`: Dimension along which to compute (for matrices)

# Returns
Signal envelope (magnitude of the analytic signal).
"""
function envelope_detection(signal::AbstractVector)
    N = length(signal)
    S = fft(signal)

    # Create one-sided spectrum multiplier
    h = zeros(Float64, N)
    if iseven(N)
        h[1] = 1.0          # DC
        h[N÷2+1] = 1.0      # Nyquist
        h[2:N÷2] .= 2.0     # Positive frequencies
    else
        h[1] = 1.0
        h[2:(N+1)÷2] .= 2.0
    end

    analytic = ifft(S .* h)
    return abs.(analytic)
end

function envelope_detection(signal::AbstractMatrix; dim::Int=1)
    if dim == 1
        result = similar(signal, Float64)
        for j in axes(signal, 2)
            result[:, j] = envelope_detection(signal[:, j])
        end
        return result
    else
        result = similar(signal, Float64)
        for i in axes(signal, 1)
            result[i, :] = envelope_detection(signal[i, :])
        end
        return result
    end
end

"""
    gradient_fd(f, dx)

Compute the gradient of a field using second-order finite differences.

# Arguments
- `f`: Input field (1D, 2D, or 3D)
- `dx`: Grid spacing(s) — scalar or tuple

# Returns
For 1D: gradient vector. For 2D/3D: tuple of gradient arrays.
"""
function gradient_fd(f::AbstractVector, dx::Real)
    N = length(f)
    df = similar(f, Float64)
    # Central differences for interior, forward/backward at edges
    if N >= 3
        df[1] = (f[2] - f[1]) / dx
        df[N] = (f[N] - f[N-1]) / dx
        for i in 2:N-1
            df[i] = (f[i+1] - f[i-1]) / (2 * dx)
        end
    elseif N == 2
        df[1] = (f[2] - f[1]) / dx
        df[2] = df[1]
    else
        df[1] = 0.0
    end
    return df
end

function gradient_fd(f::AbstractMatrix, dx::NTuple{2,<:Real})
    Nx, Ny = size(f)
    dfx = similar(f, Float64)
    dfy = similar(f, Float64)

    for j in 1:Ny
        dfx[:, j] = gradient_fd(f[:, j], dx[1])
    end
    for i in 1:Nx
        dfy[i, :] = gradient_fd(f[i, :], dx[2])
    end

    return dfx, dfy
end

"""
    gradient_spect(f, kgrid)

Compute the gradient of a field using spectral (FFT) differentiation.

# Arguments
- `f`: Input field array
- `kgrid`: KWaveGrid

# Returns
For 1D: gradient vector. For 2D: tuple `(df/dx, df/dy)`.
"""
function gradient_spect(f::AbstractVector, kgrid::KWaveGrid1D)
    F = fft(f)
    dfdx = real.(ifft(im .* kgrid.kx_vec .* F))
    return dfdx
end

function gradient_spect(f::AbstractMatrix, kgrid::KWaveGrid2D)
    F = fft(f)
    kx = kgrid.kx_vec
    ky = reshape(kgrid.ky_vec, 1, :)
    dfdx = real.(ifft(im .* kx .* F))
    dfdy = real.(ifft(im .* ky .* F))
    return dfdx, dfdy
end

"""
    filter_time_series(kgrid, medium, signal; filter_type=:lowpass, zero_phase=false, order=64)

Filter a time series signal based on the CFL number and grid properties.
Uses the maximum supported frequency from the grid as the cutoff.

# Arguments
- `kgrid`: KWaveGrid
- `medium`: KWaveMedium (for sound speed)
- `signal`: Input time series
- `filter_type`: `:lowpass` (default), `:highpass`, or `:bandpass`
- `zero_phase`: Apply zero-phase filtering (default: false)
- `order`: FIR filter order (default: 64)

# Returns
Filtered signal.
"""
function filter_time_series(kgrid::AbstractKWaveGrid, medium::KWaveMedium,
                            signal::AbstractVector;
                            filter_type::Symbol=:lowpass,
                            zero_phase::Bool=false, order::Int=64)
    dt = kgrid.dt[]
    Fs = 1.0 / dt
    c_max = medium.sound_speed isa Real ? medium.sound_speed : maximum(medium.sound_speed)
    # Maximum supported frequency from grid
    f_max = c_max * k_max(kgrid) / (2π)
    cutoff = min(f_max, Fs / 2 * 0.95)  # Stay below Nyquist
    return apply_filter(signal, Fs, cutoff, filter_type; zero_phase=zero_phase, order=order)
end

"""
    spect(signal, Fs; plot_spectrum=false)

Compute the single-sided amplitude and phase spectrum of a signal.

# Arguments
- `signal`: Input time-domain signal
- `Fs`: Sample frequency [Hz]

# Returns
`(freqs, amplitude, phase)` — frequency vector, amplitude spectrum, and phase spectrum.
"""
function spect(signal::AbstractVector, Fs::Real; plot_spectrum::Bool=false)
    N = length(signal)
    S = fft(signal)

    # Single-sided spectrum
    if iseven(N)
        n_freq = N ÷ 2 + 1
    else
        n_freq = (N + 1) ÷ 2
    end

    freqs = (0:n_freq-1) .* (Fs / N)
    amplitude = abs.(S[1:n_freq]) .* (2.0 / N)
    amplitude[1] /= 2  # DC component
    if iseven(N)
        amplitude[end] /= 2  # Nyquist component
    end
    phase = angle.(S[1:n_freq])

    return freqs, amplitude, phase
end
