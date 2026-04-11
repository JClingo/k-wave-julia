# ============================================================================
# KWave.jl — Filtering and window functions
# ============================================================================

"""
    smooth(matrix; restore_max=false)

Smooth an N-dimensional matrix using frequency-domain windowed smoothing.

Applies a Blackman window in the frequency domain along each dimension
to suppress high-frequency content (Gibbs artifacts from sharp boundaries).

# Arguments
- `matrix`: Input array (1D, 2D, or 3D)
- `restore_max`: If true, rescale output to preserve original maximum value

# Returns
Smoothed array of the same size.
"""
function smooth(matrix::AbstractArray; restore_max::Bool=false)
    nd = ndims(matrix)
    result = complex(float(copy(matrix)))
    dims = size(matrix)

    # Apply windowed smoothing in frequency domain
    # Transform, window, inverse transform along each dimension
    for d in 1:nd
        N = dims[d]
        # Create Blackman window for this dimension
        win = _blackman_window(N)

        # Reshape window for broadcasting along dimension d
        shape = ones(Int, nd)
        shape[d] = N
        win_shaped = reshape(win, shape...)

        # Forward FFT along dimension d
        result = fft(result, d)

        # Apply window (centered on zero frequency via fftshift pattern)
        # The window needs to be applied in the same order as FFT output
        win_fft = fftshift(win_shaped, d)
        result .*= win_fft

        # Inverse FFT along dimension d
        result = ifft(result, d)
    end

    smoothed = real.(result)

    if restore_max
        orig_max = maximum(abs.(matrix))
        smooth_max = maximum(abs.(smoothed))
        if smooth_max > 0
            smoothed .*= orig_max / smooth_max
        end
    end

    return smoothed
end

"""
    _blackman_window(N)

Generate a Blackman window of length N.
"""
function _blackman_window(N::Int)
    if N == 1
        return [1.0]
    end
    n = 0:N-1
    return @. 0.42 - 0.5 * cos(2π * n / (N - 1)) + 0.08 * cos(4π * n / (N - 1))
end

"""
    apply_filter(signal, Fs, cutoff_freq, filter_type; zero_phase=false, order=64)

Apply a FIR frequency filter to a time-domain signal.

# Arguments
- `signal`: Input signal vector
- `Fs`: Sample frequency [Hz]
- `cutoff_freq`: Cutoff frequency [Hz], or (low, high) tuple for bandpass
- `filter_type`: `:lowpass`, `:highpass`, or `:bandpass`
- `zero_phase`: If true, apply filter forward and backward (zero phase distortion)
- `order`: FIR filter order (default: 64)

# Returns
Filtered signal vector.
"""
function apply_filter(signal::AbstractVector, Fs::Real, cutoff_freq, filter_type::Symbol;
                      zero_phase::Bool=false, order::Int=64)
    # Normalize frequency (0 to 1, where 1 = Nyquist)
    nyquist = Fs / 2

    if filter_type == :lowpass
        fc_norm = cutoff_freq / nyquist
        responsetype = DSP.Lowpass(fc_norm)
    elseif filter_type == :highpass
        fc_norm = cutoff_freq / nyquist
        responsetype = DSP.Highpass(fc_norm)
    elseif filter_type == :bandpass
        fc_low = cutoff_freq[1] / nyquist
        fc_high = cutoff_freq[2] / nyquist
        responsetype = DSP.Bandpass(fc_low, fc_high)
    else
        error("Unknown filter type: $filter_type. Use :lowpass, :highpass, or :bandpass")
    end

    # Design FIR filter using a Kaiser window
    designmethod = DSP.FIRWindow(DSP.kaiser(order + 1, 5.0))
    filt = digitalfilter(responsetype, designmethod)

    if zero_phase
        return filtfilt(filt, signal)
    else
        return DSP.filt(filt, signal)
    end
end

"""
    gaussian_filter(signal, Fs, center_freq; bandwidth=nothing)

Apply a Gaussian-shaped frequency filter to a signal.

# Arguments
- `signal`: Input signal vector
- `Fs`: Sample frequency [Hz]
- `center_freq`: Center frequency of Gaussian [Hz]
- `bandwidth`: Gaussian bandwidth [Hz] (default: center_freq / 2)

# Returns
Filtered signal vector.
"""
function gaussian_filter(signal::AbstractVector, Fs::Real, center_freq::Real;
                         bandwidth::Union{Nothing, Real}=nothing)
    N = length(signal)
    bw = bandwidth === nothing ? center_freq / 2 : bandwidth

    # Compute FFT
    S = fft(signal)

    # Frequency vector
    freqs = fftfreq(N, Fs)

    # Gaussian window in frequency domain
    gauss_win = exp.(-((freqs .- center_freq).^2 .+ (freqs .+ center_freq).^2) ./ (2 * bw^2))

    # Apply filter and inverse transform
    return real.(ifft(S .* gauss_win))
end

"""
    fftfreq(N, Fs)

Generate frequency vector for FFT output of length N at sample rate Fs.
"""
function fftfreq(N::Int, Fs::Real)
    if iseven(N)
        f = [collect(0:N÷2-1); collect(-N÷2:-1)] .* (Fs / N)
    else
        f = [collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)] .* (Fs / N)
    end
    return f
end

"""
    get_win(N, type; param=nothing, symmetric=true)

Generate a window function of length `N`.

# Arguments
- `N`: Window length
- `type`: Window type — `:hann`, `:hamming`, `:blackman`, `:kaiser`,
  `:tukey`, `:gaussian`, `:bartlett`, `:rectangular`
- `param`: Window parameter (required for `:kaiser` (beta), `:tukey` (alpha),
  `:gaussian` (sigma))
- `symmetric`: If true (default), symmetric window. If false, periodic.

# Returns
`Vector{Float64}` of window values.
"""
function get_win(N::Int, type::Symbol; param=nothing, symmetric::Bool=true)
    if N <= 0
        return Float64[]
    end
    if N == 1
        return [1.0]
    end

    # For periodic windows, compute length N+1 and drop the last point
    M = symmetric ? N : N + 1
    n = collect(0:M-1)

    win = if type == :hann
        @. 0.5 * (1 - cos(2π * n / (M - 1)))
    elseif type == :hamming
        @. 0.54 - 0.46 * cos(2π * n / (M - 1))
    elseif type == :blackman
        @. 0.42 - 0.5 * cos(2π * n / (M - 1)) + 0.08 * cos(4π * n / (M - 1))
    elseif type == :bartlett
        @. 1 - abs(2 * n / (M - 1) - 1)
    elseif type == :rectangular
        ones(Float64, M)
    elseif type == :kaiser
        beta = param === nothing ? 5.0 : Float64(param)
        _kaiser_window(M, beta)
    elseif type == :tukey
        alpha = param === nothing ? 0.5 : Float64(param)
        _tukey_window(M, alpha)
    elseif type == :gaussian
        sigma = param === nothing ? 0.4 : Float64(param)
        @. exp(-0.5 * ((n - (M - 1) / 2) / (sigma * (M - 1) / 2))^2)
    else
        error("Unknown window type: $type. Supported: :hann, :hamming, :blackman, :bartlett, :rectangular, :kaiser, :tukey, :gaussian")
    end

    # For periodic windows, drop the last point
    return symmetric ? win : win[1:N]
end

"""
    _kaiser_window(N, beta)

Generate a Kaiser window of length N with parameter beta.
"""
function _kaiser_window(N::Int, beta::Float64)
    # I₀(x) - Modified Bessel function of the first kind, order 0
    # Julia doesn't have this in base, use series approximation
    n = collect(0:N-1)
    alpha = (N - 1) / 2
    win = [_besseli0(beta * sqrt(1 - ((i - alpha) / alpha)^2)) / _besseli0(beta) for i in n]
    return win
end

"""
Modified Bessel function I₀(x) via series expansion.
"""
function _besseli0(x::Real)
    sum_val = 1.0
    term = 1.0
    for k in 1:25
        term *= (x / (2k))^2
        sum_val += term
        if term < 1e-16 * sum_val
            break
        end
    end
    return sum_val
end

"""
    _tukey_window(N, alpha)

Generate a Tukey (tapered cosine) window.
"""
function _tukey_window(N::Int, alpha::Float64)
    if alpha <= 0
        return ones(Float64, N)
    elseif alpha >= 1
        return get_win(N, :hann)
    end

    win = ones(Float64, N)
    width = floor(Int, alpha * (N - 1) / 2)
    for i in 0:width
        val = 0.5 * (1 + cos(π * (2i / (alpha * (N - 1)) - 1)))
        win[i + 1] = val
        win[N - i] = val
    end
    return win
end
