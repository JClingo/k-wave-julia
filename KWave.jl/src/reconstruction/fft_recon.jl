# ============================================================================
# KWave.jl — FFT-based reconstruction functions
# ============================================================================

"""
    kspace_line_recon(p, dy, dt; c=1500.0, interp_fac=1, pos_cond=true)

1D FFT-based reconstruction from a planar (line) sensor.

Reconstructs the initial pressure distribution from time-series data
recorded by a line of sensors, using the k-space method.

Based on the algorithm from Treeby et al., "Photoacoustic tomography in
absorbing acoustic media using time reversal" (2010).

# Arguments
- `p`: Sensor data matrix (num_sensors × num_time_steps)
- `dy`: Sensor spacing [m]
- `dt`: Time step [s]
- `c`: Sound speed [m/s] (default: 1500)
- `interp_fac`: Interpolation factor for output (default: 1)
- `pos_cond`: Apply positivity condition (default: true)

# Returns
Reconstructed pressure distribution (2D matrix: depth × lateral).
"""
function kspace_line_recon(p::AbstractMatrix, dy::Real, dt::Real;
                           c::Real=1500.0, interp_fac::Int=1,
                           pos_cond::Bool=true)
    Ny, Nt = size(p)
    c = Float64(c)
    dy = Float64(dy)
    dt = Float64(dt)

    # Compute wavenumber and angular frequency vectors
    dky = 2π / (Ny * dy)
    dw = 2π / (Nt * dt)

    ky = _fft_freq_vec(Ny) .* dky
    w = _fft_freq_vec(Nt) .* dw

    # 2D FFT of sensor data
    P = fft(p)

    # Compute kx from dispersion relation: kx = sqrt((w/c)^2 - ky^2)
    # Using broadcasting for the 2D arrays
    ky_2d = reshape(ky, :, 1)
    w_2d = reshape(w, 1, :)

    kx_sq = (w_2d ./ c).^2 .- ky_2d.^2

    # Only propagating modes (kx² > 0)
    kx = sqrt.(max.(kx_sq, 0.0))

    # Phase shift to reconstruct at x = 0 (sensor plane)
    # The reconstruction maps temporal frequency to spatial frequency
    # Scale by kx for the Jacobian of the coordinate transform
    scaling = zeros(ComplexF64, Ny, Nt)
    for j in 1:Nt, i in 1:Ny
        if kx[i, j] > 0
            scaling[i, j] = kx[i, j] * c^2
        end
    end

    P_recon = P .* scaling

    # Inverse FFT to get reconstruction
    p_recon = real.(ifft(P_recon))

    # Interpolation
    if interp_fac > 1
        Ny_new = Ny * interp_fac
        Nt_new = Nt * interp_fac
        p_recon = resize_array(p_recon, (Ny_new, Nt_new))
    end

    # Apply positivity condition
    if pos_cond
        p_recon = max.(p_recon, 0.0)
    end

    return p_recon
end

"""
    kspace_plane_recon(p, dy, dz, dt; c=1500.0, interp_fac=1, pos_cond=true)

2D FFT-based reconstruction from a planar sensor.

Reconstructs the initial pressure distribution from time-series data
recorded by a 2D planar array of sensors.

# Arguments
- `p`: Sensor data (Ny × Nz × Nt) 3D array
- `dy`: Sensor spacing in y [m]
- `dz`: Sensor spacing in z [m]
- `dt`: Time step [s]
- `c`: Sound speed [m/s] (default: 1500)
- `interp_fac`: Interpolation factor (default: 1)
- `pos_cond`: Apply positivity condition (default: true)

# Returns
Reconstructed pressure distribution (3D array: depth × Ny × Nz).
"""
function kspace_plane_recon(p::AbstractArray{<:Real, 3}, dy::Real, dz::Real, dt::Real;
                            c::Real=1500.0, interp_fac::Int=1,
                            pos_cond::Bool=true)
    Ny, Nz, Nt = size(p)
    c = Float64(c)
    dy = Float64(dy)
    dz = Float64(dz)
    dt = Float64(dt)

    # Wavenumber and angular frequency vectors
    dky = 2π / (Ny * dy)
    dkz = 2π / (Nz * dz)
    dw = 2π / (Nt * dt)

    ky = _fft_freq_vec(Ny) .* dky
    kz = _fft_freq_vec(Nz) .* dkz
    w = _fft_freq_vec(Nt) .* dw

    # 3D FFT of sensor data
    P = fft(p)

    # Compute kx from dispersion relation: kx = sqrt((w/c)^2 - ky^2 - kz^2)
    ky_3d = reshape(ky, :, 1, 1)
    kz_3d = reshape(kz, 1, :, 1)
    w_3d = reshape(w, 1, 1, :)

    kx_sq = (w_3d ./ c).^2 .- ky_3d.^2 .- kz_3d.^2
    kx = sqrt.(max.(kx_sq, 0.0))

    # Scale by kx for Jacobian
    scaling = zeros(ComplexF64, Ny, Nz, Nt)
    for k in 1:Nt, j in 1:Nz, i in 1:Ny
        if kx[i, j, k] > 0
            scaling[i, j, k] = kx[i, j, k] * c^2
        end
    end

    P_recon = P .* scaling
    p_recon = real.(ifft(P_recon))

    if pos_cond
        p_recon = max.(p_recon, 0.0)
    end

    return p_recon
end

# ============================================================================
# Helper
# ============================================================================

"""Generate FFT frequency index vector [0, 1, ..., N/2-1, -N/2, ..., -1]."""
function _fft_freq_vec(N::Int)
    if iseven(N)
        return Float64[collect(0:N÷2-1); collect(-N÷2:-1)]
    else
        return Float64[collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)]
    end
end
