# ============================================================================
# KWave.jl — Time-stepping operations for the k-space first-order solver
# ============================================================================

# ============================================================================
# Absorption parameters struct
# ============================================================================

"""
    AbsorptionParams

Pre-computed absorption and dispersion parameters for the equation of state.
Following the k-Wave formulation (Treeby & Cox, 2010):
- `absorb_tau`: scalar absorption prefactor
- `absorb_eta`: scalar dispersion prefactor
- `absorb_nabla1`: k-space operator k^y (fractional Laplacian for absorption)
- `absorb_nabla2`: k-space operator k^(y-1) (fractional Laplacian for dispersion)
"""
struct AbsorptionParams
    absorb_tau::Float64           # -2 * alpha_nepers * c_ref^(y-1)
    absorb_eta::Float64           # 2 * alpha_nepers * c_ref^(y-1) * tan(π*y/2)
    absorb_nabla1::AbstractArray  # k^y
    absorb_nabla2::AbstractArray  # k^(y-1)
    mode::Symbol                  # :no_absorption, :no_dispersion, :stokes
end

# ============================================================================
# Staggered grid density interpolation
# ============================================================================

"""
    stagger_density(rho0, dim, N)

Interpolate density to staggered grid points by averaging adjacent values.
Uses periodic boundary conditions (circular shift).
"""
function _stagger_density_1d(rho0::AbstractVector)
    return (rho0 .+ circshift(rho0, -1)) ./ 2
end

function _stagger_density_2d(rho0::AbstractMatrix, dim::Int)
    shift = dim == 1 ? (-1, 0) : (0, -1)
    return (rho0 .+ circshift(rho0, shift)) ./ 2
end

function _stagger_density_3d(rho0::AbstractArray{<:Real, 3}, dim::Int)
    shift = ntuple(i -> i == dim ? -1 : 0, 3)
    return (rho0 .+ circshift(rho0, shift)) ./ 2
end

# ============================================================================
# 1D: Initialize from p0
# ============================================================================

function initialize_p0_1d!(p::Vector{Float64}, ux::Vector{Float64},
                           rhox::Vector{Float64},
                           source::KWaveSource, kgrid::KWaveGrid1D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::Vector{ComplexF64},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    p0 = Float64.(source.p0)

    if smooth_p0 && kgrid.Nx > 4
        p0 = smooth(p0; restore_max=true)
    end

    p .= p0

    c0_sq = medium.sound_speed isa Real ? medium.sound_speed^2 : medium.sound_speed.^2
    @. rhox = p0 / c0_sq

    dt = kgrid.dt[]
    rho0 = medium.density

    dpdx = similar(p)
    spectral_gradient!(dpdx, p0, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 1)

    if rho0 isa Real
        @. ux = -dt / (2 * rho0) * dpdx
    else
        rho0_sgx = _stagger_density_1d(rho0)
        @. ux = -dt / (2 * rho0_sgx) * dpdx
    end
end

# ============================================================================
# 2D: Initialize from p0
# ============================================================================

function initialize_p0_2d!(p::Matrix{Float64}, ux::Matrix{Float64}, uy::Matrix{Float64},
                           rhox::Matrix{Float64}, rhoy::Matrix{Float64},
                           source::KWaveSource, kgrid::KWaveGrid2D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::Matrix{ComplexF64},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    p0 = Float64.(source.p0)

    if smooth_p0 && kgrid.Nx > 1 && kgrid.Ny > 1
        p0 = smooth(p0; restore_max=true)
    end

    p .= p0

    c0_sq = medium.sound_speed isa Real ? medium.sound_speed^2 : medium.sound_speed.^2

    @. rhox = p0 / (2 * c0_sq)
    @. rhoy = p0 / (2 * c0_sq)

    dt = kgrid.dt[]
    rho0 = medium.density

    dpdx = similar(p)
    dpdy = similar(p)
    spectral_gradient!(dpdx, p0, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 2)
    spectral_gradient!(dpdy, p0, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 2)

    if rho0 isa Real
        @. ux = -dt / (2 * rho0) * dpdx
        @. uy = -dt / (2 * rho0) * dpdy
    else
        rho0_sgx = _stagger_density_2d(rho0, 1)
        rho0_sgy = _stagger_density_2d(rho0, 2)
        @. ux = -dt / (2 * rho0_sgx) * dpdx
        @. uy = -dt / (2 * rho0_sgy) * dpdy
    end
end

# ============================================================================
# 3D: Initialize from p0
# ============================================================================

function initialize_p0_3d!(p::Array{Float64,3}, ux::Array{Float64,3},
                           uy::Array{Float64,3}, uz::Array{Float64,3},
                           rhox::Array{Float64,3}, rhoy::Array{Float64,3}, rhoz::Array{Float64,3},
                           source::KWaveSource, kgrid::KWaveGrid3D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::Array{ComplexF64,3},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    p0 = Float64.(source.p0)

    if smooth_p0 && kgrid.Nx > 1 && kgrid.Ny > 1 && kgrid.Nz > 1
        p0 = smooth(p0; restore_max=true)
    end

    p .= p0

    c0_sq = medium.sound_speed isa Real ? medium.sound_speed^2 : medium.sound_speed.^2

    @. rhox = p0 / (3 * c0_sq)
    @. rhoy = p0 / (3 * c0_sq)
    @. rhoz = p0 / (3 * c0_sq)

    dt = kgrid.dt[]
    rho0 = medium.density

    dpdx = similar(p)
    dpdy = similar(p)
    dpdz = similar(p)
    spectral_gradient!(dpdx, p0, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 3)
    spectral_gradient!(dpdy, p0, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 3)
    spectral_gradient!(dpdz, p0, kgrid.kz_vec, kgrid.ddz_k_shift_pos, scratch, plans, 3, 3)

    if rho0 isa Real
        @. ux = -dt / (2 * rho0) * dpdx
        @. uy = -dt / (2 * rho0) * dpdy
        @. uz = -dt / (2 * rho0) * dpdz
    else
        rho0_sgx = _stagger_density_3d(rho0, 1)
        rho0_sgy = _stagger_density_3d(rho0, 2)
        rho0_sgz = _stagger_density_3d(rho0, 3)
        @. ux = -dt / (2 * rho0_sgx) * dpdx
        @. uy = -dt / (2 * rho0_sgy) * dpdy
        @. uz = -dt / (2 * rho0_sgz) * dpdz
    end
end

# ============================================================================
# Equation of state (shared logic)
# ============================================================================

"""
Compute pressure from equation of state with k-space correction,
optional absorption, and optional nonlinearity.

For the absorbing case, uses the k-Wave formulation (Treeby & Cox, 2010)
applied in the frequency domain:
  p_hat = c₀² · κ · (1 + absorb_nabla1) · ρ_hat
where absorb_nabla1 = τ·k^y (pre-multiplied in precomputation).

For dispersion, additionally subtracts:
  c₀² · κ · absorb_nabla2 · (ρ_hat - ρ_hat_prev) / dt

`rho_total_prev` is the total density from the previous time step (for dispersion).
"""
function _compute_pressure!(
    p::AbstractArray,
    rho_total::AbstractArray,
    rho_total_prev::Union{Nothing, AbstractArray},
    scratch1::AbstractArray{<:Complex},
    scratch2::AbstractArray{<:Complex},
    kappa::AbstractArray,
    c0,     # scalar or array
    rho0,   # scalar or array (density for nonlinear term)
    BonA,   # nothing or scalar or array
    absorb::Union{Nothing, AbsorptionParams},
    plans::FFTPlans,
    dt::Float64,
)
    is_absorbing = absorb !== nothing && absorb.mode != :no_absorption
    is_nonlinear = BonA !== nothing

    # Choose the effective density (with nonlinear correction if needed)
    rho_eff = if is_nonlinear
        @. rho_total * (1 + BonA / 2 * rho_total / rho0)
    else
        rho_total
    end

    if c0 isa Real
        c_sq = c0^2

        # FFT of effective density
        scratch1 .= complex.(rho_eff)
        plans.forward * scratch1

        if !is_absorbing
            # Lossless: p = c0^2 * IFFT(kappa * FFT(rho))
            @. scratch1 = c_sq * kappa * scratch1
        elseif absorb.mode == :no_dispersion || rho_total_prev === nothing
            # Absorption only (no dispersion):
            # p = c0^2 * IFFT(kappa * (1 + nabla1) * FFT(rho))
            @. scratch1 = c_sq * kappa * (1 + absorb.absorb_nabla1) * scratch1
        else
            # Full absorption + dispersion:
            # p = c0^2 * IFFT(kappa * ((1 + nabla1) * FFT(rho) - nabla2 * FFT(drho/dt)))
            rho_eff_prev = if is_nonlinear
                @. rho_total_prev * (1 + BonA / 2 * rho_total_prev / rho0)
            else
                rho_total_prev
            end
            scratch2 .= complex.((rho_eff .- rho_eff_prev) ./ dt)
            plans.forward * scratch2
            @. scratch1 = c_sq * kappa * (
                (1 + absorb.absorb_nabla1) * scratch1
                - absorb.absorb_nabla2 * scratch2
            )
        end

        plans.inverse * scratch1
        @. p = real(scratch1)

    else
        c0_sq = c0.^2
        # For heterogeneous media: multiply by c^2 in spatial domain first
        scratch1 .= complex.(c0_sq .* rho_eff)
        plans.forward * scratch1

        if !is_absorbing
            @. scratch1 = kappa * scratch1
        else
            # Use mean c_ref for absorption (already in absorb_nabla1 via precomputation)
            @. scratch1 = kappa * (1 + absorb.absorb_nabla1) * scratch1
        end

        plans.inverse * scratch1
        @. p = real(scratch1)
    end
end

# ============================================================================
# 1D Time stepping
# ============================================================================

function time_step_1d!(
    p::Vector{Float64},
    ux::Vector{Float64},
    rhox::Vector{Float64},
    scratch1::Vector{ComplexF64},
    scratch2::Vector{ComplexF64},
    kgrid::KWaveGrid1D,
    medium::KWaveMedium,
    source::KWaveSource,
    pml_x::Vector{Float64},
    pml_x_sgx::Vector{Float64},
    kappa::Vector{Float64},
    plans::FFTPlans,
    t_index::Int,
    absorb::Union{Nothing, AbsorptionParams},
    rho_total_prev::Union{Nothing, Vector{Float64}}=nothing,
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density

    # === STEP 1: Pressure gradient ===
    dpdx = similar(p)
    spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 1)

    # === STEP 2: Velocity update with PML ===
    if rho0 isa Real
        @. ux = pml_x_sgx * (pml_x_sgx * ux - dt / rho0 * dpdx)
    else
        rho0_sgx = _stagger_density_1d(rho0)
        @. ux = pml_x_sgx * (pml_x_sgx * ux - dt / rho0_sgx * dpdx)
    end

    # === STEP 3: Velocity sources ===
    if has_velocity_source(source)
        _inject_velocity_source_1d!(ux, source, t_index)
    end

    # === STEP 4: Velocity divergence ===
    duxdx = similar(p)
    spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 1)

    # === STEP 5: Density update with PML ===
    if rho0 isa Real
        @. rhox = pml_x * (pml_x * rhox - dt * rho0 * duxdx)
    else
        @. rhox = pml_x * (pml_x * rhox - dt * rho0 * duxdx)
    end

    # === STEP 6: Pressure sources ===
    if has_pressure_source(source)
        _inject_pressure_source_1d!(rhox, source, medium, t_index)
    end

    # === STEP 7: Equation of state ===
    rho_total = rhox

    _compute_pressure!(p, rho_total, rho_total_prev,
                       scratch1, scratch2, kappa,
                       c0, rho0, medium.BonA, absorb, plans, dt)

    # Update rho_total_prev for dispersion
    if rho_total_prev !== nothing
        rho_total_prev .= rho_total
    end
end

# ============================================================================
# 2D Time stepping
# ============================================================================

function time_step_2d!(
    p::Matrix{Float64},
    ux::Matrix{Float64},
    uy::Matrix{Float64},
    rhox::Matrix{Float64},
    rhoy::Matrix{Float64},
    scratch1::Matrix{ComplexF64},
    scratch2::Matrix{ComplexF64},
    kgrid::KWaveGrid2D,
    medium::KWaveMedium,
    source::KWaveSource,
    pml_x::Vector{Float64},
    pml_y::Vector{Float64},
    pml_x_sgx::Vector{Float64},
    pml_y_sgy::Vector{Float64},
    kappa::Matrix{Float64},
    plans::FFTPlans,
    t_index::Int,
    absorb::Union{Nothing, AbsorptionParams}=nothing,
    rho_total_prev::Union{Nothing, Matrix{Float64}}=nothing,
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density

    # Reshape PML vectors for 2D broadcasting
    pml_x_col = reshape(pml_x, :, 1)
    pml_y_row = reshape(pml_y, 1, :)
    pml_x_sgx_col = reshape(pml_x_sgx, :, 1)
    pml_y_sgy_row = reshape(pml_y_sgy, 1, :)

    # === STEP 1: Pressure gradient via FFT ===
    dpdx = similar(p)
    dpdy = similar(p)
    spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 2)
    spectral_gradient!(dpdy, p, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch1, plans, 2, 2)

    # === STEP 2: Velocity update with PML ===
    if rho0 isa Real
        @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0 * dpdx)
        @. uy = pml_y_sgy_row * (pml_y_sgy_row * uy - dt / rho0 * dpdy)
    else
        rho0_sgx = _stagger_density_2d(rho0, 1)
        rho0_sgy = _stagger_density_2d(rho0, 2)
        @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0_sgx * dpdx)
        @. uy = pml_y_sgy_row * (pml_y_sgy_row * uy - dt / rho0_sgy * dpdy)
    end

    # === STEP 3: Add velocity sources ===
    if has_velocity_source(source)
        _inject_velocity_source_2d!(ux, uy, source, t_index)
    end

    # === STEP 4: Velocity divergence via FFT ===
    duxdx = similar(p)
    duydy = similar(p)
    spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 2)
    spectral_gradient!(duydy, uy, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch1, plans, 2, 2)

    # === STEP 5: Density update with split-field PML ===
    if rho0 isa Real
        @. rhox = pml_x_col * (pml_x_col * rhox - dt * rho0 * duxdx)
        @. rhoy = pml_y_row * (pml_y_row * rhoy - dt * rho0 * duydy)
    else
        @. rhox = pml_x_col * (pml_x_col * rhox - dt * rho0 * duxdx)
        @. rhoy = pml_y_row * (pml_y_row * rhoy - dt * rho0 * duydy)
    end

    # === STEP 6: Add pressure/mass sources ===
    if has_pressure_source(source)
        _inject_pressure_source_2d!(rhox, rhoy, source, medium, t_index)
    end

    # === STEP 7: Equation of state ===
    rho_total = rhox .+ rhoy

    _compute_pressure!(p, rho_total, rho_total_prev,
                      scratch1, scratch2, kappa,
                      c0, rho0, medium.BonA, absorb, plans, dt)

    # Update rho_total_prev for dispersion
    if rho_total_prev !== nothing
        rho_total_prev .= rho_total
    end
end

# ============================================================================
# 3D Time stepping
# ============================================================================

function time_step_3d!(
    p::Array{Float64,3},
    ux::Array{Float64,3},
    uy::Array{Float64,3},
    uz::Array{Float64,3},
    rhox::Array{Float64,3},
    rhoy::Array{Float64,3},
    rhoz::Array{Float64,3},
    scratch1::Array{ComplexF64,3},
    scratch2::Array{ComplexF64,3},
    kgrid::KWaveGrid3D,
    medium::KWaveMedium,
    source::KWaveSource,
    pml_x::Vector{Float64},
    pml_y::Vector{Float64},
    pml_z::Vector{Float64},
    pml_x_sgx::Vector{Float64},
    pml_y_sgy::Vector{Float64},
    pml_z_sgz::Vector{Float64},
    kappa::Array{Float64,3},
    plans::FFTPlans,
    t_index::Int,
    absorb::Union{Nothing, AbsorptionParams}=nothing,
    rho_total_prev::Union{Nothing, Array{Float64,3}}=nothing,
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density
    Nx, Ny, Nz = kgrid.Nx, kgrid.Ny, kgrid.Nz

    # Reshape PML vectors for 3D broadcasting
    pml_x_r = reshape(pml_x, :, 1, 1)
    pml_y_r = reshape(pml_y, 1, :, 1)
    pml_z_r = reshape(pml_z, 1, 1, :)
    pml_x_sgx_r = reshape(pml_x_sgx, :, 1, 1)
    pml_y_sgy_r = reshape(pml_y_sgy, 1, :, 1)
    pml_z_sgz_r = reshape(pml_z_sgz, 1, 1, :)

    # === STEP 1: Pressure gradient via FFT ===
    dpdx = similar(p)
    dpdy = similar(p)
    dpdz = similar(p)
    spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 3)
    spectral_gradient!(dpdy, p, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch1, plans, 2, 3)
    spectral_gradient!(dpdz, p, kgrid.kz_vec, kgrid.ddz_k_shift_pos, scratch1, plans, 3, 3)

    # === STEP 2: Velocity update with PML ===
    if rho0 isa Real
        @. ux = pml_x_sgx_r * (pml_x_sgx_r * ux - dt / rho0 * dpdx)
        @. uy = pml_y_sgy_r * (pml_y_sgy_r * uy - dt / rho0 * dpdy)
        @. uz = pml_z_sgz_r * (pml_z_sgz_r * uz - dt / rho0 * dpdz)
    else
        rho0_sgx = _stagger_density_3d(rho0, 1)
        rho0_sgy = _stagger_density_3d(rho0, 2)
        rho0_sgz = _stagger_density_3d(rho0, 3)
        @. ux = pml_x_sgx_r * (pml_x_sgx_r * ux - dt / rho0_sgx * dpdx)
        @. uy = pml_y_sgy_r * (pml_y_sgy_r * uy - dt / rho0_sgy * dpdy)
        @. uz = pml_z_sgz_r * (pml_z_sgz_r * uz - dt / rho0_sgz * dpdz)
    end

    # === STEP 3: Add velocity sources ===
    if has_velocity_source(source)
        _inject_velocity_source_3d!(ux, uy, uz, source, t_index)
    end

    # === STEP 4: Velocity divergence via FFT ===
    duxdx = similar(p)
    duydy = similar(p)
    duzdz = similar(p)
    spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 3)
    spectral_gradient!(duydy, uy, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch1, plans, 2, 3)
    spectral_gradient!(duzdz, uz, kgrid.kz_vec, kgrid.ddz_k_shift_neg, scratch1, plans, 3, 3)

    # === STEP 5: Density update with split-field PML ===
    if rho0 isa Real
        @. rhox = pml_x_r * (pml_x_r * rhox - dt * rho0 * duxdx)
        @. rhoy = pml_y_r * (pml_y_r * rhoy - dt * rho0 * duydy)
        @. rhoz = pml_z_r * (pml_z_r * rhoz - dt * rho0 * duzdz)
    else
        @. rhox = pml_x_r * (pml_x_r * rhox - dt * rho0 * duxdx)
        @. rhoy = pml_y_r * (pml_y_r * rhoy - dt * rho0 * duydy)
        @. rhoz = pml_z_r * (pml_z_r * rhoz - dt * rho0 * duzdz)
    end

    # === STEP 6: Add pressure/mass sources ===
    if has_pressure_source(source)
        _inject_pressure_source_3d!(rhox, rhoy, rhoz, source, medium, t_index)
    end

    # === STEP 7: Equation of state ===
    rho_total = rhox .+ rhoy .+ rhoz

    _compute_pressure!(p, rho_total, rho_total_prev,
                      scratch1, scratch2, kappa,
                      c0, rho0, medium.BonA, absorb, plans, dt)

    # Update rho_total_prev for dispersion
    if rho_total_prev !== nothing
        rho_total_prev .= rho_total
    end
end

# ============================================================================
# Source injection helpers — 1D
# ============================================================================

function _inject_velocity_source_1d!(ux, source::KWaveSource, t_index::Int)
    if source.u_mask === nothing
        return
    end
    mask_indices = findall(source.u_mask)

    if source.ux !== nothing
        t_col = min(t_index, size(source.ux, 2))
        if source.u_mode == Dirichlet
            for (j, idx) in enumerate(mask_indices)
                ux[idx] = source.ux[j, t_col]
            end
        else
            for (j, idx) in enumerate(mask_indices)
                ux[idx] += source.ux[j, t_col]
            end
        end
    end
end

function _inject_pressure_source_1d!(rhox, source::KWaveSource,
                                     medium::KWaveMedium, t_index::Int)
    if source.p_mask === nothing || source.p === nothing
        return
    end
    mask_indices = findall(source.p_mask)
    t_col = min(t_index, size(source.p, 2))
    c0 = medium.sound_speed

    for (j, idx) in enumerate(mask_indices)
        p_val = source.p[j, t_col]
        c_local = c0 isa Real ? c0 : c0[idx]
        rho_val = p_val / c_local^2
        if source.p_mode == Dirichlet
            rhox[idx] = rho_val
        else
            rhox[idx] += rho_val
        end
    end
end

# ============================================================================
# Source injection helpers — 2D
# ============================================================================

function _inject_velocity_source_2d!(ux, uy, source::KWaveSource, t_index::Int)
    if source.u_mask === nothing
        return
    end
    mask_indices = findall(source.u_mask)

    if source.ux !== nothing
        t_col = min(t_index, size(source.ux, 2))
        if source.u_mode == Dirichlet
            for (j, idx) in enumerate(mask_indices)
                ux[idx] = source.ux[j, t_col]
            end
        else
            for (j, idx) in enumerate(mask_indices)
                ux[idx] += source.ux[j, t_col]
            end
        end
    end

    if source.uy !== nothing
        t_col = min(t_index, size(source.uy, 2))
        if source.u_mode == Dirichlet
            for (j, idx) in enumerate(mask_indices)
                uy[idx] = source.uy[j, t_col]
            end
        else
            for (j, idx) in enumerate(mask_indices)
                uy[idx] += source.uy[j, t_col]
            end
        end
    end
end

function _inject_pressure_source_2d!(rhox, rhoy, source::KWaveSource,
                                     medium::KWaveMedium, t_index::Int)
    if source.p_mask === nothing || source.p === nothing
        return
    end
    mask_indices = findall(source.p_mask)
    t_col = min(t_index, size(source.p, 2))
    c0 = medium.sound_speed

    for (j, idx) in enumerate(mask_indices)
        p_val = source.p[j, t_col]
        c_local = c0 isa Real ? c0 : c0[idx]
        rho_val = p_val / (2 * c_local^2)
        if source.p_mode == Dirichlet
            rhox[idx] = rho_val
            rhoy[idx] = rho_val
        else
            rhox[idx] += rho_val
            rhoy[idx] += rho_val
        end
    end
end

# ============================================================================
# Source injection helpers — 3D
# ============================================================================

function _inject_velocity_source_3d!(ux, uy, uz, source::KWaveSource, t_index::Int)
    if source.u_mask === nothing
        return
    end
    mask_indices = findall(source.u_mask)

    if source.ux !== nothing
        t_col = min(t_index, size(source.ux, 2))
        if source.u_mode == Dirichlet
            for (j, idx) in enumerate(mask_indices)
                ux[idx] = source.ux[j, t_col]
            end
        else
            for (j, idx) in enumerate(mask_indices)
                ux[idx] += source.ux[j, t_col]
            end
        end
    end

    if source.uy !== nothing
        t_col = min(t_index, size(source.uy, 2))
        if source.u_mode == Dirichlet
            for (j, idx) in enumerate(mask_indices)
                uy[idx] = source.uy[j, t_col]
            end
        else
            for (j, idx) in enumerate(mask_indices)
                uy[idx] += source.uy[j, t_col]
            end
        end
    end

    if source.uz !== nothing
        t_col = min(t_index, size(source.uz, 2))
        if source.u_mode == Dirichlet
            for (j, idx) in enumerate(mask_indices)
                uz[idx] = source.uz[j, t_col]
            end
        else
            for (j, idx) in enumerate(mask_indices)
                uz[idx] += source.uz[j, t_col]
            end
        end
    end
end

function _inject_pressure_source_3d!(rhox, rhoy, rhoz, source::KWaveSource,
                                     medium::KWaveMedium, t_index::Int)
    if source.p_mask === nothing || source.p === nothing
        return
    end
    mask_indices = findall(source.p_mask)
    t_col = min(t_index, size(source.p, 2))
    c0 = medium.sound_speed

    for (j, idx) in enumerate(mask_indices)
        p_val = source.p[j, t_col]
        c_local = c0 isa Real ? c0 : c0[idx]
        rho_val = p_val / (3 * c_local^2)
        if source.p_mode == Dirichlet
            rhox[idx] = rho_val
            rhoy[idx] = rho_val
            rhoz[idx] = rho_val
        else
            rhox[idx] += rho_val
            rhoy[idx] += rho_val
            rhoz[idx] += rho_val
        end
    end
end

# ============================================================================
# Sensor recording — generic for all dimensions
# ============================================================================

"""
    record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices, t_index, Nt, dt, rho0)

Record sensor data for the current time step. Works for 1D, 2D, and 3D.

`velocities` is a NamedTuple like `(ux=ux,)`, `(ux=ux, uy=uy)`, or `(ux=ux, uy=uy, uz=uz)`.
"""
function record_sensor_data!(sensor_data::Dict{Symbol, AbstractArray},
                             p::AbstractArray, velocities::NamedTuple,
                             sensor::KWaveSensor, mask_indices::Vector,
                             t_index::Int, Nt::Int, dt::Float64, rho0)
    for field in sensor.record
        if field == :p
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:p][j, t_index] = p[idx]
            end
        elseif field == :p_max
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:p_max][j] = max(sensor_data[:p_max][j], p[idx])
            end
        elseif field == :p_min
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:p_min][j] = min(sensor_data[:p_min][j], p[idx])
            end
        elseif field == :p_rms
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:p_rms][j] += p[idx]^2
            end
        elseif field == :p_final && t_index == Nt
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:p_final][j] = p[idx]
            end
        elseif field == :ux && haskey(velocities, :ux)
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:ux][j, t_index] = velocities.ux[idx]
            end
        elseif field == :uy && haskey(velocities, :uy)
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:uy][j, t_index] = velocities.uy[idx]
            end
        elseif field == :uz && haskey(velocities, :uz)
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:uz][j, t_index] = velocities.uz[idx]
            end
        elseif field == :u_max
            for (j, idx) in enumerate(mask_indices)
                u_mag = _velocity_magnitude(velocities, idx)
                sensor_data[:u_max][j] = max(sensor_data[:u_max][j], u_mag)
            end
        elseif field == :u_rms
            for (j, idx) in enumerate(mask_indices)
                u_mag = _velocity_magnitude(velocities, idx)
                sensor_data[:u_rms][j] += u_mag^2
            end
        elseif field == :u_final && t_index == Nt
            for (j, idx) in enumerate(mask_indices)
                u_mag = _velocity_magnitude(velocities, idx)
                sensor_data[:u_final][j] = u_mag
            end
        elseif field == :I_avg
            for (j, idx) in enumerate(mask_indices)
                # Time-averaged acoustic intensity: <p * u> along each velocity component
                # Simplified: use magnitude of p * u_x for 1D, sum of components for nD
                I_inst = _acoustic_intensity(p[idx], velocities, idx)
                sensor_data[:I_avg][j] += I_inst / Nt
            end
        elseif field == :I_max
            for (j, idx) in enumerate(mask_indices)
                I_inst = _acoustic_intensity(p[idx], velocities, idx)
                sensor_data[:I_max][j] = max(sensor_data[:I_max][j], I_inst)
            end
        end
    end
end

function _velocity_magnitude(velocities::NamedTuple, idx)
    mag_sq = 0.0
    for v in values(velocities)
        mag_sq += v[idx]^2
    end
    return sqrt(mag_sq)
end

function _acoustic_intensity(p_val::Float64, velocities::NamedTuple, idx)
    # Acoustic intensity magnitude: |p * u|
    I_sq = 0.0
    for v in values(velocities)
        I_sq += (p_val * v[idx])^2
    end
    return sqrt(I_sq)
end

"""
    finalize_sensor_data!(sensor_data, Nt)

Post-process sensor data after the time loop.
"""
function finalize_sensor_data!(sensor_data::Dict{Symbol, AbstractArray}, Nt::Int)
    if haskey(sensor_data, :p_rms)
        @. sensor_data[:p_rms] = sqrt(sensor_data[:p_rms] / Nt)
    end
    if haskey(sensor_data, :u_rms)
        @. sensor_data[:u_rms] = sqrt(sensor_data[:u_rms] / Nt)
    end
end

# ============================================================================
# Sensor directivity (frequency-domain filtering)
# ============================================================================

"""
    apply_sensor_directivity!(sensor_data, sensor, kgrid)

Apply sensor directivity to recorded pressure data using frequency-domain filtering.
"""
function apply_sensor_directivity!(sensor_data::Dict{Symbol, AbstractArray},
                                   sensor::KWaveSensor, kgrid::AbstractKWaveGrid)
    if sensor.directivity_angle === nothing || sensor.directivity_size === nothing
        return
    end
    if !haskey(sensor_data, :p)
        return
    end

    # Apply directivity as a frequency-domain angular filter
    # Directivity pattern: sinc(k * d * sin(theta) / 2)
    # where d is the element size and theta is the directivity angle
    p_data = sensor_data[:p]
    n_sensor = size(p_data, 1)
    Nt = size(p_data, 2)
    dt = kgrid.dt[]
    Fs = 1.0 / dt

    for j in 1:n_sensor
        angle = sensor.directivity_angle isa Real ? sensor.directivity_angle : sensor.directivity_angle[j]
        if angle ≈ 0
            continue
        end

        signal = p_data[j, :]
        freq_content = fft(signal)
        freqs = fftfreq(Nt, Fs)

        for f_idx in 1:Nt
            k_val = 2π * abs(freqs[f_idx]) / (kgrid isa KWaveGrid1D ? kgrid.dx : minimum(grid_spacing(kgrid)))
            arg = k_val * sensor.directivity_size * sin(angle) / 2
            directivity = arg ≈ 0 ? 1.0 : sin(arg) / arg
            freq_content[f_idx] *= directivity
        end

        p_data[j, :] .= real.(ifft(freq_content))
    end
end

"""
    apply_frequency_response!(sensor_data, sensor, kgrid)

Apply sensor frequency response (bandpass filter) to recorded data.
"""
function apply_frequency_response!(sensor_data::Dict{Symbol, AbstractArray},
                                    sensor::KWaveSensor, kgrid::AbstractKWaveGrid)
    if sensor.frequency_response === nothing
        return
    end
    if !haskey(sensor_data, :p)
        return
    end

    center_freq, bandwidth = sensor.frequency_response
    p_data = sensor_data[:p]
    n_sensor = size(p_data, 1)
    Nt = size(p_data, 2)
    dt = kgrid.dt[]
    Fs = 1.0 / dt

    for j in 1:n_sensor
        signal = p_data[j, :]
        freq_content = fft(signal)
        freqs = fftfreq(Nt, Fs)

        for f_idx in 1:Nt
            f = abs(freqs[f_idx])
            # Gaussian frequency response
            response = exp(-((f - center_freq)^2) / (2 * bandwidth^2))
            freq_content[f_idx] *= response
        end

        p_data[j, :] .= real.(ifft(freq_content))
    end
end
