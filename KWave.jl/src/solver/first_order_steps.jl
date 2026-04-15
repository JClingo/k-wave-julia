# ============================================================================
# KWave.jl — Time-stepping operations for the k-space first-order solver
# ============================================================================

# ============================================================================
# Absorption parameters struct
# ============================================================================

"""
    AbsorptionParams{T,A}

Pre-computed absorption and dispersion parameters for the equation of state.
Following the k-Wave formulation (Treeby & Cox, 2010):
- `absorb_tau`: scalar absorption prefactor
- `absorb_eta`: scalar dispersion prefactor
- `absorb_nabla1`: k-space operator k^y  (may be full-grid or rfft-shaped)
- `absorb_nabla2`: k-space operator k^(y-1)

`T` is the floating-point working precision; `A` is the array type of the
nabla operators (either full-grid or rfft-truncated).
"""
struct AbsorptionParams{T<:AbstractFloat, A<:AbstractArray{T}}
    absorb_tau::T                 # -2 * alpha_nepers * c_ref^(y-1)
    absorb_eta::T                 # 2 * alpha_nepers * c_ref^(y-1) * tan(π*y/2)
    absorb_nabla1::A              # k^y
    absorb_nabla2::A              # k^(y-1)
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

function initialize_p0_1d!(p::AbstractVector, ux::AbstractVector,
                           rhox::AbstractVector,
                           source::KWaveSource, kgrid::KWaveGrid1D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::AbstractVector{<:Complex},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    FT = eltype(p)
    p0 = FT.(source.p0)

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

function initialize_p0_2d!(p::AbstractMatrix, ux::AbstractMatrix, uy::AbstractMatrix,
                           rhox::AbstractMatrix, rhoy::AbstractMatrix,
                           source::KWaveSource, kgrid::KWaveGrid2D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::AbstractMatrix{<:Complex},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    FT = eltype(p)
    p0 = FT.(source.p0)

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

function initialize_p0_3d!(p::AbstractArray{<:Real,3}, ux::AbstractArray{<:Real,3},
                           uy::AbstractArray{<:Real,3}, uz::AbstractArray{<:Real,3},
                           rhox::AbstractArray{<:Real,3}, rhoy::AbstractArray{<:Real,3},
                           rhoz::AbstractArray{<:Real,3},
                           source::KWaveSource, kgrid::KWaveGrid3D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::AbstractArray{<:Complex,3},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    FT = eltype(p)
    p0 = FT.(source.p0)

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

Uses rfft/irfft (real-to-complex): `scratch1` and `scratch2` are complex
arrays of shape `plans.rfft_dims` (first dim = Nx÷2+1), and `kappa_r` /
`absorb.absorb_nabla1/2` must already be rfft-truncated to that shape.

For the absorbing case, uses the k-Wave formulation (Treeby & Cox, 2010):
  p_hat = c₀² · κ · (1 + absorb_nabla1) · rfft(ρ)
For dispersion, additionally subtracts:
  c₀² · κ · absorb_nabla2 · rfft((ρ - ρ_prev) / dt)

`tmp_real` is a pre-allocated real scratch of the same shape as `rho_total`
(full grid size). It is used for staging heterogeneous and dispersion
computations without allocating inside the hot loop.
"""
function _compute_pressure!(
    p::AbstractArray{T},
    rho_total::AbstractArray{T},
    rho_total_prev::Union{Nothing, AbstractArray{T}},
    scratch1::AbstractArray{Complex{T}},   # rfft-shaped: (Nx÷2+1, …)
    scratch2::AbstractArray{Complex{T}},   # rfft-shaped: (Nx÷2+1, …)
    tmp_real::AbstractArray{T},            # full grid shape, real scratch
    kappa_r::AbstractArray{T},             # rfft-shaped kappa
    c0,     # scalar or array
    rho0,   # scalar or array (density for nonlinear term)
    BonA,   # nothing or scalar or array
    absorb::Union{Nothing, AbsorptionParams},
    plans::FFTPlans,
    dt::Real,
) where T <: AbstractFloat
    is_absorbing = absorb !== nothing && absorb.mode != :no_absorption
    is_nonlinear = BonA !== nothing

    if c0 isa Real
        c_sq = T(c0)^2

        # rfft of effective density → scratch1
        if is_nonlinear
            @. tmp_real = rho_total * (1 + BonA / 2 * rho_total / rho0)
            mul!(scratch1, plans.forward, tmp_real)
        else
            mul!(scratch1, plans.forward, rho_total)
        end

        if !is_absorbing
            # Lossless: p = irfft( c0² · κ · rfft(ρ) )
            @. scratch1 = c_sq * kappa_r * scratch1

        elseif absorb.mode == :no_dispersion || rho_total_prev === nothing
            # Absorption only (no dispersion):
            # p = irfft( c0² · κ · (1 + nabla1) · rfft(ρ) )
            @. scratch1 = c_sq * kappa_r * (1 + absorb.absorb_nabla1) * scratch1

        else
            # Full absorption + dispersion:
            # p = irfft( c0² · κ · ((1+nabla1)·rfft(ρ) − nabla2·rfft((ρ−ρ_prev)/dt)) )
            if is_nonlinear
                # tmp_real currently has rho_eff; compute (rho_eff - rho_eff_prev)/dt
                # rho_eff_prev allocation is unavoidable in the nonlinear+dispersion case
                rho_eff_prev = @. rho_total_prev * (1 + BonA / 2 * rho_total_prev / rho0)
                @. tmp_real = (tmp_real - rho_eff_prev) / dt
            else
                @. tmp_real = (rho_total - rho_total_prev) / dt
            end
            mul!(scratch2, plans.forward, tmp_real)
            @. scratch1 = c_sq * kappa_r * (
                (1 + absorb.absorb_nabla1) * scratch1
                - absorb.absorb_nabla2 * scratch2
            )
        end

        # irfft writes directly to p (real output, full grid shape)
        mul!(p, plans.inverse, scratch1)

    else
        # Heterogeneous media: fuse c² · ρ_eff into tmp_real, then rfft
        if is_nonlinear
            @. tmp_real = c0^2 * rho_total * (1 + BonA / 2 * rho_total / rho0)
        else
            @. tmp_real = c0^2 * rho_total
        end
        mul!(scratch1, plans.forward, tmp_real)

        if !is_absorbing
            @. scratch1 = kappa_r * scratch1
        else
            @. scratch1 = kappa_r * (1 + absorb.absorb_nabla1) * scratch1
        end

        mul!(p, plans.inverse, scratch1)
    end
end

# ============================================================================
# 1D Time stepping
# ============================================================================

function time_step_1d!(
    p::AbstractVector,
    ux::AbstractVector,
    rhox::AbstractVector,
    scratch1::AbstractVector{<:Complex},
    scratch2::AbstractVector{<:Complex},
    dpdx::AbstractVector,                               # pre-allocated pressure gradient
    duxdx::AbstractVector,                              # pre-allocated velocity divergence
    kgrid::KWaveGrid1D,
    medium::KWaveMedium,
    source::KWaveSource,
    pml_x::AbstractVector,
    pml_x_sgx::AbstractVector,
    kappa_r::AbstractVector,                            # rfft-shaped kappa
    plans::FFTPlans,
    t_index::Int,
    absorb::Union{Nothing, AbsorptionParams},
    tmp_real::AbstractVector,                           # full-grid real scratch
    rho_total_prev::Union{Nothing, AbstractVector}=nothing,
    rho0_sgx::Union{Nothing, AbstractVector}=nothing,  # pre-computed staggered density
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density

    # === STEP 1: Pressure gradient ===
    spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 1)

    # === STEP 2: Velocity update with PML ===
    if rho0_sgx !== nothing
        @. ux = pml_x_sgx * (pml_x_sgx * ux - dt / rho0_sgx * dpdx)
    else
        @. ux = pml_x_sgx * (pml_x_sgx * ux - dt / rho0 * dpdx)
    end

    # === STEP 3: Velocity sources ===
    if has_velocity_source(source)
        _inject_velocity_source_1d!(ux, source, t_index)
    end

    # === STEP 4: Velocity divergence ===
    spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 1)

    # === STEP 5: Density update with PML ===
    @. rhox = pml_x * (pml_x * rhox - dt * rho0 * duxdx)

    # === STEP 6: Pressure sources ===
    if has_pressure_source(source)
        _inject_pressure_source_1d!(rhox, source, medium, t_index)
    end

    # === STEP 7: Equation of state ===
    rho_total = rhox

    _compute_pressure!(p, rho_total, rho_total_prev,
                       scratch1, scratch2, tmp_real, kappa_r,
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
    p::AbstractMatrix,
    ux::AbstractMatrix,
    uy::AbstractMatrix,
    rhox::AbstractMatrix,
    rhoy::AbstractMatrix,
    scratch1::AbstractMatrix{<:Complex},
    scratch2::AbstractMatrix{<:Complex},
    dpdx::AbstractMatrix,                               # pre-allocated pressure gradient x
    dpdy::AbstractMatrix,                               # pre-allocated pressure gradient y
    duxdx::AbstractMatrix,                              # pre-allocated velocity divergence x
    duydy::AbstractMatrix,                              # pre-allocated velocity divergence y
    rho_total::AbstractMatrix,                          # pre-allocated density sum buffer
    kgrid::KWaveGrid2D,
    medium::KWaveMedium,
    source::KWaveSource,
    pml_x_col,                                          # pre-reshaped PML (Nx×1)
    pml_y_row,                                          # pre-reshaped PML (1×Ny)
    pml_x_sgx_col,                                      # pre-reshaped staggered PML (Nx×1)
    pml_y_sgy_row,                                      # pre-reshaped staggered PML (1×Ny)
    kappa_r::AbstractMatrix,                            # rfft-shaped kappa
    plans::FFTPlans,
    t_index::Int,
    absorb::Union{Nothing, AbsorptionParams}=nothing,
    tmp_real::Union{Nothing, AbstractMatrix}=nothing,   # full-grid real scratch
    rho_total_prev::Union{Nothing, AbstractMatrix}=nothing,
    rho0_sgx::Union{Nothing, AbstractMatrix}=nothing,  # pre-computed staggered density x
    rho0_sgy::Union{Nothing, AbstractMatrix}=nothing,  # pre-computed staggered density y
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density

    # === STEP 1: Pressure gradient via FFT ===
    spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 2)
    spectral_gradient!(dpdy, p, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch1, plans, 2, 2)

    # === STEP 2: Velocity update with PML ===
    if rho0_sgx !== nothing
        @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0_sgx * dpdx)
        @. uy = pml_y_sgy_row * (pml_y_sgy_row * uy - dt / rho0_sgy * dpdy)
    else
        @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0 * dpdx)
        @. uy = pml_y_sgy_row * (pml_y_sgy_row * uy - dt / rho0 * dpdy)
    end

    # === STEP 3: Add velocity sources ===
    if has_velocity_source(source)
        _inject_velocity_source_2d!(ux, uy, source, t_index)
    end

    # === STEP 4: Velocity divergence via FFT ===
    spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 2)
    spectral_gradient!(duydy, uy, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch1, plans, 2, 2)

    # === STEP 5: Density update with split-field PML ===
    @. rhox = pml_x_col * (pml_x_col * rhox - dt * rho0 * duxdx)
    @. rhoy = pml_y_row * (pml_y_row * rhoy - dt * rho0 * duydy)

    # === STEP 6: Add pressure/mass sources ===
    if has_pressure_source(source)
        _inject_pressure_source_2d!(rhox, rhoy, source, medium, t_index)
    end

    # === STEP 7: Equation of state ===
    @. rho_total = rhox + rhoy

    _compute_pressure!(p, rho_total, rho_total_prev,
                      scratch1, scratch2, tmp_real, kappa_r,
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
    p::AbstractArray{<:Real,3},
    ux::AbstractArray{<:Real,3},
    uy::AbstractArray{<:Real,3},
    uz::AbstractArray{<:Real,3},
    rhox::AbstractArray{<:Real,3},
    rhoy::AbstractArray{<:Real,3},
    rhoz::AbstractArray{<:Real,3},
    scratch1::AbstractArray{<:Complex,3},
    scratch2::AbstractArray{<:Complex,3},
    dpdx::AbstractArray{<:Real,3},                     # pre-allocated pressure gradient x
    dpdy::AbstractArray{<:Real,3},                     # pre-allocated pressure gradient y
    dpdz::AbstractArray{<:Real,3},                     # pre-allocated pressure gradient z
    duxdx::AbstractArray{<:Real,3},                    # pre-allocated velocity divergence x
    duydy::AbstractArray{<:Real,3},                    # pre-allocated velocity divergence y
    duzdz::AbstractArray{<:Real,3},                    # pre-allocated velocity divergence z
    rho_total::AbstractArray{<:Real,3},                # pre-allocated density sum buffer
    kgrid::KWaveGrid3D,
    medium::KWaveMedium,
    source::KWaveSource,
    pml_x_r,                                           # pre-reshaped PML (Nx×1×1)
    pml_y_r,                                           # pre-reshaped PML (1×Ny×1)
    pml_z_r,                                           # pre-reshaped PML (1×1×Nz)
    pml_x_sgx_r,                                       # pre-reshaped staggered PML (Nx×1×1)
    pml_y_sgy_r,                                       # pre-reshaped staggered PML (1×Ny×1)
    pml_z_sgz_r,                                       # pre-reshaped staggered PML (1×1×Nz)
    kappa_r::AbstractArray{<:Real,3},                  # rfft-shaped kappa
    plans::FFTPlans,
    t_index::Int,
    absorb::Union{Nothing, AbsorptionParams}=nothing,
    tmp_real::Union{Nothing, AbstractArray{<:Real,3}}=nothing,  # full-grid real scratch
    rho_total_prev::Union{Nothing, AbstractArray{<:Real,3}}=nothing,
    rho0_sgx::Union{Nothing, AbstractArray{<:Real,3}}=nothing, # pre-computed staggered density x
    rho0_sgy::Union{Nothing, AbstractArray{<:Real,3}}=nothing, # pre-computed staggered density y
    rho0_sgz::Union{Nothing, AbstractArray{<:Real,3}}=nothing, # pre-computed staggered density z
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density

    # === STEP 1: Pressure gradient via FFT ===
    spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 3)
    spectral_gradient!(dpdy, p, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch1, plans, 2, 3)
    spectral_gradient!(dpdz, p, kgrid.kz_vec, kgrid.ddz_k_shift_pos, scratch1, plans, 3, 3)

    # === STEP 2: Velocity update with PML ===
    if rho0_sgx !== nothing
        @. ux = pml_x_sgx_r * (pml_x_sgx_r * ux - dt / rho0_sgx * dpdx)
        @. uy = pml_y_sgy_r * (pml_y_sgy_r * uy - dt / rho0_sgy * dpdy)
        @. uz = pml_z_sgz_r * (pml_z_sgz_r * uz - dt / rho0_sgz * dpdz)
    else
        @. ux = pml_x_sgx_r * (pml_x_sgx_r * ux - dt / rho0 * dpdx)
        @. uy = pml_y_sgy_r * (pml_y_sgy_r * uy - dt / rho0 * dpdy)
        @. uz = pml_z_sgz_r * (pml_z_sgz_r * uz - dt / rho0 * dpdz)
    end

    # === STEP 3: Add velocity sources ===
    if has_velocity_source(source)
        _inject_velocity_source_3d!(ux, uy, uz, source, t_index)
    end

    # === STEP 4: Velocity divergence via FFT ===
    spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 3)
    spectral_gradient!(duydy, uy, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch1, plans, 2, 3)
    spectral_gradient!(duzdz, uz, kgrid.kz_vec, kgrid.ddz_k_shift_neg, scratch1, plans, 3, 3)

    # === STEP 5: Density update with split-field PML ===
    @. rhox = pml_x_r * (pml_x_r * rhox - dt * rho0 * duxdx)
    @. rhoy = pml_y_r * (pml_y_r * rhoy - dt * rho0 * duydy)
    @. rhoz = pml_z_r * (pml_z_r * rhoz - dt * rho0 * duzdz)

    # === STEP 6: Add pressure/mass sources ===
    if has_pressure_source(source)
        _inject_pressure_source_3d!(rhox, rhoy, rhoz, source, medium, t_index)
    end

    # === STEP 7: Equation of state ===
    @. rho_total = rhox + rhoy + rhoz

    _compute_pressure!(p, rho_total, rho_total_prev,
                      scratch1, scratch2, tmp_real, kappa_r,
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
