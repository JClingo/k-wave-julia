# ============================================================================
# KWave.jl — Time-stepping operations for the k-space first-order solver
# ============================================================================

"""
    initialize_p0_2d!(p, ux, uy, rhox, rhoy, source, kgrid, medium, plans, scratch, smooth_p0)

Initialize fields from initial pressure distribution (photoacoustic source).

For an initial value problem, the initial pressure is used to set:
1. The pressure field directly
2. Split density fields (rhox, rhoy) from p = c² * (rhox + rhoy)
3. Initial velocity fields from the pressure gradient
   (offset by dt/2 for leap-frog staggering)
"""
function initialize_p0_2d!(p::Matrix{Float64}, ux::Matrix{Float64}, uy::Matrix{Float64},
                           rhox::Matrix{Float64}, rhoy::Matrix{Float64},
                           source::KWaveSource, kgrid::KWaveGrid2D, medium::KWaveMedium,
                           plans::FFTPlans, scratch::Matrix{ComplexF64},
                           smooth_p0::Bool)
    if !has_p0(source)
        return
    end

    p0 = Float64.(source.p0)

    # Smooth initial pressure if requested (reduces Gibbs ringing)
    if smooth_p0 && kgrid.Nx > 1 && kgrid.Ny > 1
        p0 = smooth(p0; restore_max=true)
    end

    # Set initial pressure
    p .= p0

    # Compute sound speed squared (scalar or array)
    c0_sq = if medium.sound_speed isa Real
        medium.sound_speed^2
    else
        medium.sound_speed.^2
    end

    # Split density equally between x and y components
    # From p = c² * (rhox + rhoy), with equal split:
    @. rhox = p0 / (2 * c0_sq)
    @. rhoy = p0 / (2 * c0_sq)

    # Compute initial particle velocity from pressure gradient
    # ux = -dt / (2 * rho0) * dp0/dx
    # The factor of 2 accounts for the half time-step offset in leap-frog
    dt = kgrid.dt[]
    rho0 = medium.density

    # dp0/dx using spectral gradient on staggered grid
    dpdx = similar(p)
    dpdy = similar(p)
    spectral_gradient!(dpdx, p0, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 2)
    spectral_gradient!(dpdy, p0, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 2)

    if rho0 isa Real
        @. ux = -dt / (2 * rho0) * dpdx
        @. uy = -dt / (2 * rho0) * dpdy
    else
        @. ux = -dt / (2 * rho0) * dpdx
        @. uy = -dt / (2 * rho0) * dpdy
    end
end

"""
    time_step_2d!(p, ux, uy, rhox, rhoy, scratch1, scratch2,
                  kgrid, medium, source, pml_x, pml_y, pml_x_sgx, pml_y_sgy,
                  kappa, plans, t_index)

Execute one time step of the 2D k-space first-order solver.

Algorithm:
1. Compute pressure gradient via FFT (dp/dx, dp/dy on staggered grid)
2. Update velocity fields with PML absorption
3. Inject velocity sources (if any)
4. Compute velocity divergence via FFT (dux/dx, duy/dy)
5. Update split density fields with PML absorption
6. Inject pressure/mass sources (if any)
7. Compute pressure from equation of state
"""
function time_step_2d!(
    # Field arrays (mutated each step)
    p::Matrix{Float64},
    ux::Matrix{Float64},
    uy::Matrix{Float64},
    rhox::Matrix{Float64},
    rhoy::Matrix{Float64},
    # Complex scratch arrays
    scratch1::Matrix{ComplexF64},
    scratch2::Matrix{ComplexF64},
    # Grid and physics
    kgrid::KWaveGrid2D,
    medium::KWaveMedium,
    source::KWaveSource,
    # PML arrays (1D, broadcast)
    pml_x::Vector{Float64},
    pml_y::Vector{Float64},
    pml_x_sgx::Vector{Float64},
    pml_y_sgy::Vector{Float64},
    # k-space correction
    kappa::Matrix{Float64},
    # FFT plans
    plans::FFTPlans,
    # Current time index
    t_index::Int,
)
    dt = kgrid.dt[]
    c0 = medium.sound_speed
    rho0 = medium.density

    # Reshape PML vectors for 2D broadcasting
    pml_x_col = reshape(pml_x, :, 1)         # (Nx, 1)
    pml_y_row = reshape(pml_y, 1, :)         # (1, Ny)
    pml_x_sgx_col = reshape(pml_x_sgx, :, 1) # (Nx, 1)
    pml_y_sgy_row = reshape(pml_y_sgy, 1, :) # (1, Ny)

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
        # Heterogeneous density: need staggered grid interpolation
        # For now, use direct density (TODO: proper staggered grid interpolation)
        @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0 * dpdx)
        @. uy = pml_y_sgy_row * (pml_y_sgy_row * uy - dt / rho0 * dpdy)
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

    # === STEP 7: Pressure from equation of state ===
    if is_lossless(medium)
        # Lossless case: p = c0^2 * (rhox + rhoy)
        # Apply k-space correction in frequency domain
        rho_total = rhox .+ rhoy
        scratch1 .= complex.(rho_total)
        plans.forward * scratch1

        if c0 isa Real
            @. scratch1 = c0^2 * kappa * scratch1
        else
            c0_sq = c0.^2
            scratch2 .= complex.(c0_sq)
            # For heterogeneous media: FFT(c^2 * rho) approximation
            # In the simplest case, just use pointwise multiplication
            scratch1 .= complex.(c0_sq .* rho_total)
            plans.forward * scratch1
            @. scratch1 = kappa * scratch1
        end

        plans.inverse * scratch1
        @. p = real(scratch1)
    else
        # Absorbing case: TODO in Phase 2
        # For now, fall back to lossless
        if c0 isa Real
            @. p = c0^2 * (rhox + rhoy)
        else
            @. p = c0^2 * (rhox + rhoy)
        end
    end
end

# ============================================================================
# Source injection helpers
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
        # Convert pressure to density and split equally
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
# Sensor recording
# ============================================================================

"""
    record_sensor_data!(sensor_data, p, ux, uy, sensor, mask_indices, t_index, Nt)

Record sensor data for the current time step.
"""
function record_sensor_data!(sensor_data::Dict{Symbol, AbstractArray},
                             p::Matrix{Float64}, ux::Matrix{Float64}, uy::Matrix{Float64},
                             sensor::KWaveSensor, mask_indices::Vector{CartesianIndex{2}},
                             t_index::Int, Nt::Int)
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
        elseif field == :ux
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:ux][j, t_index] = ux[idx]
            end
        elseif field == :uy
            for (j, idx) in enumerate(mask_indices)
                sensor_data[:uy][j, t_index] = uy[idx]
            end
        end
    end
end

"""
    finalize_sensor_data!(sensor_data, Nt)

Post-process sensor data after the time loop (e.g., finalize RMS values).
"""
function finalize_sensor_data!(sensor_data::Dict{Symbol, AbstractArray}, Nt::Int)
    if haskey(sensor_data, :p_rms)
        @. sensor_data[:p_rms] = sqrt(sensor_data[:p_rms] / Nt)
    end
    if haskey(sensor_data, :u_rms)
        @. sensor_data[:u_rms] = sqrt(sensor_data[:u_rms] / Nt)
    end
end
