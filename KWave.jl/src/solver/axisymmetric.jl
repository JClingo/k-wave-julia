# ============================================================================
# KWave.jl — Axisymmetric solver (kspace_first_order_as)
# ============================================================================

"""
    kspace_first_order_as(kgrid, medium, source, sensor; kwargs...)

Run an axisymmetric k-space pseudospectral time-domain simulation.

Uses a 2D grid where the first dimension (x) is the axial direction and
the second dimension (y) is the radial direction. The solution is axially
symmetric about the x-axis, enabling efficient 3D-equivalent simulations
on a 2D grid.

The axisymmetric formulation uses a modified wave equation that accounts
for the 1/r geometric spreading term in cylindrical coordinates.

# Arguments
Same as `kspace_first_order` for 2D, with the convention that:
- Dimension 1 (x) = axial direction
- Dimension 2 (y) = radial direction (r ≥ 0)

# Returns
`SimulationOutput` containing recorded sensor data.
"""
function kspace_first_order_as(
    kgrid::KWaveGrid2D,
    medium::KWaveMedium,
    source::KWaveSource,
    sensor::KWaveSensor;
    pml_inside::Bool=true,
    pml_size::Union{Int, Tuple{Int, Int}}=20,
    pml_alpha::Union{Float64, Tuple{Float64, Float64}}=2.0,
    smooth_p0::Bool=true,
    smooth_c0::Bool=false,
    smooth_rho0::Bool=false,
    data_cast::Type{T}=Float64,
    plot_sim::Bool=false,
    plot_layout::Symbol=:default,
    plot_scale::Union{Symbol, Tuple}=:auto,
    record_movie::Union{Nothing, String}=nothing,
    progress_callback::Union{Nothing, Function}=nothing,
) where T <: AbstractFloat

    # ================================================================
    # 1. Input validation
    # ================================================================
    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running.")
    end
    validate_record_fields(sensor)

    Nx, Nr = kgrid.Nx, kgrid.Ny  # x = axial, y = radial
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    pml_x_size, pml_r_size = pml_size isa Tuple ? pml_size : (pml_size, pml_size)
    pml_x_alpha, pml_r_alpha = pml_alpha isa Tuple ? pml_alpha : (pml_alpha, pml_alpha)

    # ================================================================
    # 2. Pre-compute PML
    # ================================================================
    pml_x = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt)
    pml_r = get_pml(Nr, kgrid.dy, pml_r_size, pml_r_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt; staggered=true)
    pml_r_sgr = get_pml(Nr, kgrid.dy, pml_r_size, pml_r_alpha, dt; staggered=true)

    # ================================================================
    # 3. Pre-compute k-space correction and absorption
    # ================================================================
    c_ref = medium.sound_speed isa Real ? medium.sound_speed : mean(medium.sound_speed)
    absorb = _precompute_absorption(medium, kgrid.k, c_ref)

    # ================================================================
    # 4. Create FFT plans
    # ================================================================
    plans = create_fft_plans(kgrid; data_cast=T)

    # rfft-shaped kappa and absorption operators: slice first dim to Nx÷2+1
    n1 = plans.rfft_dims[1]
    kappa = T.(_compute_kappa(kgrid.k, c_ref, dt))[1:n1, :]
    absorb = if absorb !== nothing
        AbsorptionParams(
            T(absorb.absorb_tau),
            T(absorb.absorb_eta),
            T.(absorb.absorb_nabla1)[1:n1, :],
            T.(absorb.absorb_nabla2)[1:n1, :],
            absorb.mode,
        )
    else
        nothing
    end

    # ================================================================
    # 5. Radial coordinate array for geometric correction
    # ================================================================
    # r_vec: radial distance from axis (y=0 maps to r=0 at grid center)
    # For axisymmetric, the radial axis goes from r=0 to r_max
    # Use absolute value of y coordinates as radial distance
    r_vec = abs.(kgrid.y_vec)
    # Avoid division by zero at r=0
    r_vec_safe = copy(r_vec)
    r_vec_safe[r_vec .== 0] .= kgrid.dy / 2  # half grid point

    # ================================================================
    # 6. Allocate field arrays
    # ================================================================
    p = zeros(Float64, Nx, Nr)
    ux = zeros(Float64, Nx, Nr)   # axial velocity
    ur = zeros(Float64, Nx, Nr)   # radial velocity
    rhox = zeros(Float64, Nx, Nr)
    rhor = zeros(Float64, Nx, Nr)
    scratch1 = zeros(ComplexF64, plans.rfft_dims...)
    scratch2 = zeros(ComplexF64, plans.rfft_dims...)
    tmp_real  = zeros(T, Nx, Nr)   # real staging buffer for _compute_pressure!

    rho_total_prev = (absorb !== nothing && absorb.mode != :no_dispersion) ?
        zeros(Float64, Nx, Nr) : nothing

    # ================================================================
    # 7. Sensor data
    # ================================================================
    sensor_data = _create_sensor_data(sensor, (Nx, Nr), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{2}[]
    end

    # ================================================================
    # 8. Initialize from p0
    # ================================================================
    initialize_p0_2d!(p, ux, ur, rhox, rhor, source, kgrid, medium, plans, scratch1, smooth_p0)

    # ================================================================
    # 9. Time loop
    # ================================================================
    prog = Progress(Nt; desc="k-Wave AS: ", enabled=true)

    pml_x_col = reshape(pml_x, :, 1)
    pml_r_row = reshape(pml_r, 1, :)
    pml_x_sgx_col = reshape(pml_x_sgx, :, 1)
    pml_r_sgr_row = reshape(pml_r_sgr, 1, :)
    r_row = reshape(r_vec_safe, 1, :)

    c0 = medium.sound_speed
    rho0 = medium.density

    for t_index in 1:Nt
        # === STEP 1: Pressure gradient ===
        dpdx = similar(p)
        dpdr = similar(p)
        spectral_gradient!(dpdx, p, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch1, plans, 1, 2)
        spectral_gradient!(dpdr, p, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch1, plans, 2, 2)

        # === STEP 2: Velocity update with PML ===
        if rho0 isa Real
            @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0 * dpdx)
            @. ur = pml_r_sgr_row * (pml_r_sgr_row * ur - dt / rho0 * dpdr)
        else
            rho0_sgx = _stagger_density_2d(rho0, 1)
            rho0_sgr = _stagger_density_2d(rho0, 2)
            @. ux = pml_x_sgx_col * (pml_x_sgx_col * ux - dt / rho0_sgx * dpdx)
            @. ur = pml_r_sgr_row * (pml_r_sgr_row * ur - dt / rho0_sgr * dpdr)
        end

        # === STEP 3: Velocity sources ===
        if has_velocity_source(source)
            _inject_velocity_source_2d!(ux, ur, source, t_index)
        end

        # === STEP 4: Velocity divergence with axisymmetric correction ===
        # div(u) = dux/dx + dur/dr + ur/r  (cylindrical coordinates)
        duxdx = similar(p)
        durdr = similar(p)
        spectral_gradient!(duxdx, ux, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch1, plans, 1, 2)
        spectral_gradient!(durdr, ur, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch1, plans, 2, 2)

        # === STEP 5: Density update with axisymmetric term ===
        if rho0 isa Real
            @. rhox = pml_x_col * (pml_x_col * rhox - dt * rho0 * duxdx)
            # Radial component includes the 1/r geometric term
            @. rhor = pml_r_row * (pml_r_row * rhor - dt * rho0 * (durdr + ur / r_row))
        else
            @. rhox = pml_x_col * (pml_x_col * rhox - dt * rho0 * duxdx)
            @. rhor = pml_r_row * (pml_r_row * rhor - dt * rho0 * (durdr + ur / r_row))
        end

        # === STEP 6: Pressure sources ===
        if has_pressure_source(source)
            _inject_pressure_source_2d!(rhox, rhor, source, medium, t_index)
        end

        # === STEP 7: Equation of state ===
        rho_total = rhox .+ rhor

        _compute_pressure!(p, rho_total, rho_total_prev,
                          scratch1, scratch2, tmp_real, kappa,
                          c0, rho0, medium.BonA, absorb, plans, dt)

        if rho_total_prev !== nothing
            rho_total_prev .= rho_total
        end

        # === Record sensor data ===
        if !isempty(mask_indices)
            velocities = (ux=ux, uy=ur)
            record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices,
                               t_index, Nt, dt, medium.density)
        end

        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    # ================================================================
    # 10. Finalize
    # ================================================================
    finalize_sensor_data!(sensor_data, Nt)
    apply_sensor_directivity!(sensor_data, sensor, kgrid)
    apply_frequency_response!(sensor_data, sensor, kgrid)

    return SimulationOutput(sensor_data)
end
