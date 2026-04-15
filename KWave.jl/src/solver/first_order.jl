# ============================================================================
# KWave.jl — kspace_first_order: Main simulation entry points (1D, 2D, 3D)
# ============================================================================

# ============================================================================
# Shared helpers
# ============================================================================

"""
Pre-compute absorption parameters from medium properties.
Returns `nothing` if medium is lossless, otherwise an `AbsorptionParams`.
"""
function _precompute_absorption(medium::KWaveMedium, k_grid::AbstractArray, c_ref::Float64)
    if is_lossless(medium)
        return nothing
    end

    y = medium.alpha_power
    alpha_nepers = db2neper(medium.alpha_coeff, y)

    # Absorption and dispersion prefactors (Treeby & Cox, 2010)
    absorb_tau = -2 * alpha_nepers * c_ref^(y - 1)
    absorb_eta = 2 * alpha_nepers * c_ref^(y - 1) * tan(π * y / 2)

    # Pre-multiply into k-space operators (used directly in frequency domain)
    absorb_nabla1 = absorb_tau .* k_grid.^y          # absorption: tau * k^y
    absorb_nabla2 = absorb_eta .* k_grid.^(y - 1)    # dispersion: eta * k^(y-1)

    return AbsorptionParams(absorb_tau, absorb_eta, absorb_nabla1, absorb_nabla2, medium.alpha_mode)
end

"""
Compute k-space correction factor: kappa = sinc(c_ref * k * dt / 2) (unnormalized).
"""
function _compute_kappa(k_grid::AbstractArray, c_ref::Float64, dt::Float64)
    kappa = similar(k_grid)
    for i in eachindex(k_grid)
        arg = c_ref * k_grid[i] * dt / 2
        kappa[i] = arg ≈ 0 ? 1.0 : sin(arg) / arg
    end
    return kappa
end

# ============================================================================
# Solvers
# ============================================================================


"""
    kspace_first_order(kgrid, medium, source, sensor; kwargs...)

Run a k-space pseudospectral time-domain simulation.
Dispatches on grid dimensionality (1D, 2D, 3D).

# Arguments
- `kgrid`: `KWaveGrid1D`, `KWaveGrid2D`, or `KWaveGrid3D`
- `medium`: `KWaveMedium` with sound speed, density, and optional absorption
- `source`: `KWaveSource` with initial pressure and/or time-varying sources
- `sensor`: `KWaveSensor` defining recording locations and fields

# Keyword Arguments
- `pml_inside`: Place PML inside the grid (true) or expand grid (false). Default: `true`
- `pml_size`: PML thickness in grid points. Default: `20`
- `pml_alpha`: PML absorption coefficient. Default: `2.0`
- `smooth_p0`: Smooth initial pressure distribution. Default: `true`
- `smooth_c0`: Smooth sound speed field. Default: `false`
- `smooth_rho0`: Smooth density field. Default: `false`
- `data_cast`: Floating point type. Default: `Float64`
- `save_to_disk`: Save input to HDF5 file at this path instead of running. Default: `nothing`
- `plot_sim`: Display real-time simulation progress. Default: `false`
- `plot_layout`: Layout for simulation display. Default: `:default`
- `plot_scale`: Color scale for display. Default: `:auto`
- `record_movie`: Path to save movie file, or `nothing`. Default: `nothing`
- `progress_callback`: Optional callback `f(t_index, Nt, p_field)` called each step

# Returns
`SimulationOutput` containing recorded sensor data.

# Notes
- Call [`make_time!`](@ref) on `kgrid` before this function.
- Grid dimensionality is dispatched automatically from `kgrid` type.
- Use `data_cast=Float32` with GPU backends (CUDA / Metal / AMDGPU).
- Pass `save_to_disk="/path/input.h5"` to write an HDF5 input file for the k-Wave C++ binary
  instead of running the simulation.

# See Also
[`kspace_first_order_as`](@ref), [`acoustic_field_propagator`](@ref),
[`KWaveGrid`](@ref), [`make_time!`](@ref), [`KWaveMedium`](@ref),
[`KWaveSource`](@ref), [`KWaveSensor`](@ref), [`SimulationOutput`](@ref),
[`get_optimal_pml_size`](@ref)
"""
function kspace_first_order(
    kgrid::KWaveGrid1D,
    medium::KWaveMedium,
    source::KWaveSource,
    sensor::KWaveSensor;
    pml_inside::Bool=true,
    pml_size::Int=20,
    pml_alpha::Float64=2.0,
    smooth_p0::Bool=true,
    smooth_c0::Bool=false,
    smooth_rho0::Bool=false,
    data_cast::Type{T}=Float64,
    save_to_disk::Union{Nothing, String}=nothing,
    plot_sim::Bool=false,
    plot_layout::Symbol=:default,
    plot_scale::Union{Symbol, Tuple}=:auto,
    record_movie::Union{Nothing, String}=nothing,
    progress_callback::Union{Nothing, Function}=nothing,
    show_progress::Bool=true,
) where T <: AbstractFloat

    # ================================================================
    # 1. Input validation
    # ================================================================
    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running.")
    end
    validate_record_fields(sensor)
    _check_time_reversal = is_time_reversal(sensor)

    Nx = kgrid.Nx
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    # ================================================================
    # 2. Save-to-disk mode
    # ================================================================
    if save_to_disk !== nothing
        write_attributes(save_to_disk)
        write_grid(save_to_disk, kgrid; pml_size=(pml_size,), pml_alpha=(pml_alpha,))
        if medium.sound_speed isa AbstractArray
            write_matrix(save_to_disk, "c0", medium.sound_speed)
        else
            write_matrix(save_to_disk, "c_ref", medium.sound_speed)
        end
        if has_p0(source)
            write_matrix(save_to_disk, "p0", source.p0)
        end
        return SimulationOutput(Dict{Symbol, AbstractArray}())
    end

    # ================================================================
    # 3. Pre-compute PML
    # ================================================================
    pml_x = get_pml(Nx, kgrid.dx, pml_size, pml_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_size, pml_alpha, dt; staggered=true)

    # ================================================================
    # 4. Create FFT plans and pre-compute spectral operators
    # ================================================================
    plans = create_fft_plans(kgrid; data_cast=T)
    sops  = create_spectral_ops(kgrid; data_cast=T)
    n1    = rfft_first_dim(plans)   # = Nx÷2+1

    c_ref = medium.sound_speed isa Real ? medium.sound_speed : mean(medium.sound_speed)

    # rfft-shaped, T-typed kappa (spectral k-space correction factor)
    kappa_r = T.(_compute_kappa(kgrid.k, c_ref, dt))[1:n1]

    # rfft-sliced, T-typed absorption params (or nothing for lossless media)
    _absorb_f64 = _precompute_absorption(medium, kgrid.k, c_ref)
    absorb_r = if _absorb_f64 !== nothing
        AbsorptionParams(T(_absorb_f64.absorb_tau), T(_absorb_f64.absorb_eta),
                         T.(_absorb_f64.absorb_nabla1)[1:n1],
                         T.(_absorb_f64.absorb_nabla2)[1:n1],
                         _absorb_f64.mode)
    else
        nothing
    end

    # ================================================================
    # 5. Allocate field arrays (all in working precision T)
    # ================================================================
    p    = zeros(T, Nx)
    ux   = zeros(T, Nx)
    rhox = zeros(T, Nx)

    # rfft scratch: shape (Nx÷2+1,) — ~half the memory of a complex full-size array
    scratch1 = zeros(Complex{T}, n1)
    scratch2 = zeros(Complex{T}, n1)

    # Pre-allocated scratch buffers — eliminates per-step allocations in hot loop
    dpdx     = zeros(T, Nx)
    duxdx    = zeros(T, Nx)
    tmp_real = zeros(T, Nx)   # real staging buffer for _compute_pressure!

    # Pre-compute staggered density for heterogeneous media (time-invariant)
    rho0_sgx_1d = medium.density isa AbstractArray ?
        T.(_stagger_density_1d(medium.density)) : nothing

    # For dispersion: track previous total density
    rho_total_prev = (absorb_r !== nothing && absorb_r.mode != :no_dispersion) ?
        zeros(T, Nx) : nothing

    # ================================================================
    # 7. Sensor data
    # ================================================================
    sensor_data = _create_sensor_data(sensor, (Nx,), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{1}[]
    end

    # ================================================================
    # 8. Time-reversal setup
    # ================================================================
    tr_data = nothing
    tr_mask_indices = CartesianIndex{1}[]
    if _check_time_reversal
        tr_data = sensor.time_reversal_boundary_data
        Nt = size(tr_data, 2)
        kgrid.Nt[] = Nt
        # TR injection mask: sensor.mask defines where boundary data is injected
        tr_mask_indices = sensor.mask isa AbstractArray{Bool} ? findall(sensor.mask) : CartesianIndex{1}[]
    end

    # ================================================================
    # 9. Initialize from p0
    # ================================================================
    if !_check_time_reversal
        initialize_p0_1d!(p, ux, rhox, source, kgrid, medium, plans, scratch1, smooth_p0)
    end


    # ================================================================
    # 10. Time loop
    # ================================================================
    prog = Progress(Nt; desc="k-Wave 1D: ", enabled=show_progress)

    for t_index in 1:Nt
        # Time-reversal: inject boundary data at sensor locations (reversed in time)
        if _check_time_reversal
            tr_t = Nt - t_index + 1
            c0 = medium.sound_speed
            for (j, idx) in enumerate(tr_mask_indices)
                c_local = c0 isa Real ? c0 : c0[idx]
                p_val = tr_data[j, tr_t]
                rhox[idx] = p_val / c_local^2
            end
        end

        time_step_1d!(
            p, ux, rhox,
            scratch1, scratch2,
            dpdx, duxdx,
            kgrid, medium, source,
            pml_x, pml_x_sgx,
            kappa_r, plans, sops, t_index,
            absorb_r, tmp_real, rho_total_prev,
            rho0_sgx_1d,
        )

        if !isempty(mask_indices)
            velocities = (ux=ux,)
            record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices,
                               t_index, Nt, dt, medium.density)
        end

        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    # ================================================================
    # 11. Finalize
    # ================================================================
    finalize_sensor_data!(sensor_data, Nt)

    # In time-reversal mode, always return full-field p_final
    if _check_time_reversal
        sensor_data[:p_final] = copy(p)
    end
    apply_sensor_directivity!(sensor_data, sensor, kgrid)
    apply_frequency_response!(sensor_data, sensor, kgrid)

    return SimulationOutput(sensor_data)
end

# ============================================================================
# 2D Solver
# ============================================================================

function kspace_first_order(
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
    save_to_disk::Union{Nothing, String}=nothing,
    plot_sim::Bool=false,
    plot_layout::Symbol=:default,
    plot_scale::Union{Symbol, Tuple}=:auto,
    record_movie::Union{Nothing, String}=nothing,
    progress_callback::Union{Nothing, Function}=nothing,
    show_progress::Bool=true,
) where T <: AbstractFloat

    # ================================================================
    # 1. Input validation
    # ================================================================
    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running.")
    end
    validate_record_fields(sensor)
    _check_time_reversal = is_time_reversal(sensor)

    Nx, Ny = kgrid.Nx, kgrid.Ny
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    # Normalize PML parameters
    pml_x_size, pml_y_size = pml_size isa Tuple ? pml_size : (pml_size, pml_size)
    pml_x_alpha, pml_y_alpha = pml_alpha isa Tuple ? pml_alpha : (pml_alpha, pml_alpha)

    # ================================================================
    # 2. Save-to-disk mode
    # ================================================================
    if save_to_disk !== nothing
        write_attributes(save_to_disk)
        write_grid(save_to_disk, kgrid; pml_size=(pml_x_size, pml_y_size),
                   pml_alpha=(pml_x_alpha, pml_y_alpha))
        if medium.sound_speed isa AbstractArray
            write_matrix(save_to_disk, "c0", medium.sound_speed)
        else
            write_matrix(save_to_disk, "c_ref", medium.sound_speed)
        end
        if medium.density isa AbstractArray
            write_matrix(save_to_disk, "rho0", medium.density)
        end
        if has_p0(source)
            write_matrix(save_to_disk, "p0", source.p0)
        end
        return SimulationOutput(Dict{Symbol, AbstractArray}())
    end

    # ================================================================
    # 3. Pre-compute PML
    # ================================================================
    pml_x = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt)
    pml_y = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt; staggered=true)
    pml_y_sgy = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt; staggered=true)

    # ================================================================
    # 4. Create FFT plans and pre-compute spectral operators
    # ================================================================
    plans = create_fft_plans(kgrid; data_cast=T)
    sops  = create_spectral_ops(kgrid; data_cast=T)
    n1    = rfft_first_dim(plans)   # = Nx÷2+1

    c_ref = medium.sound_speed isa Real ? medium.sound_speed : mean(medium.sound_speed)

    # rfft-shaped, T-typed kappa
    kappa_r = T.(_compute_kappa(kgrid.k, c_ref, dt))[1:n1, :]

    # rfft-sliced, T-typed absorption params
    _absorb_f64 = _precompute_absorption(medium, kgrid.k, c_ref)
    absorb_r = if _absorb_f64 !== nothing
        AbsorptionParams(T(_absorb_f64.absorb_tau), T(_absorb_f64.absorb_eta),
                         T.(_absorb_f64.absorb_nabla1)[1:n1, :],
                         T.(_absorb_f64.absorb_nabla2)[1:n1, :],
                         _absorb_f64.mode)
    else
        nothing
    end

    # ================================================================
    # 5. Allocate field arrays (all in working precision T)
    # ================================================================
    p    = zeros(T, Nx, Ny)
    ux   = zeros(T, Nx, Ny)
    uy   = zeros(T, Nx, Ny)
    rhox = zeros(T, Nx, Ny)
    rhoy = zeros(T, Nx, Ny)

    # rfft scratch: shape (Nx÷2+1, Ny) — saves ~half the complex memory
    scratch1  = zeros(Complex{T}, n1, Ny)
    scratch2  = zeros(Complex{T}, n1, Ny)

    # Pre-allocated scratch buffers — eliminates per-step allocations in hot loop
    dpdx      = zeros(T, Nx, Ny)
    dpdy      = zeros(T, Nx, Ny)
    duxdx     = zeros(T, Nx, Ny)
    duydy     = zeros(T, Nx, Ny)
    rho_total = zeros(T, Nx, Ny)
    tmp_real  = zeros(T, Nx, Ny)   # real staging buffer for _compute_pressure!

    # Pre-reshape PML vectors — eliminates per-step wrapper object allocations
    pml_x_col     = reshape(pml_x,     :, 1)
    pml_y_row     = reshape(pml_y,     1, :)
    pml_x_sgx_col = reshape(pml_x_sgx, :, 1)
    pml_y_sgy_row = reshape(pml_y_sgy, 1, :)

    # Pre-compute staggered density for heterogeneous media (time-invariant)
    rho0_sgx_2d = medium.density isa AbstractArray ?
        T.(_stagger_density_2d(medium.density, 1)) : nothing
    rho0_sgy_2d = medium.density isa AbstractArray ?
        T.(_stagger_density_2d(medium.density, 2)) : nothing

    rho_total_prev = (absorb_r !== nothing && absorb_r.mode != :no_dispersion) ?
        zeros(T, Nx, Ny) : nothing

    # ================================================================
    # 7. Sensor data
    # ================================================================
    sensor_data = _create_sensor_data(sensor, (Nx, Ny), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{2}[]
    end

    # ================================================================
    # 8. Time-reversal setup
    # ================================================================
    tr_data = nothing
    tr_mask_indices = CartesianIndex{2}[]
    if _check_time_reversal
        tr_data = sensor.time_reversal_boundary_data
        Nt = size(tr_data, 2)
        kgrid.Nt[] = Nt
        tr_mask_indices = sensor.mask isa AbstractArray{Bool} ? findall(sensor.mask) : CartesianIndex{2}[]
    end

    # ================================================================
    # 9. Initialize from p0
    # ================================================================
    if !_check_time_reversal
        initialize_p0_2d!(p, ux, uy, rhox, rhoy, source, kgrid, medium, plans, scratch1, smooth_p0)
    end


    # ================================================================
    # 10. Time loop
    # ================================================================
    prog = Progress(Nt; desc="k-Wave 2D: ", enabled=show_progress)

    for t_index in 1:Nt
        # Time-reversal: inject boundary data at sensor locations (reversed in time)
        if _check_time_reversal
            tr_t = Nt - t_index + 1
            c0 = medium.sound_speed
            for (j, idx) in enumerate(tr_mask_indices)
                c_local = c0 isa Real ? c0 : c0[idx]
                p_val = tr_data[j, tr_t]
                rho_val = p_val / (2 * c_local^2)
                rhox[idx] = rho_val
                rhoy[idx] = rho_val
            end
        end

        time_step_2d!(
            p, ux, uy, rhox, rhoy,
            scratch1, scratch2,
            dpdx, dpdy, duxdx, duydy, rho_total,
            kgrid, medium, source,
            pml_x_col, pml_y_row, pml_x_sgx_col, pml_y_sgy_row,
            kappa_r, plans, sops, t_index,
            absorb_r, tmp_real, rho_total_prev,
            rho0_sgx_2d, rho0_sgy_2d,
        )

        if !isempty(mask_indices)
            velocities = (ux=ux, uy=uy)
            record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices,
                               t_index, Nt, dt, medium.density)
        end

        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    # ================================================================
    # 11. Finalize
    # ================================================================
    finalize_sensor_data!(sensor_data, Nt)

    # In time-reversal mode, always return full-field p_final
    if _check_time_reversal
        sensor_data[:p_final] = reshape(copy(p), :)
    end
    apply_sensor_directivity!(sensor_data, sensor, kgrid)
    apply_frequency_response!(sensor_data, sensor, kgrid)

    return SimulationOutput(sensor_data)
end

# ============================================================================
# 3D Solver
# ============================================================================

function kspace_first_order(
    kgrid::KWaveGrid3D,
    medium::KWaveMedium,
    source::KWaveSource,
    sensor::KWaveSensor;
    pml_inside::Bool=true,
    pml_size::Union{Int, NTuple{3, Int}}=20,
    pml_alpha::Union{Float64, NTuple{3, Float64}}=2.0,
    smooth_p0::Bool=true,
    smooth_c0::Bool=false,
    smooth_rho0::Bool=false,
    data_cast::Type{T}=Float64,
    save_to_disk::Union{Nothing, String}=nothing,
    plot_sim::Bool=false,
    plot_layout::Symbol=:default,
    plot_scale::Union{Symbol, Tuple}=:auto,
    record_movie::Union{Nothing, String}=nothing,
    progress_callback::Union{Nothing, Function}=nothing,
    show_progress::Bool=true,
) where T <: AbstractFloat

    # ================================================================
    # 1. Input validation
    # ================================================================
    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running.")
    end
    validate_record_fields(sensor)
    _check_time_reversal = is_time_reversal(sensor)

    Nx, Ny, Nz = kgrid.Nx, kgrid.Ny, kgrid.Nz
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    # Normalize PML parameters
    pml_x_size, pml_y_size, pml_z_size = pml_size isa Tuple ? pml_size : (pml_size, pml_size, pml_size)
    pml_x_alpha, pml_y_alpha, pml_z_alpha = pml_alpha isa Tuple ? pml_alpha : (pml_alpha, pml_alpha, pml_alpha)

    # ================================================================
    # 2. Save-to-disk mode
    # ================================================================
    if save_to_disk !== nothing
        write_attributes(save_to_disk)
        write_grid(save_to_disk, kgrid; pml_size=(pml_x_size, pml_y_size, pml_z_size),
                   pml_alpha=(pml_x_alpha, pml_y_alpha, pml_z_alpha))
        if medium.sound_speed isa AbstractArray
            write_matrix(save_to_disk, "c0", medium.sound_speed)
        else
            write_matrix(save_to_disk, "c_ref", medium.sound_speed)
        end
        if medium.density isa AbstractArray
            write_matrix(save_to_disk, "rho0", medium.density)
        end
        if has_p0(source)
            write_matrix(save_to_disk, "p0", source.p0)
        end
        return SimulationOutput(Dict{Symbol, AbstractArray}())
    end

    # ================================================================
    # 3. Pre-compute PML
    # ================================================================
    pml_x = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt)
    pml_y = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt)
    pml_z = get_pml(Nz, kgrid.dz, pml_z_size, pml_z_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt; staggered=true)
    pml_y_sgy = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt; staggered=true)
    pml_z_sgz = get_pml(Nz, kgrid.dz, pml_z_size, pml_z_alpha, dt; staggered=true)

    # ================================================================
    # 4. Create FFT plans and pre-compute spectral operators
    # ================================================================
    plans = create_fft_plans(kgrid; data_cast=T)
    sops  = create_spectral_ops(kgrid; data_cast=T)
    n1    = rfft_first_dim(plans)   # = Nx÷2+1

    c_ref = medium.sound_speed isa Real ? medium.sound_speed : mean(medium.sound_speed)

    # rfft-shaped, T-typed kappa
    kappa_r = T.(_compute_kappa(kgrid.k, c_ref, dt))[1:n1, :, :]

    # rfft-sliced, T-typed absorption params
    _absorb_f64 = _precompute_absorption(medium, kgrid.k, c_ref)
    absorb_r = if _absorb_f64 !== nothing
        AbsorptionParams(T(_absorb_f64.absorb_tau), T(_absorb_f64.absorb_eta),
                         T.(_absorb_f64.absorb_nabla1)[1:n1, :, :],
                         T.(_absorb_f64.absorb_nabla2)[1:n1, :, :],
                         _absorb_f64.mode)
    else
        nothing
    end

    # ================================================================
    # 5. Allocate field arrays (all in working precision T)
    # ================================================================
    p    = zeros(T, Nx, Ny, Nz)
    ux   = zeros(T, Nx, Ny, Nz)
    uy   = zeros(T, Nx, Ny, Nz)
    uz   = zeros(T, Nx, Ny, Nz)
    rhox = zeros(T, Nx, Ny, Nz)
    rhoy = zeros(T, Nx, Ny, Nz)
    rhoz = zeros(T, Nx, Ny, Nz)

    # rfft scratch: shape (Nx÷2+1, Ny, Nz) — saves ~half the complex memory
    scratch1  = zeros(Complex{T}, n1, Ny, Nz)
    scratch2  = zeros(Complex{T}, n1, Ny, Nz)

    # Pre-allocated scratch buffers — eliminates per-step allocations in hot loop
    dpdx      = zeros(T, Nx, Ny, Nz)
    dpdy      = zeros(T, Nx, Ny, Nz)
    dpdz      = zeros(T, Nx, Ny, Nz)
    duxdx     = zeros(T, Nx, Ny, Nz)
    duydy     = zeros(T, Nx, Ny, Nz)
    duzdz     = zeros(T, Nx, Ny, Nz)
    rho_total = zeros(T, Nx, Ny, Nz)
    tmp_real  = zeros(T, Nx, Ny, Nz)   # real staging buffer for _compute_pressure!

    # Pre-reshape PML vectors — eliminates per-step wrapper object allocations
    pml_x_r     = reshape(pml_x,     :, 1, 1)
    pml_y_r     = reshape(pml_y,     1, :, 1)
    pml_z_r     = reshape(pml_z,     1, 1, :)
    pml_x_sgx_r = reshape(pml_x_sgx, :, 1, 1)
    pml_y_sgy_r = reshape(pml_y_sgy, 1, :, 1)
    pml_z_sgz_r = reshape(pml_z_sgz, 1, 1, :)

    # Pre-compute staggered density for heterogeneous media (time-invariant)
    rho0_sgx_3d = medium.density isa AbstractArray ?
        T.(_stagger_density_3d(medium.density, 1)) : nothing
    rho0_sgy_3d = medium.density isa AbstractArray ?
        T.(_stagger_density_3d(medium.density, 2)) : nothing
    rho0_sgz_3d = medium.density isa AbstractArray ?
        T.(_stagger_density_3d(medium.density, 3)) : nothing

    rho_total_prev = (absorb_r !== nothing && absorb_r.mode != :no_dispersion) ?
        zeros(T, Nx, Ny, Nz) : nothing

    # ================================================================
    # 7. Sensor data
    # ================================================================
    sensor_data = _create_sensor_data(sensor, (Nx, Ny, Nz), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{3}[]
    end

    # ================================================================
    # 8. Time-reversal setup
    # ================================================================
    tr_data = nothing
    tr_mask_indices = CartesianIndex{3}[]
    if _check_time_reversal
        tr_data = sensor.time_reversal_boundary_data
        Nt = size(tr_data, 2)
        kgrid.Nt[] = Nt
        tr_mask_indices = sensor.mask isa AbstractArray{Bool} ? findall(sensor.mask) : CartesianIndex{3}[]
    end

    # ================================================================
    # 9. Initialize from p0
    # ================================================================
    if !_check_time_reversal
        initialize_p0_3d!(p, ux, uy, uz, rhox, rhoy, rhoz,
                          source, kgrid, medium, plans, scratch1, smooth_p0)
    end


    # ================================================================
    # 10. Time loop
    # ================================================================
    prog = Progress(Nt; desc="k-Wave 3D: ", enabled=show_progress)

    for t_index in 1:Nt
        if _check_time_reversal
            tr_t = Nt - t_index + 1
            c0 = medium.sound_speed
            for (j, idx) in enumerate(tr_mask_indices)
                c_local = c0 isa Real ? c0 : c0[idx]
                p_val = tr_data[j, tr_t]
                rho_val = p_val / (3 * c_local^2)
                rhox[idx] = rho_val
                rhoy[idx] = rho_val
                rhoz[idx] = rho_val
            end
        end

        time_step_3d!(
            p, ux, uy, uz, rhox, rhoy, rhoz,
            scratch1, scratch2,
            dpdx, dpdy, dpdz, duxdx, duydy, duzdz, rho_total,
            kgrid, medium, source,
            pml_x_r, pml_y_r, pml_z_r, pml_x_sgx_r, pml_y_sgy_r, pml_z_sgz_r,
            kappa_r, plans, sops, t_index,
            absorb_r, tmp_real, rho_total_prev,
            rho0_sgx_3d, rho0_sgy_3d, rho0_sgz_3d,
        )

        if !isempty(mask_indices)
            velocities = (ux=ux, uy=uy, uz=uz)
            record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices,
                               t_index, Nt, dt, medium.density)
        end

        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    # ================================================================
    # 11. Finalize
    # ================================================================
    finalize_sensor_data!(sensor_data, Nt)

    # In time-reversal mode, always return full-field p_final
    if _check_time_reversal
        sensor_data[:p_final] = reshape(copy(p), :)
    end
    apply_sensor_directivity!(sensor_data, sensor, kgrid)
    apply_frequency_response!(sensor_data, sensor, kgrid)

    return SimulationOutput(sensor_data)
end
