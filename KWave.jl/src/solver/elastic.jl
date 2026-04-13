# ============================================================================
# KWave.jl — Elastic wave solvers (pstd_elastic_2d, pstd_elastic_3d)
# ============================================================================
# Pseudospectral time-domain elastic wave propagation solver.
# Supports compressional + shear waves, viscoelastic absorption,
# and stress tensor sources.
# ============================================================================

"""
    ElasticMedium{T}

Elastic medium properties for PSTD elastic simulations.

# Fields
- `sound_speed_compression`: Compressional wave speed [m/s]
- `sound_speed_shear`: Shear wave speed [m/s]
- `density`: Mass density [kg/m³]
- `alpha_coeff_compression`: Absorption for compressional waves [dB/(MHz^y cm)]
- `alpha_coeff_shear`: Absorption for shear waves [dB/(MHz^y cm)]
- `alpha_power`: Absorption power law exponent
"""
struct ElasticMedium{T<:AbstractFloat}
    sound_speed_compression::Union{T, AbstractArray{T}}
    sound_speed_shear::Union{T, AbstractArray{T}}
    density::Union{T, AbstractArray{T}}
    alpha_coeff_compression::Union{Nothing, T, AbstractArray{T}}
    alpha_coeff_shear::Union{Nothing, T, AbstractArray{T}}
    alpha_power::Union{Nothing, T}
end

function ElasticMedium(;
    sound_speed_compression::Union{Real, AbstractArray{<:Real}},
    sound_speed_shear::Union{Real, AbstractArray{<:Real}},
    density::Union{Real, AbstractArray{<:Real}}=1000.0,
    alpha_coeff_compression::Union{Nothing, Real, AbstractArray{<:Real}}=nothing,
    alpha_coeff_shear::Union{Nothing, Real, AbstractArray{<:Real}}=nothing,
    alpha_power::Union{Nothing, Real}=nothing)

    _f64 = x -> x === nothing ? nothing : (x isa AbstractArray ? Float64.(x) : Float64(x))
    return ElasticMedium{Float64}(
        _f64(sound_speed_compression),
        _f64(sound_speed_shear),
        _f64(density),
        _f64(alpha_coeff_compression),
        _f64(alpha_coeff_shear),
        alpha_power === nothing ? nothing : Float64(alpha_power),
    )
end

"""
    ElasticSource{T}

Source definitions for elastic simulations.

Supports stress tensor sources (`sxx`, `syy`, `sxy`, etc.) and velocity sources.
"""
struct ElasticSource{T<:AbstractFloat}
    # Stress sources
    s_mask::Union{Nothing, AbstractArray{Bool}}
    sxx::Union{Nothing, AbstractArray{T}}
    syy::Union{Nothing, AbstractArray{T}}
    szz::Union{Nothing, AbstractArray{T}}
    sxy::Union{Nothing, AbstractArray{T}}
    sxz::Union{Nothing, AbstractArray{T}}
    syz::Union{Nothing, AbstractArray{T}}
    s_mode::SourceMode
    # Velocity sources
    u_mask::Union{Nothing, AbstractArray{Bool}}
    ux::Union{Nothing, AbstractArray{T}}
    uy::Union{Nothing, AbstractArray{T}}
    uz::Union{Nothing, AbstractArray{T}}
    u_mode::SourceMode
    # Initial pressure
    p0::Union{Nothing, AbstractArray{T}}
end

function ElasticSource(;
    s_mask::Union{Nothing, AbstractArray{Bool}}=nothing,
    sxx::Union{Nothing, AbstractArray{<:Real}}=nothing,
    syy::Union{Nothing, AbstractArray{<:Real}}=nothing,
    szz::Union{Nothing, AbstractArray{<:Real}}=nothing,
    sxy::Union{Nothing, AbstractArray{<:Real}}=nothing,
    sxz::Union{Nothing, AbstractArray{<:Real}}=nothing,
    syz::Union{Nothing, AbstractArray{<:Real}}=nothing,
    s_mode::SourceMode=Additive,
    u_mask::Union{Nothing, AbstractArray{Bool}}=nothing,
    ux::Union{Nothing, AbstractArray{<:Real}}=nothing,
    uy::Union{Nothing, AbstractArray{<:Real}}=nothing,
    uz::Union{Nothing, AbstractArray{<:Real}}=nothing,
    u_mode::SourceMode=Additive,
    p0::Union{Nothing, AbstractArray{<:Real}}=nothing)

    _f64 = x -> x === nothing ? nothing : Float64.(x)
    return ElasticSource{Float64}(
        s_mask, _f64(sxx), _f64(syy), _f64(szz),
        _f64(sxy), _f64(sxz), _f64(syz), s_mode,
        u_mask, _f64(ux), _f64(uy), _f64(uz), u_mode,
        _f64(p0),
    )
end

has_stress_source(s::ElasticSource) = s.s_mask !== nothing
has_velocity_source(s::ElasticSource) = s.u_mask !== nothing
has_p0(s::ElasticSource) = s.p0 !== nothing

# ============================================================================
# 2D Elastic Solver
# ============================================================================

"""
    pstd_elastic_2d(kgrid, medium, source, sensor; kwargs...)

Simulate 2D elastic wave propagation using the pseudospectral time-domain method.

Solves the coupled compressional and shear wave equations on a staggered grid.
The stress-velocity formulation uses the full elastic wave equation with
optional viscoelastic absorption (Kelvin-Voigt model).

# Arguments
- `kgrid::KWaveGrid2D`: 2D computational grid
- `medium::ElasticMedium`: Elastic medium properties
- `source::ElasticSource`: Source definitions (stress tensor and/or velocity)
- `sensor::KWaveSensor`: Sensor/detector definitions

# Keyword Arguments
- `pml_size`: PML thickness in grid points (default: 20)
- `pml_alpha`: PML absorption coefficient (default: 2.0)
- `data_cast`: Float type for computation (default: Float64)

# Returns
`SimulationOutput` containing recorded sensor data.
"""
function pstd_elastic_2d(
    kgrid::KWaveGrid2D,
    medium::ElasticMedium,
    source::ElasticSource,
    sensor::KWaveSensor;
    pml_size::Union{Int, Tuple{Int, Int}}=20,
    pml_alpha::Union{Float64, Tuple{Float64, Float64}}=2.0,
    data_cast::Type{T}=Float64,
    plot_sim::Bool=false,
    progress_callback::Union{Nothing, Function}=nothing,
) where T <: AbstractFloat

    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running.")
    end
    validate_record_fields(sensor)

    Nx, Ny = kgrid.Nx, kgrid.Ny
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    pml_x_size, pml_y_size = pml_size isa Tuple ? pml_size : (pml_size, pml_size)
    pml_x_alpha, pml_y_alpha = pml_alpha isa Tuple ? pml_alpha : (pml_alpha, pml_alpha)

    # PML profiles
    pml_x = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt)
    pml_y = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt; staggered=true)
    pml_y_sgy = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt; staggered=true)

    pml_x_col = reshape(pml_x, :, 1)
    pml_y_row = reshape(pml_y, 1, :)
    pml_x_sgx_col = reshape(pml_x_sgx, :, 1)
    pml_y_sgy_row = reshape(pml_y_sgy, 1, :)

    # Material parameters
    cp = medium.sound_speed_compression
    cs = medium.sound_speed_shear
    rho0 = medium.density

    # Lamé parameters: λ = ρ(cp² - 2cs²), μ = ρcs²
    lambda = rho0 isa Real ?
        rho0 * (cp^2 - 2 * cs^2) :
        rho0 .* (cp.^2 .- 2 .* cs.^2)
    mu = rho0 isa Real ? rho0 * cs^2 : rho0 .* cs.^2

    # Viscoelastic absorption (Kelvin-Voigt)
    if medium.alpha_coeff_compression !== nothing && medium.alpha_power !== nothing
        c_ref = cp isa Real ? cp : mean(cp)
        y = medium.alpha_power
        alpha_p_nepers = db2neper(medium.alpha_coeff_compression, y)
        chi_p = 2 * alpha_p_nepers * c_ref^(y - 1)
    else
        chi_p = 0.0
    end

    if medium.alpha_coeff_shear !== nothing && medium.alpha_power !== nothing
        cs_ref = cs isa Real ? cs : mean(cs)
        y = medium.alpha_power
        alpha_s_nepers = db2neper(medium.alpha_coeff_shear, y)
        chi_s = 2 * alpha_s_nepers * cs_ref^(y - 1)
    else
        chi_s = 0.0
    end

    # FFT plans
    plans = create_fft_plans(kgrid; data_cast=T)
    scratch = zeros(ComplexF64, Nx, Ny)

    # Field arrays — stress-velocity formulation
    vx = zeros(Float64, Nx, Ny)      # x-velocity
    vy = zeros(Float64, Nx, Ny)      # y-velocity
    sxx = zeros(Float64, Nx, Ny)     # normal stress xx
    syy = zeros(Float64, Nx, Ny)     # normal stress yy
    sxy = zeros(Float64, Nx, Ny)     # shear stress xy

    # Sensor data
    sensor_data = _create_sensor_data(sensor, (Nx, Ny), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{2}[]
    end

    # Initialize from p0 (isotropic initial pressure → sxx = syy = -p0)
    if has_p0(source)
        p0 = source.p0
        @. sxx = -p0
        @. syy = -p0
    end

    prog = Progress(Nt; desc="PSTD Elastic 2D: ", enabled=true)

    for t_index in 1:Nt
        # === STEP 1: Compute stress gradients ===
        dsxx_dx = similar(sxx)
        dsxy_dy = similar(sxy)
        dsxy_dx = similar(sxy)
        dsyy_dy = similar(syy)

        spectral_gradient!(dsxx_dx, sxx, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 2)
        spectral_gradient!(dsxy_dy, sxy, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 2)
        spectral_gradient!(dsxy_dx, sxy, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 2)
        spectral_gradient!(dsyy_dy, syy, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 2)

        # === STEP 2: Velocity update (momentum equation) ===
        if rho0 isa Real
            @. vx = pml_x_sgx_col * (pml_x_sgx_col * vx + dt / rho0 * (dsxx_dx + dsxy_dy))
            @. vy = pml_y_sgy_row * (pml_y_sgy_row * vy + dt / rho0 * (dsxy_dx + dsyy_dy))
        else
            rho0_sgx = _stagger_density_2d(rho0, 1)
            rho0_sgy = _stagger_density_2d(rho0, 2)
            @. vx = pml_x_sgx_col * (pml_x_sgx_col * vx + dt / rho0_sgx * (dsxx_dx + dsxy_dy))
            @. vy = pml_y_sgy_row * (pml_y_sgy_row * vy + dt / rho0_sgy * (dsxy_dx + dsyy_dy))
        end

        # === STEP 3: Velocity source injection ===
        if has_velocity_source(source)
            u_indices = findall(source.u_mask)
            if source.ux !== nothing
                sig_x = size(source.ux, 2) == 1 ? source.ux[:, 1] : source.ux[:, min(t_index, size(source.ux, 2))]
                for (j, idx) in enumerate(u_indices)
                    if source.u_mode == Dirichlet
                        vx[idx] = sig_x[j]
                    else
                        vx[idx] += sig_x[j]
                    end
                end
            end
            if source.uy !== nothing
                sig_y = size(source.uy, 2) == 1 ? source.uy[:, 1] : source.uy[:, min(t_index, size(source.uy, 2))]
                for (j, idx) in enumerate(u_indices)
                    if source.u_mode == Dirichlet
                        vy[idx] = sig_y[j]
                    else
                        vy[idx] += sig_y[j]
                    end
                end
            end
        end

        # === STEP 4: Velocity gradients ===
        dvx_dx = similar(vx)
        dvy_dy = similar(vy)
        dvx_dy = similar(vx)
        dvy_dx = similar(vy)

        spectral_gradient!(dvx_dx, vx, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch, plans, 1, 2)
        spectral_gradient!(dvy_dy, vy, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch, plans, 2, 2)
        spectral_gradient!(dvx_dy, vx, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch, plans, 2, 2)
        spectral_gradient!(dvy_dx, vy, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch, plans, 1, 2)

        # === STEP 5: Stress update (constitutive relation) ===
        # σxx = (λ + 2μ) dvx/dx + λ dvy/dy  (+ viscous terms)
        # σyy = λ dvx/dx + (λ + 2μ) dvy/dy
        # σxy = μ (dvx/dy + dvy/dx)
        if lambda isa Real
            @. sxx = pml_x_col * (pml_x_col * sxx + dt * (
                (lambda + 2*mu) * dvx_dx + lambda * dvy_dy +
                chi_p * ((lambda + 2*mu) * dvx_dx + lambda * dvy_dy)
            ))
            @. syy = pml_y_row * (pml_y_row * syy + dt * (
                lambda * dvx_dx + (lambda + 2*mu) * dvy_dy +
                chi_p * (lambda * dvx_dx + (lambda + 2*mu) * dvy_dy)
            ))
            @. sxy = pml_x_col * pml_y_row * (sxy + dt * (
                mu * (dvx_dy + dvy_dx) + chi_s * mu * (dvx_dy + dvy_dx)
            ))
        else
            @. sxx = pml_x_col * (pml_x_col * sxx + dt * (
                (lambda + 2*mu) * dvx_dx + lambda * dvy_dy
            ))
            @. syy = pml_y_row * (pml_y_row * syy + dt * (
                lambda * dvx_dx + (lambda + 2*mu) * dvy_dy
            ))
            @. sxy = pml_x_col * pml_y_row * (sxy + dt * (
                mu .* (dvx_dy + dvy_dx)
            ))
        end

        # === STEP 6: Stress source injection ===
        if has_stress_source(source)
            s_indices = findall(source.s_mask)
            if source.sxx !== nothing
                sig = size(source.sxx, 2) == 1 ? source.sxx[:, 1] : source.sxx[:, min(t_index, size(source.sxx, 2))]
                for (j, idx) in enumerate(s_indices)
                    if source.s_mode == Dirichlet
                        sxx[idx] = sig[j]
                    else
                        sxx[idx] += sig[j]
                    end
                end
            end
            if source.syy !== nothing
                sig = size(source.syy, 2) == 1 ? source.syy[:, 1] : source.syy[:, min(t_index, size(source.syy, 2))]
                for (j, idx) in enumerate(s_indices)
                    if source.s_mode == Dirichlet
                        syy[idx] = sig[j]
                    else
                        syy[idx] += sig[j]
                    end
                end
            end
            if source.sxy !== nothing
                sig = size(source.sxy, 2) == 1 ? source.sxy[:, 1] : source.sxy[:, min(t_index, size(source.sxy, 2))]
                for (j, idx) in enumerate(s_indices)
                    if source.s_mode == Dirichlet
                        sxy[idx] = sig[j]
                    else
                        sxy[idx] += sig[j]
                    end
                end
            end
        end

        # === Record sensor data ===
        # Pressure = -(sxx + syy) / 2
        p = @. -(sxx + syy) / 2
        if !isempty(mask_indices)
            velocities = (ux=vx, uy=vy)
            record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices,
                               t_index, Nt, dt, medium.density)
        end

        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    finalize_sensor_data!(sensor_data, Nt)
    return SimulationOutput(sensor_data)
end

# ============================================================================
# 3D Elastic Solver
# ============================================================================

"""
    pstd_elastic_3d(kgrid, medium, source, sensor; kwargs...)

Simulate 3D elastic wave propagation using the pseudospectral time-domain method.

# Arguments
- `kgrid::KWaveGrid3D`: 3D computational grid
- `medium::ElasticMedium`: Elastic medium properties
- `source::ElasticSource`: Source definitions
- `sensor::KWaveSensor`: Sensor/detector definitions

# Returns
`SimulationOutput` containing recorded sensor data.
"""
function pstd_elastic_3d(
    kgrid::KWaveGrid3D,
    medium::ElasticMedium,
    source::ElasticSource,
    sensor::KWaveSensor;
    pml_size::Union{Int, NTuple{3, Int}}=10,
    pml_alpha::Union{Float64, NTuple{3, Float64}}=2.0,
    data_cast::Type{T}=Float64,
    plot_sim::Bool=false,
    progress_callback::Union{Nothing, Function}=nothing,
) where T <: AbstractFloat

    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running.")
    end
    validate_record_fields(sensor)

    Nx, Ny, Nz = kgrid.Nx, kgrid.Ny, kgrid.Nz
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    pml_x_size, pml_y_size, pml_z_size = pml_size isa Tuple ? pml_size : (pml_size, pml_size, pml_size)
    pml_x_alpha, pml_y_alpha, pml_z_alpha = pml_alpha isa Tuple ? pml_alpha : (pml_alpha, pml_alpha, pml_alpha)

    # PML
    pml_x = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt)
    pml_y = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt)
    pml_z = get_pml(Nz, kgrid.dz, pml_z_size, pml_z_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt; staggered=true)
    pml_y_sgy = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt; staggered=true)
    pml_z_sgz = get_pml(Nz, kgrid.dz, pml_z_size, pml_z_alpha, dt; staggered=true)

    pml_x_3d = reshape(pml_x, :, 1, 1)
    pml_y_3d = reshape(pml_y, 1, :, 1)
    pml_z_3d = reshape(pml_z, 1, 1, :)
    pml_x_sgx_3d = reshape(pml_x_sgx, :, 1, 1)
    pml_y_sgy_3d = reshape(pml_y_sgy, 1, :, 1)
    pml_z_sgz_3d = reshape(pml_z_sgz, 1, 1, :)

    # Material
    cp = medium.sound_speed_compression
    cs = medium.sound_speed_shear
    rho0 = medium.density

    lambda = rho0 isa Real ?
        rho0 * (cp^2 - 2 * cs^2) :
        rho0 .* (cp.^2 .- 2 .* cs.^2)
    mu = rho0 isa Real ? rho0 * cs^2 : rho0 .* cs.^2

    # FFT plans
    plans = create_fft_plans(kgrid; data_cast=T)
    scratch = zeros(ComplexF64, Nx, Ny, Nz)

    # Field arrays
    vx = zeros(Float64, Nx, Ny, Nz)
    vy = zeros(Float64, Nx, Ny, Nz)
    vz = zeros(Float64, Nx, Ny, Nz)
    sxx_f = zeros(Float64, Nx, Ny, Nz)
    syy_f = zeros(Float64, Nx, Ny, Nz)
    szz_f = zeros(Float64, Nx, Ny, Nz)
    sxy_f = zeros(Float64, Nx, Ny, Nz)
    sxz_f = zeros(Float64, Nx, Ny, Nz)
    syz_f = zeros(Float64, Nx, Ny, Nz)

    # Sensor data
    sensor_data = _create_sensor_data(sensor, (Nx, Ny, Nz), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{3}[]
    end

    # Initialize from p0
    if has_p0(source)
        p0 = source.p0
        @. sxx_f = -p0
        @. syy_f = -p0
        @. szz_f = -p0
    end

    prog = Progress(Nt; desc="PSTD Elastic 3D: ", enabled=true)

    for t_index in 1:Nt
        # === Stress gradients ===
        dsxx_dx = similar(sxx_f); dsxy_dy = similar(sxy_f); dsxz_dz = similar(sxz_f)
        dsxy_dx = similar(sxy_f); dsyy_dy = similar(syy_f); dsyz_dz = similar(syz_f)
        dsxz_dx = similar(sxz_f); dsyz_dy = similar(syz_f); dszz_dz = similar(szz_f)

        spectral_gradient!(dsxx_dx, sxx_f, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 3)
        spectral_gradient!(dsxy_dy, sxy_f, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 3)
        spectral_gradient!(dsxz_dz, sxz_f, kgrid.kz_vec, kgrid.ddz_k_shift_pos, scratch, plans, 3, 3)

        spectral_gradient!(dsxy_dx, sxy_f, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 3)
        spectral_gradient!(dsyy_dy, syy_f, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 3)
        spectral_gradient!(dsyz_dz, syz_f, kgrid.kz_vec, kgrid.ddz_k_shift_pos, scratch, plans, 3, 3)

        spectral_gradient!(dsxz_dx, sxz_f, kgrid.kx_vec, kgrid.ddx_k_shift_pos, scratch, plans, 1, 3)
        spectral_gradient!(dsyz_dy, syz_f, kgrid.ky_vec, kgrid.ddy_k_shift_pos, scratch, plans, 2, 3)
        spectral_gradient!(dszz_dz, szz_f, kgrid.kz_vec, kgrid.ddz_k_shift_pos, scratch, plans, 3, 3)

        # === Velocity update ===
        if rho0 isa Real
            @. vx = pml_x_sgx_3d * (pml_x_sgx_3d * vx + dt / rho0 * (dsxx_dx + dsxy_dy + dsxz_dz))
            @. vy = pml_y_sgy_3d * (pml_y_sgy_3d * vy + dt / rho0 * (dsxy_dx + dsyy_dy + dsyz_dz))
            @. vz = pml_z_sgz_3d * (pml_z_sgz_3d * vz + dt / rho0 * (dsxz_dx + dsyz_dy + dszz_dz))
        else
            rho0_sgx = _stagger_density_3d(rho0, 1)
            rho0_sgy = _stagger_density_3d(rho0, 2)
            rho0_sgz = _stagger_density_3d(rho0, 3)
            @. vx = pml_x_sgx_3d * (pml_x_sgx_3d * vx + dt / rho0_sgx * (dsxx_dx + dsxy_dy + dsxz_dz))
            @. vy = pml_y_sgy_3d * (pml_y_sgy_3d * vy + dt / rho0_sgy * (dsxy_dx + dsyy_dy + dsyz_dz))
            @. vz = pml_z_sgz_3d * (pml_z_sgz_3d * vz + dt / rho0_sgz * (dsxz_dx + dsyz_dy + dszz_dz))
        end

        # === Velocity source injection ===
        if has_velocity_source(source)
            u_indices = findall(source.u_mask)
            for (field_arr, src_data) in [(vx, source.ux), (vy, source.uy), (vz, source.uz)]
                if src_data !== nothing
                    sig = size(src_data, 2) == 1 ? src_data[:, 1] : src_data[:, min(t_index, size(src_data, 2))]
                    for (j, idx) in enumerate(u_indices)
                        if source.u_mode == Dirichlet
                            field_arr[idx] = sig[j]
                        else
                            field_arr[idx] += sig[j]
                        end
                    end
                end
            end
        end

        # === Velocity gradients ===
        dvx_dx = similar(vx); dvx_dy = similar(vx); dvx_dz = similar(vx)
        dvy_dx = similar(vy); dvy_dy = similar(vy); dvy_dz = similar(vy)
        dvz_dx = similar(vz); dvz_dy = similar(vz); dvz_dz = similar(vz)

        spectral_gradient!(dvx_dx, vx, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch, plans, 1, 3)
        spectral_gradient!(dvx_dy, vx, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch, plans, 2, 3)
        spectral_gradient!(dvx_dz, vx, kgrid.kz_vec, kgrid.ddz_k_shift_neg, scratch, plans, 3, 3)

        spectral_gradient!(dvy_dx, vy, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch, plans, 1, 3)
        spectral_gradient!(dvy_dy, vy, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch, plans, 2, 3)
        spectral_gradient!(dvy_dz, vy, kgrid.kz_vec, kgrid.ddz_k_shift_neg, scratch, plans, 3, 3)

        spectral_gradient!(dvz_dx, vz, kgrid.kx_vec, kgrid.ddx_k_shift_neg, scratch, plans, 1, 3)
        spectral_gradient!(dvz_dy, vz, kgrid.ky_vec, kgrid.ddy_k_shift_neg, scratch, plans, 2, 3)
        spectral_gradient!(dvz_dz, vz, kgrid.kz_vec, kgrid.ddz_k_shift_neg, scratch, plans, 3, 3)

        # === Stress update ===
        if lambda isa Real
            @. sxx_f = pml_x_3d * (pml_x_3d * sxx_f + dt * ((lambda + 2*mu) * dvx_dx + lambda * (dvy_dy + dvz_dz)))
            @. syy_f = pml_y_3d * (pml_y_3d * syy_f + dt * ((lambda + 2*mu) * dvy_dy + lambda * (dvx_dx + dvz_dz)))
            @. szz_f = pml_z_3d * (pml_z_3d * szz_f + dt * ((lambda + 2*mu) * dvz_dz + lambda * (dvx_dx + dvy_dy)))
            @. sxy_f = pml_x_3d * pml_y_3d * (sxy_f + dt * mu * (dvx_dy + dvy_dx))
            @. sxz_f = pml_x_3d * pml_z_3d * (sxz_f + dt * mu * (dvx_dz + dvz_dx))
            @. syz_f = pml_y_3d * pml_z_3d * (syz_f + dt * mu * (dvy_dz + dvz_dy))
        else
            @. sxx_f = pml_x_3d * (pml_x_3d * sxx_f + dt * ((lambda + 2*mu) * dvx_dx + lambda * (dvy_dy + dvz_dz)))
            @. syy_f = pml_y_3d * (pml_y_3d * syy_f + dt * ((lambda + 2*mu) * dvy_dy + lambda * (dvx_dx + dvz_dz)))
            @. szz_f = pml_z_3d * (pml_z_3d * szz_f + dt * ((lambda + 2*mu) * dvz_dz + lambda * (dvx_dx + dvy_dy)))
            @. sxy_f = pml_x_3d * pml_y_3d * (sxy_f + dt * mu .* (dvx_dy + dvy_dx))
            @. sxz_f = pml_x_3d * pml_z_3d * (sxz_f + dt * mu .* (dvx_dz + dvz_dx))
            @. syz_f = pml_y_3d * pml_z_3d * (syz_f + dt * mu .* (dvy_dz + dvz_dy))
        end

        # === Stress source injection ===
        if has_stress_source(source)
            s_indices = findall(source.s_mask)
            for (field_arr, src_data) in [
                (sxx_f, source.sxx), (syy_f, source.syy), (szz_f, source.szz),
                (sxy_f, source.sxy), (sxz_f, source.sxz), (syz_f, source.syz)
            ]
                if src_data !== nothing
                    sig = size(src_data, 2) == 1 ? src_data[:, 1] : src_data[:, min(t_index, size(src_data, 2))]
                    for (j, idx) in enumerate(s_indices)
                        if source.s_mode == Dirichlet
                            field_arr[idx] = sig[j]
                        else
                            field_arr[idx] += sig[j]
                        end
                    end
                end
            end
        end

        # === Record sensor data (pressure = -trace(σ)/3) ===
        p = @. -(sxx_f + syy_f + szz_f) / 3
        if !isempty(mask_indices)
            velocities = (ux=vx, uy=vy, uz=vz)
            record_sensor_data!(sensor_data, p, velocities, sensor, mask_indices,
                               t_index, Nt, dt, medium.density)
        end

        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    finalize_sensor_data!(sensor_data, Nt)
    return SimulationOutput(sensor_data)
end
