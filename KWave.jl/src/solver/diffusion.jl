# ============================================================================
# KWave.jl — Thermal simulation solvers
# ============================================================================
# kwave_diffusion: Pennes' bioheat equation solver (1D/2D/3D)
# bioheat_exact: Analytical solution for homogeneous bioheat equation
# ============================================================================

"""
    ThermalMedium{T}

Thermal medium properties for bioheat simulations.

# Fields
- `thermal_conductivity`: Thermal conductivity [W/(m·K)]
- `density`: Mass density [kg/m³]
- `specific_heat`: Specific heat capacity [J/(kg·K)]
- `perfusion_rate`: Blood perfusion rate [kg/(m³·s)] (optional)
- `blood_temperature`: Arterial blood temperature [°C] (default: 37.0)
- `blood_specific_heat`: Blood specific heat [J/(kg·K)] (default: 3617.0)
- `metabolic_rate`: Metabolic heat generation [W/m³] (optional)
"""
struct ThermalMedium{T<:AbstractFloat}
    thermal_conductivity::Union{T, AbstractArray{T}}
    density::Union{T, AbstractArray{T}}
    specific_heat::Union{T, AbstractArray{T}}
    perfusion_rate::Union{Nothing, T, AbstractArray{T}}
    blood_temperature::T
    blood_specific_heat::T
    metabolic_rate::Union{Nothing, T, AbstractArray{T}}
end

function ThermalMedium(;
    thermal_conductivity::Union{Real, AbstractArray{<:Real}},
    density::Union{Real, AbstractArray{<:Real}},
    specific_heat::Union{Real, AbstractArray{<:Real}},
    perfusion_rate::Union{Nothing, Real, AbstractArray{<:Real}}=nothing,
    blood_temperature::Real=37.0,
    blood_specific_heat::Real=3617.0,
    metabolic_rate::Union{Nothing, Real, AbstractArray{<:Real}}=nothing)

    _f64 = x -> x === nothing ? nothing : (x isa AbstractArray ? Float64.(x) : Float64(x))
    return ThermalMedium{Float64}(
        _f64(thermal_conductivity), _f64(density), _f64(specific_heat),
        _f64(perfusion_rate), Float64(blood_temperature),
        Float64(blood_specific_heat), _f64(metabolic_rate),
    )
end

"""
    ThermalSource{T}

Heat source for thermal simulations.

# Fields
- `Q`: Volumetric heat source [W/m³] — array matching grid, or matrix (grid_points × Nt)
- `mode`: `:steady` or `:time_varying`
"""
struct ThermalSource{T<:AbstractFloat}
    Q::Union{Nothing, AbstractArray{T}}
    mode::Symbol  # :steady or :time_varying
end

function ThermalSource(;
    Q::Union{Nothing, AbstractArray{<:Real}}=nothing,
    mode::Symbol=:steady)
    _f64 = x -> x === nothing ? nothing : Float64.(x)
    return ThermalSource{Float64}(_f64(Q), mode)
end

# ============================================================================
# kwave_diffusion: Pennes' bioheat equation solver
# ============================================================================

"""
    kwave_diffusion(kgrid, medium, source, T0, Nt; kwargs...)

Solve the Pennes' bioheat equation using a pseudospectral method.

The bioheat equation:
```
ρ·c_p ∂T/∂t = κ∇²T - w_b·c_b·(T - T_a) + Q_m + Q_ext
```

where:
- `κ`: thermal conductivity [W/(m·K)]
- `w_b`: blood perfusion rate [kg/(m³·s)]
- `c_b`: blood specific heat [J/(kg·K)]
- `T_a`: arterial blood temperature [°C]
- `Q_m`: metabolic heat [W/m³]
- `Q_ext`: external heat source [W/m³]

# Arguments
- `kgrid`: KWaveGrid (1D, 2D, or 3D)
- `medium`: ThermalMedium with tissue properties
- `source`: ThermalSource with heat deposition
- `T0`: Initial temperature distribution [°C] (scalar or array)
- `Nt`: Number of time steps

# Keyword Arguments
- `record_every`: Record temperature every N steps (default: 1)

# Returns
`(T_final, T_history)` — final temperature field and history array.
"""
function kwave_diffusion(
    kgrid::KWaveGrid1D,
    medium::ThermalMedium,
    source::ThermalSource,
    T0::Union{Real, AbstractVector},
    Nt::Int;
    record_every::Int=1,
)
    Nx = kgrid.Nx
    dt = kgrid.dt[]

    T_field = T0 isa Real ? fill(Float64(T0), Nx) : Float64.(copy(T0))

    n_records = div(Nt, record_every) + 1
    T_history = zeros(Float64, Nx, n_records)
    T_history[:, 1] = T_field
    rec_idx = 2

    kappa = medium.thermal_conductivity
    rho_cp = medium.density isa Real && medium.specific_heat isa Real ?
        medium.density * medium.specific_heat :
        medium.density .* medium.specific_heat

    plans = create_fft_plans(kgrid)
    scratch = zeros(ComplexF64, Nx)

    prog = Progress(Nt; desc="Bioheat 1D: ", enabled=true)

    for t in 1:Nt
        # Laplacian via spectral method: ∇²T = ifft(-k² * fft(T))
        T_hat = plans.forward * complex.(T_field)
        k_sq = kgrid.kx_vec.^2
        laplacian = real.(plans.inverse * (-k_sq .* T_hat))

        # Diffusion term
        if kappa isa Real
            dTdt = kappa ./ rho_cp .* laplacian
        else
            dTdt = kappa ./ rho_cp .* laplacian
        end

        # Perfusion term: -w_b * c_b * (T - T_a) / (ρ * c_p)
        if medium.perfusion_rate !== nothing
            wb_cb = medium.perfusion_rate isa Real ?
                medium.perfusion_rate * medium.blood_specific_heat :
                medium.perfusion_rate .* medium.blood_specific_heat
            @. dTdt -= wb_cb / rho_cp * (T_field - medium.blood_temperature)
        end

        # Metabolic heat: Q_m / (ρ * c_p)
        if medium.metabolic_rate !== nothing
            @. dTdt += medium.metabolic_rate / rho_cp
        end

        # External heat source
        if source.Q !== nothing
            if source.mode == :time_varying && ndims(source.Q) == 2
                Q_t = source.Q[:, min(t, size(source.Q, 2))]
            else
                Q_t = source.Q
            end
            @. dTdt += Q_t / rho_cp
        end

        # Forward Euler time step
        @. T_field += dt * dTdt

        if t % record_every == 0 && rec_idx <= n_records
            T_history[:, rec_idx] = T_field
            rec_idx += 1
        end

        next!(prog)
    end

    return T_field, T_history[:, 1:min(rec_idx-1, n_records)]
end

function kwave_diffusion(
    kgrid::KWaveGrid2D,
    medium::ThermalMedium,
    source::ThermalSource,
    T0::Union{Real, AbstractMatrix},
    Nt::Int;
    record_every::Int=1,
)
    Nx, Ny = kgrid.Nx, kgrid.Ny
    dt = kgrid.dt[]

    T_field = T0 isa Real ? fill(Float64(T0), Nx, Ny) : Float64.(copy(T0))

    n_records = div(Nt, record_every) + 1
    T_history = zeros(Float64, Nx, Ny, n_records)
    T_history[:, :, 1] = T_field
    rec_idx = 2

    kappa = medium.thermal_conductivity
    rho_cp = medium.density isa Real && medium.specific_heat isa Real ?
        medium.density * medium.specific_heat :
        medium.density .* medium.specific_heat

    plans = create_fft_plans(kgrid)

    prog = Progress(Nt; desc="Bioheat 2D: ", enabled=true)

    kx_2d = kgrid.kx_vec
    ky_2d = reshape(kgrid.ky_vec, 1, :)
    k_sq = kx_2d.^2 .+ ky_2d.^2

    for t in 1:Nt
        T_hat = plans.forward * complex.(T_field)
        laplacian = real.(plans.inverse * (-k_sq .* T_hat))

        if kappa isa Real
            dTdt = kappa ./ rho_cp .* laplacian
        else
            dTdt = kappa ./ rho_cp .* laplacian
        end

        if medium.perfusion_rate !== nothing
            wb_cb = medium.perfusion_rate isa Real ?
                medium.perfusion_rate * medium.blood_specific_heat :
                medium.perfusion_rate .* medium.blood_specific_heat
            @. dTdt -= wb_cb / rho_cp * (T_field - medium.blood_temperature)
        end

        if medium.metabolic_rate !== nothing
            @. dTdt += medium.metabolic_rate / rho_cp
        end

        if source.Q !== nothing
            if source.mode == :time_varying && ndims(source.Q) == 3
                Q_t = source.Q[:, :, min(t, size(source.Q, 3))]
            else
                Q_t = source.Q
            end
            @. dTdt += Q_t / rho_cp
        end

        @. T_field += dt * dTdt

        if t % record_every == 0 && rec_idx <= n_records
            T_history[:, :, rec_idx] = T_field
            rec_idx += 1
        end

        next!(prog)
    end

    return T_field, T_history[:, :, 1:min(rec_idx-1, n_records)]
end

function kwave_diffusion(
    kgrid::KWaveGrid3D,
    medium::ThermalMedium,
    source::ThermalSource,
    T0::Union{Real, AbstractArray{<:Real, 3}},
    Nt::Int;
    record_every::Int=1,
)
    Nx, Ny, Nz = kgrid.Nx, kgrid.Ny, kgrid.Nz
    dt = kgrid.dt[]

    T_field = T0 isa Real ? fill(Float64(T0), Nx, Ny, Nz) : Float64.(copy(T0))

    n_records = div(Nt, record_every) + 1
    T_history = zeros(Float64, Nx, Ny, Nz, n_records)
    T_history[:, :, :, 1] = T_field
    rec_idx = 2

    kappa = medium.thermal_conductivity
    rho_cp = medium.density isa Real && medium.specific_heat isa Real ?
        medium.density * medium.specific_heat :
        medium.density .* medium.specific_heat

    plans = create_fft_plans(kgrid)

    kx_3d = kgrid.kx_vec
    ky_3d = reshape(kgrid.ky_vec, 1, :, 1)
    kz_3d = reshape(kgrid.kz_vec, 1, 1, :)
    k_sq = kx_3d.^2 .+ ky_3d.^2 .+ kz_3d.^2

    prog = Progress(Nt; desc="Bioheat 3D: ", enabled=true)

    for t in 1:Nt
        T_hat = plans.forward * complex.(T_field)
        laplacian = real.(plans.inverse * (-k_sq .* T_hat))

        if kappa isa Real
            dTdt = kappa ./ rho_cp .* laplacian
        else
            dTdt = kappa ./ rho_cp .* laplacian
        end

        if medium.perfusion_rate !== nothing
            wb_cb = medium.perfusion_rate isa Real ?
                medium.perfusion_rate * medium.blood_specific_heat :
                medium.perfusion_rate .* medium.blood_specific_heat
            @. dTdt -= wb_cb / rho_cp * (T_field - medium.blood_temperature)
        end

        if medium.metabolic_rate !== nothing
            @. dTdt += medium.metabolic_rate / rho_cp
        end

        if source.Q !== nothing
            if source.mode == :time_varying && ndims(source.Q) == 4
                Q_t = source.Q[:, :, :, min(t, size(source.Q, 4))]
            else
                Q_t = source.Q
            end
            @. dTdt += Q_t / rho_cp
        end

        @. T_field += dt * dTdt

        if t % record_every == 0 && rec_idx <= n_records
            T_history[:, :, :, rec_idx] = T_field
            rec_idx += 1
        end

        next!(prog)
    end

    return T_field, T_history[:, :, :, 1:min(rec_idx-1, n_records)]
end

# ============================================================================
# bioheat_exact: Analytical solution for homogeneous bioheat
# ============================================================================

"""
    bioheat_exact(T0, Q, medium, t; x=nothing, kgrid=nothing)

Compute the exact solution to the Pennes' bioheat equation for
homogeneous media with a spatially uniform initial condition.

For an initial Gaussian temperature distribution T0(r) centered at the origin,
the exact solution uses the Green's function approach.

# Simple case (uniform initial temperature + uniform source):
```
T(t) = T_a + (T0 - T_a) * exp(-P*t) + (Q + Q_m) / P * (1 - exp(-P*t))
```
where `P = w_b * c_b / (ρ * c_p)` is the perfusion coefficient.

# Arguments
- `T0`: Initial temperature [°C] (scalar for uniform, or array for distribution)
- `Q`: External heat source [W/m³] (scalar)
- `medium`: ThermalMedium (homogeneous properties only)
- `t`: Time point(s) [s] — scalar or vector

# Returns
Temperature at time t.
"""
function bioheat_exact(T0::Real, Q::Real, medium::ThermalMedium, t::Union{Real, AbstractVector})
    kappa = medium.thermal_conductivity isa Real ? medium.thermal_conductivity :
        error("bioheat_exact requires homogeneous medium (scalar properties)")
    rho = medium.density isa Real ? medium.density :
        error("bioheat_exact requires homogeneous medium")
    cp = medium.specific_heat isa Real ? medium.specific_heat :
        error("bioheat_exact requires homogeneous medium")

    T_a = medium.blood_temperature
    wb = medium.perfusion_rate !== nothing ? medium.perfusion_rate : 0.0
    wb isa Real || error("bioheat_exact requires homogeneous medium")
    cb = medium.blood_specific_heat

    Qm = medium.metabolic_rate !== nothing ? medium.metabolic_rate : 0.0
    Qm isa Real || error("bioheat_exact requires homogeneous medium")

    rho_cp = rho * cp

    if wb > 0
        P = wb * cb / rho_cp
        Q_total = (Q + Qm) / rho_cp

        if t isa Real
            return T_a + (Float64(T0) - T_a) * exp(-P * t) + Q_total / P * (1 - exp(-P * t))
        else
            return [T_a + (Float64(T0) - T_a) * exp(-P * ti) + Q_total / P * (1 - exp(-P * ti)) for ti in t]
        end
    else
        # No perfusion: pure diffusion with source
        Q_total = (Q + Qm) / rho_cp

        if t isa Real
            return Float64(T0) + Q_total * t
        else
            return [Float64(T0) + Q_total * ti for ti in t]
        end
    end
end
