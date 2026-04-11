# ============================================================================
# KWave.jl — kspace_first_order: Main simulation entry point
# ============================================================================

"""
    kspace_first_order(kgrid, medium, source, sensor; kwargs...)

Run a k-space pseudospectral time-domain simulation.

This is the primary simulation function, dispatching on the grid dimensionality
to run 1D, 2D, or 3D simulations.

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
- `progress_callback`: Optional callback `f(t_index, Nt, p_field)` called each step

# Returns
`SimulationOutput` containing recorded sensor data.

# Example
```julia
# Create grid
kgrid = KWaveGrid(128, 1e-4, 128, 1e-4)
make_time!(kgrid, 1500.0)

# Define medium
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)

# Create initial pressure source (disc)
p0 = Float64.(make_disc(128, 128, 64, 64, 10))
source = KWaveSource(p0=p0)

# Record everywhere
sensor = KWaveSensor(mask=trues(128, 128))

# Run simulation
output = kspace_first_order(kgrid, medium, source, sensor)
```
"""
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
    progress_callback::Union{Nothing, Function}=nothing,
) where T <: AbstractFloat

    # ================================================================
    # 1. Input validation
    # ================================================================
    if kgrid.Nt[] == 0 || kgrid.dt[] == 0
        error("Time array not set. Call make_time!(kgrid, sound_speed) before running the simulation.")
    end
    validate_record_fields(sensor)

    Nx, Ny = kgrid.Nx, kgrid.Ny
    dt = kgrid.dt[]
    Nt = kgrid.Nt[]

    # Normalize PML parameters
    pml_x_size, pml_y_size = pml_size isa Tuple ? pml_size : (pml_size, pml_size)
    pml_x_alpha, pml_y_alpha = pml_alpha isa Tuple ? pml_alpha : (pml_alpha, pml_alpha)

    # ================================================================
    # 2. Save-to-disk mode (for C++ binary execution)
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
    # 3. Pre-compute PML absorption arrays
    # ================================================================
    pml_x = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt)
    pml_y = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt)
    pml_x_sgx = get_pml(Nx, kgrid.dx, pml_x_size, pml_x_alpha, dt; staggered=true)
    pml_y_sgy = get_pml(Ny, kgrid.dy, pml_y_size, pml_y_alpha, dt; staggered=true)

    # ================================================================
    # 4. Pre-compute k-space correction factor
    # ================================================================
    # kappa = sinc(c_ref * k * dt / 2)
    # where sinc(x) = sin(πx)/(πx) in Julia (normalized sinc)
    # k-Wave uses unnormalized sinc: sin(x)/x
    # So: kappa = sin(c_ref * k * dt / 2) / (c_ref * k * dt / 2)
    c_ref = medium.sound_speed isa Real ? medium.sound_speed : mean(medium.sound_speed)
    kappa = similar(kgrid.k)
    for j in 1:Ny
        for i in 1:Nx
            arg = c_ref * kgrid.k[i, j] * dt / 2
            kappa[i, j] = arg ≈ 0 ? 1.0 : sin(arg) / arg
        end
    end

    # ================================================================
    # 5. Create FFT plans
    # ================================================================
    plans = create_fft_plans(kgrid; data_cast=T)

    # ================================================================
    # 6. Allocate field arrays
    # ================================================================
    p = zeros(Float64, Nx, Ny)
    ux = zeros(Float64, Nx, Ny)
    uy = zeros(Float64, Nx, Ny)
    rhox = zeros(Float64, Nx, Ny)
    rhoy = zeros(Float64, Nx, Ny)
    scratch1 = zeros(ComplexF64, Nx, Ny)
    scratch2 = zeros(ComplexF64, Nx, Ny)

    # ================================================================
    # 7. Allocate sensor recording storage
    # ================================================================
    sensor_data = _create_sensor_data(sensor, (Nx, Ny), Nt)
    mask_indices = if sensor.mask isa AbstractArray{Bool}
        findall(sensor.mask)
    else
        CartesianIndex{2}[]
    end

    # ================================================================
    # 8. Initialize from p0 (photoacoustic initial value problem)
    # ================================================================
    initialize_p0_2d!(p, ux, uy, rhox, rhoy, source, kgrid, medium, plans, scratch1, smooth_p0)

    # ================================================================
    # 9. Time loop
    # ================================================================
    prog = Progress(Nt; desc="k-Wave simulation: ", enabled=!isnothing(progress_callback) || true)

    for t_index in 1:Nt
        # Execute one time step
        time_step_2d!(
            p, ux, uy, rhox, rhoy,
            scratch1, scratch2,
            kgrid, medium, source,
            pml_x, pml_y, pml_x_sgx, pml_y_sgy,
            kappa, plans, t_index,
        )

        # Record sensor data
        if !isempty(mask_indices)
            record_sensor_data!(sensor_data, p, ux, uy, sensor, mask_indices, t_index, Nt)
        end

        # Progress callback
        if progress_callback !== nothing
            progress_callback(t_index, Nt, p)
        end

        next!(prog)
    end

    # ================================================================
    # 10. Finalize and return
    # ================================================================
    finalize_sensor_data!(sensor_data, Nt)

    return SimulationOutput(sensor_data)
end

# ============================================================================
# 1D solver stub (TODO: implement in Phase 2)
# ============================================================================

function kspace_first_order(
    kgrid::KWaveGrid1D,
    medium::KWaveMedium,
    source::KWaveSource,
    sensor::KWaveSensor;
    kwargs...
)
    error("1D solver not yet implemented. Coming in Phase 2.")
end

# ============================================================================
# 3D solver stub (TODO: implement in Phase 2)
# ============================================================================

function kspace_first_order(
    kgrid::KWaveGrid3D,
    medium::KWaveMedium,
    source::KWaveSource,
    sensor::KWaveSensor;
    kwargs...
)
    error("3D solver not yet implemented. Coming in Phase 2.")
end
