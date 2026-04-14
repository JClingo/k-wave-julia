# Source & Sensor

## Source Types

### Photoacoustic / Initial Value Problem

Set `p0` to an initial pressure distribution matching the grid size:

```julia
p0 = zeros(Nx, Ny)
p0[Nx÷2, Ny÷2] = 5e4    # 50 kPa point source
source = KWaveSource(p0=p0)
```

### Time-Varying Pressure Source

Provide a binary mask and a time series matrix:

```julia
mask = falses(Nx, Ny)
mask[10, Ny÷2] = true              # single source point
signal = tone_burst(10e6, 5e6, 3)  # tone burst at 5 MHz
source = KWaveSource(
    p_mask = mask,
    p      = reshape(signal, 1, :), # (1 × Nt)
)
```

For multiple source points, the `p` matrix must be `(n_sources × Nt)`.

### Time-Varying Velocity Source

```julia
mask = falses(Nx, Ny)
mask[:, 1] .= true        # entire left edge
source = KWaveSource(
    u_mask = mask,
    ux     = velocity_x,  # (n_sources × Nt)
    uy     = velocity_y,
)
```

### Source Modes

| Mode | Effect |
|---|---|
| `Dirichlet` | Hard source — replaces field values at source points |
| `Additive` | Soft source — adds to existing field (default) |
| `AdditiveNoCorrection` | Additive without PML correction term |

## Sensor Configuration

### Binary Grid Mask

```julia
mask = falses(Nx, Ny)
mask[1, :]   .= true   # top row
sensor = KWaveSensor(mask=mask, record=[:p])
```

### Full-Field Recording

`mask=nothing` records the entire pressure field at every step. Memory-intensive for large grids.

```julia
sensor = KWaveSensor(record=[:p_final])  # record only at t=end
```

### Cartesian Sensor Positions

Pass a `(2 × n_points)` matrix of physical coordinates [m]:

```julia
positions = [0.0  1e-3  2e-3;    # x coordinates
             0.0  0.0   0.0 ]    # y coordinates
sensor = KWaveSensor(mask=positions, record=[:p])
```

### Time Reversal

For time-reversal reconstruction, set `time_reversal_boundary_data` to previously recorded sensor data:

```julia
sensor_tr = KWaveSensor(
    mask                        = sensor.mask,
    time_reversal_boundary_data = output[:p],  # from forward simulation
)
is_time_reversal(sensor_tr)  # true
```

## See Also

[`KWaveSource`](@ref), [`KWaveSensor`](@ref), [`SimulationOutput`](@ref),
[`SourceMode`](@ref), [`tone_burst`](@ref), [`cart2grid`](@ref)
