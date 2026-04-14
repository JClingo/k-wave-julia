# Source & Sensor

## Source Modes

```@docs
SourceMode
```

| Mode | Effect |
|---|---|
| `Dirichlet` | Source values replace field values (hard source) |
| `Additive` | Source values added to field (soft source, default) |
| `AdditiveNoCorrection` | Additive without PML correction |

## Source

```@docs
KWaveSource
```

**Source type summary:**

| Source type | Required fields | Use case |
|---|---|---|
| Photoacoustic / IVP | `p0` | Initial pressure distribution |
| Time-varying pressure | `p_mask`, `p` | Driven transducer, CW source |
| Time-varying velocity | `u_mask`, `ux`/`uy`/`uz` | Velocity boundary condition |

`p0` and time-varying sources can coexist.

## Sensor

```@docs
KWaveSensor
is_time_reversal
```

**Recordable fields:**

| Symbol | Recorded quantity | Shape |
|---|---|---|
| `:p` | Pressure time series | `(n_sensor, Nt)` |
| `:p_max` | Maximum pressure | `(n_sensor,)` |
| `:p_min` | Minimum pressure | `(n_sensor,)` |
| `:p_rms` | RMS pressure | `(n_sensor,)` |
| `:p_final` | Final pressure field | `(Nx[, Ny[, Nz]])` |
| `:ux`, `:uy`, `:uz` | Velocity components | `(n_sensor, Nt)` |
| `:u_max` | Maximum velocity magnitude | `(n_sensor,)` |
| `:u_rms` | RMS velocity magnitude | `(n_sensor,)` |
| `:u_final` | Final velocity field | grid-shaped |
| `:I_avg` | Time-averaged intensity | `(n_sensor,)` |
| `:I_max` | Maximum intensity | `(n_sensor,)` |

## Simulation Output

```@docs
SimulationOutput
```

**Accessing results:**

```julia
output = kspace_first_order(kgrid, medium, source, sensor)

# Access a field
p = output[:p]              # Matrix (n_sensor × Nt)

# Check what was recorded
keys(output)                # e.g. [:p, :p_max]
haskey(output, :p_max)      # true/false
```
