# Medium

## Acoustic Medium

`KWaveMedium` holds the acoustic properties of the propagation medium.

### Homogeneous Medium

Scalar sound speed and density produce a **homogeneous** (spatially uniform) medium:

```julia
medium = KWaveMedium(sound_speed=1500.0, density=1000.0)
is_homogeneous(medium)  # true
```

### Heterogeneous Medium

Pass arrays (matching grid size) for spatially varying properties:

```julia
c = 1500.0 .+ 200.0 .* rand(Nx, Ny)   # random variation ±200 m/s
medium = KWaveMedium(sound_speed=c, density=1000.0)
is_homogeneous(medium)  # false
```

### Absorbing Medium

Set `alpha_coeff` [dB/(MHz^y cm)] and `alpha_power` to enable power-law absorption:

```julia
medium = KWaveMedium(
    sound_speed = 1540.0,
    density     = 1050.0,
    alpha_coeff = 0.5,    # liver-like tissue
    alpha_power = 1.3,
    alpha_mode  = :stokes,
)
is_lossless(medium)  # false
```

### Nonlinear Medium

Set `BonA` (the nonlinearity parameter B/A) to enable cumulative nonlinear effects:

```julia
medium = KWaveMedium(
    sound_speed = 1500.0,
    density     = 1000.0,
    BonA        = 3.5,    # water at 20°C
)
is_nonlinear(medium)  # true
```

### Absorption Modes

| Mode | When to use |
|---|---|
| `:no_absorption` | Lossless simulation (fastest) |
| `:no_dispersion` | Narrowband CW; removes phase dispersion artefacts |
| `:stokes` | Broadband / pulsed simulations; accurate amplitude + phase |

## Elastic Medium

For elastic wave problems use `ElasticMedium` with `pstd_elastic_2d` / `pstd_elastic_3d`:

```julia
elastic_medium = ElasticMedium(
    sound_speed_compression = 2500.0,   # P-wave speed [m/s]
    sound_speed_shear       = 1200.0,   # S-wave speed [m/s]
    density                 = 2200.0,   # [kg/m³]
)
```

## Thermal Medium

For bioheat simulations use `ThermalMedium` with `kwave_diffusion`:

```julia
thermal_medium = ThermalMedium(
    thermal_conductivity  = 0.52,    # [W/(m·K)]
    specific_heat         = 3600.0,  # [J/(kg·K)]
    density               = 1050.0, # [kg/m³]
    perfusion_coefficient = 0.0,    # blood perfusion [1/s]
)
```

## Unitful Integration

With `using Unitful`, quantities with units are accepted:

```julia
using KWave, Unitful
medium = KWaveMedium(sound_speed=1500u"m/s", density=1000u"kg/m^3")
```

Values are converted to `Float64` internally; grid spacings still require plain `Float64`.

## See Also

[`KWaveMedium`](@ref), [`ElasticMedium`](@ref), [`ThermalMedium`](@ref),
[`is_homogeneous`](@ref), [`is_lossless`](@ref), [`is_nonlinear`](@ref),
[`db2neper`](@ref), [`water_sound_speed`](@ref)
