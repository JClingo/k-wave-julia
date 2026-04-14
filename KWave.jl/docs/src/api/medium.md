# Medium

## Acoustic Medium

```@docs
KWaveMedium
is_homogeneous
is_lossless
is_nonlinear
```

## Elastic Medium

```@docs
ElasticMedium
ElasticSource
AbsorptionParams
```

## Thermal Medium

```@docs
ThermalMedium
ThermalSource
```

## Absorption Modes

| Mode symbol | Behaviour |
|---|---|
| `:no_absorption` | Lossless (default when `alpha_coeff` not set) |
| `:no_dispersion` | Amplitude attenuation only, zero phase distortion |
| `:stokes` | Full power-law amplitude + dispersion (broadband-accurate) |

Use `:stokes` for broadband simulations; `:no_dispersion` for narrowband CW problems where phase accuracy is less critical.
