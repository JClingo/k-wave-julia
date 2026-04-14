# Unitful Integration

KWave.jl supports physical units via the [Unitful.jl](https://painterqubits.github.io/Unitful.jl/)
package extension `KWaveUnitfulExt`.

## Loading

```julia
using KWave, Unitful
```

## Usage

Pass `Unitful` quantities to `KWaveMedium`:

```julia
medium = KWaveMedium(
    sound_speed = 1500u"m/s",
    density     = 1000u"kg/m^3",
    alpha_coeff = 0.5u"dB/MHz/cm",   # dB/(MHz cm) — if supported
)
```

Values are stripped of units and converted to `Float64` (SI) internally before the simulation runs.

## Grid Spacings

Grid spacings in `KWaveGrid` currently require plain `Float64` in metres:

```julia
dx = ustrip(u"m", 0.1u"mm")   # convert: 1e-4
kgrid = KWaveGrid(256, dx, 256, dx)
```

## Limitations

- Not all `KWaveMedium` fields accept Unitful quantities — check the `KWaveUnitfulExt` source for the full list of supported conversions.
- Source and sensor arrays (`p0`, `p`, `ux`, …) must still be plain `AbstractArray{Float64}`.

## See Also

[`KWaveMedium`](@ref)
