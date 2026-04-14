# Arrays

`KWaveArray` represents a multi-element transducer array with elements specified in physical coordinates (metres).
Use it when you need distributed source/sensor signals across multiple geometric elements.

## Element Types

```@docs
ArrayElement
ArcElement
BowlElement
DiscElement
RectElement
SphereElement
```

## Array

```@docs
KWaveArray
add_arc_element!
add_bowl_element!
add_disc_element!
add_rect_element!
```

## Grid Projection & Signal Distribution

```@docs
get_element_binary_mask
get_array_binary_mask
get_distributed_source_signal
combine_sensor_data
```

## Typical Workflow

```julia
# Build a linear array of arc elements
arr = KWaveArray()
for i in 1:8
    add_arc_element!(arr; position=[i*1e-3, 0.0], radius=5e-3, diameter=2e-3, focus=[i*1e-3, 5e-3])
end

# Project onto grid
mask = get_array_binary_mask(arr, kgrid)    # BitMatrix (Nx × Ny)

# Distribute a per-element source signal
source_signal = tone_burst(10e6, 5e6, 3)   # single element signal
dist_signal = get_distributed_source_signal(arr, kgrid, source_signal)
# dist_signal: Matrix (n_source_points × Nt)

# Set up source
source = KWaveSource(p_mask=mask, p=dist_signal)
```
