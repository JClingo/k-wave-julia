# Transducer

`KWaveTransducer` models a 1D linear array transducer aligned along the y-axis.
It supports active element selection, beamforming delays, and combined source/sensor operation.

```@docs
KWaveTransducer
get_transducer_binary_mask
get_transducer_source
combine_transducer_sensor_data
```

## Notes

- Elements are indexed along the y-axis.
- Use the `active_elements` field to fire a subset (partial aperture, synthetic aperture).
- `get_transducer_source` returns a pre-multiplied source signal with the correct delays applied.
- `combine_transducer_sensor_data` sums weighted per-element receive signals to form an RF line.

## See Also

[`KWaveArray`](@ref) for arbitrary-geometry arrays.
[`beamform_delay_and_sum`](@ref) for post-hoc DAS beamforming on recorded sensor data.
