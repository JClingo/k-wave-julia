# Reconstruction

## k-Space Reconstruction

```@docs
kspace_line_recon
kspace_plane_recon
```

## Ultrasound Beamforming

```@docs
beamform_delay_and_sum
scan_conversion
```

## Notes

- `kspace_line_recon` and `kspace_plane_recon` assume the sensor is a flat line/plane perpendicular to the imaging axis and that the medium is homogeneous.
- `beamform_delay_and_sum` operates on recorded sensor data after simulation; for real-time simulation display use the `progress_callback` in [`kspace_first_order`](@ref).
- `scan_conversion` converts a delay-and-sum image in polar coordinates (radius, angle) to a Cartesian image suitable for display.
