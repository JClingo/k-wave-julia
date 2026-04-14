# API Reference

Complete reference for all exported symbols in KWave.jl, organized by category.

| Category | Contents |
|---|---|
| [Solvers](solvers.md) | `kspace_first_order`, axisymmetric, CW, elastic, thermal solvers |
| [Grid](grid.md) | `KWaveGrid`, `make_time!`, wavenumber helpers |
| [Medium](medium.md) | `KWaveMedium`, `ElasticMedium`, `ThermalMedium` |
| [Source & Sensor](source_sensor.md) | `KWaveSource`, `KWaveSensor`, `SimulationOutput`, `SourceMode` |
| [Signal Processing](signal.md) | `tone_burst`, `gaussian_pulse`, `smooth`, `filter_time_series`, `spect` |
| [Geometry](geometry.md) | `make_disc`, `make_ball`, Cartesian point sets |
| [Arrays](array.md) | `KWaveArray`, element types, `get_array_binary_mask` |
| [Transducer](transducer.md) | `KWaveTransducer`, source/sensor helpers |
| [Reconstruction](reconstruction.md) | `kspace_line_recon`, `beamform_delay_and_sum` |
| [Visualization](visualization.md) | `voxel_plot`, `max_intensity_projection`, display/movie |
| [Utilities](utils.md) | `cart2grid`, `db2neper`, `water_sound_speed`, PML helpers |
| [I/O](io.md) | HDF5 read/write for k-Wave C++ binary compatibility |
| [Reference Solutions](reference.md) | Analytical solutions for validation |
| [Python Interop](interop.md) | `python_kwave_grid`, bridge to Python k-Wave |
