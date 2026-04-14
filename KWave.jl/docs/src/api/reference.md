# Reference Solutions

Analytical and semi-analytical solutions for validation against numerical simulations.

## Focused Transducer Acoustic Field

```@docs
focused_annulus_oneil
focused_bowl_oneil
```

Both functions implement O'Neil's (1949) solution for the on-axis and off-axis pressure field of
a focused piston transducer in a lossless fluid.
The return value is a complex pressure array; use `abs.()` for amplitude and `angle.()` for phase.

## Nonlinear Propagation

```@docs
mendousse
```

`mendousse` computes the nonlinear distortion of a finite-amplitude plane wave using the
Mendousse (1953) solution (valid for weak nonlinearity).

## References

- H. T. O'Neil, "Theory of focusing radiators," *J. Acoust. Soc. Am.*, 21(5), 516–526, 1949.
- J. S. Mendousse, "Nonlinear dissipative distortion of progressive sound waves at moderate amplitudes," *J. Acoust. Soc. Am.*, 25(1), 51–54, 1953.
