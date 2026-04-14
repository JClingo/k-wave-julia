# Theory

## k-Space Pseudospectral Method

KWave.jl solves the acoustic wave equations using the **k-space pseudospectral time-domain (PSTD)** method described by Treeby & Cox (2010). Spatial derivatives are evaluated exactly in the Fourier domain; time stepping uses a first-order staggered scheme with a k-space correction factor `κ` that eliminates numerical dispersion.

### Wave Equations

For a linear, lossless, homogeneous medium the coupled first-order equations are:

```
∂ρ/∂t = -ρ₀ ∇·u
ρ₀ ∂u/∂t = -∇p
∂p/∂t = -ρ₀ c₀² ∇·u
```

where `p` is acoustic pressure, `u` is particle velocity, `ρ` is acoustic density, `ρ₀` is the ambient density, and `c₀` is the sound speed.

### Spatial Derivatives

Spatial derivatives are computed as:

```
∂f/∂x ≈ F⁻¹{ iκₓ · F{f} }
```

where `F` is the discrete Fourier transform and `κₓ` is the wavenumber vector. For staggered grids the shift operators `exp(±iκₓ dx/2)` are applied before the inverse transform.

### k-Space Correction

The k-space correction factor

```
κ = sinc(c_ref k dt / 2)
```

is applied to the wavenumber operator before each time step. This compensates for the phase error introduced by the finite-difference time stepping scheme and effectively gives **spectral accuracy in both space and time** for homogeneous media.

### Absorption Model

Power-law absorption is implemented following Treeby & Cox (2010) using two fractional Laplacian operators:

```
absorb_nabla1 = τ · k^y       # amplitude absorption
absorb_nabla2 = η · k^(y-1)   # phase dispersion
```

where `τ = -2α c_ref^(y-1)`, `η = 2α c_ref^(y-1) tan(πy/2)`, `α` is the absorption coefficient in Np/m, and `y` is the power-law exponent.

Three absorption modes are available:

| Mode | Behaviour |
|---|---|
| `:no_absorption` | Lossless (default if `alpha_coeff` not set) |
| `:no_dispersion` | Amplitude loss only, no phase distortion |
| `:stokes` | Full power-law absorption + dispersion (Stokes-like, broadband accurate) |

### Nonlinearity

The cumulative nonlinear term (Westervelt equation) is included when `medium.BonA` is set:

```
∂p/∂t += (B/A) / (2 ρ₀ c₀²) · p · ∂p/∂t
```

### PML

The Perfectly Matched Layer (PML) absorbs outgoing waves at domain boundaries. It is implemented as a frequency-shifted damping term applied as a multiplication in the Fourier domain each time step. The default profile is a 4th-order polynomial:

```
σ(x) = pml_alpha · (x / L)^4
```

where `L` is the PML thickness and `x` is the distance from the inner edge. See [PML Configuration](manual/pml.md) for tuning guidance.

## References

- B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave fields," *J. Biomed. Opt.*, 15(2), 021314, 2010.
- B. E. Treeby and B. T. Cox, "Modeling power law absorption and dispersion for acoustic propagation using the fractional Laplacian," *J. Acoust. Soc. Am.*, 127(5), 2741–2748, 2010.
- E. S. Wise, B. T. Cox, J. Jaros, and B. E. Treeby, "Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods," *IEEE T. Ultrason. Ferroelectr. Freq. Control*, 66(3), 506–521, 2019.
