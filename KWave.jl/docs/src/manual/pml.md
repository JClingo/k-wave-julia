# PML Configuration

## What the PML Does

The **Perfectly Matched Layer** (PML) absorbs outgoing waves at domain boundaries, preventing
reflections from wrapping around due to the periodic nature of the DFT.
Without a PML, energy leaving one side re-enters from the opposite side.

## Default Behaviour

By default, the PML sits **inside** the declared grid (`pml_inside=true`).
The effective simulation domain is therefore smaller than the declared grid:

```
effective domain = grid - 2 × pml_size (per dimension)
```

Set `pml_inside=false` to expand the grid instead, preserving the full declared domain at the
cost of increased memory and computation.

## Tuning

| Parameter | Default | Effect |
|---|---|---|
| `pml_size` | `20` | Thickness in grid points. Thicker → better absorption, more memory |
| `pml_alpha` | `2.0` | Absorption strength. Higher → stronger attenuation; too high → reflection |

Recommended values by grid size:

| Grid points per dimension | `pml_size` | `pml_alpha` |
|---|---|---|
| < 128 | 10–15 | 2.0 |
| 128–512 | 20 | 2.0 |
| > 512 | 40–60 | 2.0 |

Use [`get_optimal_pml_size`](@ref) for an automated estimate.

## PML Profile

The PML absorption profile is a 4th-order polynomial:

```
σ(x) = pml_alpha · (x / L)^4
```

where `x` is the distance from the inner PML boundary and `L = pml_size * dx`.

## See Also

[`get_pml`](@ref), [`get_optimal_pml_size`](@ref), [`kspace_first_order`](@ref)
