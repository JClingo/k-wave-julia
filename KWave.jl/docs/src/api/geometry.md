# Geometry

Functions for creating binary grid masks and Cartesian point sets.

**Coordinate convention:** binary mask functions accept 1-based integer grid indices.
Cartesian functions (`make_cart_*`) use physical coordinates in metres.

## Binary Grid Masks

These return `BitArray` (or `BitMatrix`) with `true` at grid points inside the shape.

```@docs
make_disc
make_circle
make_ball
make_sphere
make_arc
make_line
make_bowl
make_multi_arc
make_multi_bowl
make_spherical_section
```

## Cartesian Point Sets

These return `Matrix{Float64}` where columns are `[x; y]` or `[x; y; z]` coordinates.

```@docs
make_cart_circle
make_cart_sphere
make_cart_arc
make_cart_bowl
make_cart_disc
make_cart_rect
```
