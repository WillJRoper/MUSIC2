# Refinement regions

When `levelmin != levelmax`, MUSIC refines a region of the box down to
`levelmax`. `[setup] / region` picks its shape: `box` (default), `ellipsoid`,
`convex_hull` or `multibox`. All positions are in box units, on the unit cube
`[0,1)^3`.

## `box`

```ini
[setup]
region     = box
ref_center = 0.5, 0.5, 0.5
ref_extent = 0.2, 0.2, 0.2
# ref_offset = 0.4, 0.4, 0.4    # instead of ref_center
# ref_dims   = 128, 128, 128    # instead of ref_extent
```

Give one position and one size:

- `ref_center`: centre of the region.
- `ref_offset`: lower corner of the region. `ref_center` wins if both are
  given.
- `ref_extent`: size as a fraction of the box.
- `ref_dims`: size in `levelmax` cells. Must be even, because the finest grid
  is refined in whole 2x2x2 blocks. `ref_extent` wins if both are given.

## `ellipsoid`

Fitted to points:

```ini
[setup]
region            = ellipsoid
region_point_file = halo_lagrangian.txt
```

- `region_point_file`: the points to enclose (see [Point files](#point-files)).
  MUSIC fits the minimum-volume ellipsoid around them and logs its matrix and
  centre, so they can be reused as below.

Or given explicitly, as `(x - c)^T A (x - c) <= 1`:

```ini
[setup]
region                     = ellipsoid
region_ellipsoid_matrix[0] = 2022.3683, -727.43207, 647.41876
region_ellipsoid_matrix[1] = -727.43207, 1231.4218, 156.45541
region_ellipsoid_matrix[2] = 647.41876, 156.45541, 1635.1184
region_ellipsoid_center    = 0.56548864, 0.5168543, 0.6127409
```

- `region_ellipsoid_matrix[0..2]`: the rows of `A`.
- `region_ellipsoid_center`: the centre `c`.

This is the form [CosmICweb](https://cosmicweb.eu) produces (see
`data/ExampleConfigs/AURIGA6_ics.conf`).

Either way the finest level's ellipsoid is grown by one finest cell, and each
coarser level's by `padding` cells of that level.

## `convex_hull`

```ini
[setup]
region            = convex_hull
region_point_file = halo_lagrangian.txt
```

- `region_point_file`: the points to enclose. The finest level is their
  convex hull grown by √3 finest cells; each coarser level adds `padding`
  cells. Points straddling the periodic boundary are unwrapped relative to
  the first point.

Fits elongated regions more tightly than `ellipsoid`.

## `multibox`

```ini
[setup]
region             = multibox
region_point_file  = halo_lagrangian.txt
region_point_level = 9
```

- `region_point_file`: the points to enclose.
- `region_point_level`: the points are binned on a grid of
  `2^(region_point_level - 1)` cells per side. Default `levelmin`.

Every grid cell holding a point is refined to `levelmax`, plus a one-cell
buffer; each coarser level adds another one-cell buffer. This follows
irregular regions without an ellipsoid's overshoot. The level map is written
to `griddump.dat`. `multibox` does not wrap points across the box edge and
ignores `region_point_shift`.

## Point files

`region_point_file` is plain text, one point per line:

- `x y z`: positions in box units.
- `x y z vx vy vz`: each position is moved back along its velocity,
  `x - v / vfact`, a Zel'dovich estimate of its Lagrangian position. Use this
  for points taken from a snapshot rather than traced back to the ICs by ID.

Points are usually a target halo's particles traced back to their initial
positions in a uniform parent run.
[`tools/lagrangian_points.py`](../tools/lagrangian_points.py) does this for
SWIFT-style files; see the [Quickstart](quickstart.md#3-zoom-on-a-halo).

### Reusing points from a previous zoom

```ini
[setup]
region_point_shift    = 4, -2, 0
region_point_levelmin = 7
```

MUSIC shifts zooms so the region sits near the box centre, logging
`Domain will be shifted by (i, j, k)`. Points taken from such a run must have
that shift undone:

- `region_point_shift`: the shift logged by the previous run.
- `region_point_levelmin`: that run's `levelmin`.

`ellipsoid` and `convex_hull` only.

## Options shared by all regions

```ini
[setup]
padding     = 8
align_top   = no
no_shift    = no
force_shift = no
```

- `padding`: cells added around each level, in that level's cells.
- `align_top`: snap every level to the coarsest grid's cells. Default `no`.
- `no_shift`: do not shift the domain to centre the region. Default `no`.
- `force_shift`: shift anyway, overriding `no_shift`. Default `no`.

Uniform runs (`levelmin = levelmax`) are never shifted.
