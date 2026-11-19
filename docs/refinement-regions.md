# Refinement regions

When `levelmin != levelmax`, MUSIC refines a region of the box down to
`levelmax`. The shape of that region is chosen in `[setup]`:

```ini
[setup]
region = box          # default; or ellipsoid, convex_hull, multibox
```

All positions are in box units, i.e. on the unit cube `[0,1)^3`.

## `box` (default)

A rectangular region, set by one size option and one position option. Both
are required.

| Option | Meaning |
|---|---|
| `ref_extent = lx, ly, lz` | size as a fraction of the box |
| `ref_dims = nx, ny, nz` | size in finest-level (`levelmax`) cells; must be even |
| `ref_center = x, y, z` | centre of the region |
| `ref_offset = x, y, z` | lower corner of the region |

If both size options are given, `ref_extent` wins; likewise `ref_center` wins
over `ref_offset`. Odd `ref_dims` are rejected, because the finest grid is
refined and masked in whole 2x2x2 blocks.

```ini
ref_center = 0.5, 0.5, 0.5
ref_extent = 0.2, 0.2, 0.2
```

## `ellipsoid`

An ellipsoid, either fitted to a set of points or given explicitly.

**From points** (minimum-volume bounding ellipsoid):

```ini
region            = ellipsoid
region_point_file = halo_lagrangian.txt
```

**Explicitly**, as the matrix `A` and centre `c` of the ellipsoid
`(x - c)^T A (x - c) <= 1`:

```ini
region                     = ellipsoid
region_ellipsoid_matrix[0] = 2022.3683, -727.43207, 647.41876
region_ellipsoid_matrix[1] = -727.43207, 1231.4218, 156.45541
region_ellipsoid_matrix[2] = 647.41876, 156.45541, 1635.1184
region_ellipsoid_center    = 0.56548864, 0.5168543, 0.6127409
```

This is the form [CosmICweb](https://cosmicweb.eu) produces (see
`data/ExampleConfigs/AURIGA6_ics.conf`). When fitting from points, MUSIC logs
the fitted `A` and `c`, so they can be pasted into later configs.

The ellipsoid is grown by one finest cell, and each coarser level's ellipsoid
by `padding` cells of that level.

## `convex_hull`

The convex hull of the points in `region_point_file`, grown by √3 finest
cells. Coarser levels are grown by a further `padding` cells per level. The
region centre is the hull centroid. Points that straddle the periodic boundary
are unwrapped relative to the first point.

```ini
region            = convex_hull
region_point_file = halo_lagrangian.txt
```

## `multibox`

A blocky region built on a grid: every grid cell holding a point is refined to
`levelmax`, plus a one-cell buffer around it, and each coarser level adds a
further one-cell buffer around the level above. This follows irregular
Lagrangian regions more closely than a single box without an ellipsoid's
overshoot.

```ini
region             = multibox
region_point_file  = halo_lagrangian.txt
region_point_level = 9        # optional, default levelmin
```

The grid has `2^(region_point_level - 1)` cells per side. The resulting level
map is written to `griddump.dat`. Unlike the other point-based regions,
`multibox` does not handle points near the box edge (no periodic wrapping) and
ignores `region_point_shift`.

## Point files

`region_point_file` is a plain text file, one point per line:

- **3 columns** `x y z`: positions in box units.
- **6 columns** `x y z vx vy vz`: each position is moved back along its
  velocity, `x - v / vfact`, giving a Zel'dovich estimate of its Lagrangian
  position. Use this when the points come from a snapshot rather than from
  particle IDs mapped back to the ICs.

The usual source is the particles of a target halo traced back to their
initial (Lagrangian) positions in a uniform-resolution parent run.
[`tools/lagrangian_points.py`](../tools/lagrangian_points.py) does this for
SWIFT-style files; see the [Quickstart](quickstart.md#3-zoom-on-a-halo).

### Re-using points from a previous zoom

For zooms MUSIC shifts the whole volume so the refinement region sits near the
centre, and logs `Domain will be shifted by (i, j, k)`. If the points come from
a run that was itself shifted, undo that shift with

```ini
region_point_shift    = i, j, k   # the shift logged by the previous run
region_point_levelmin = 7         # levelmin of the previous run
```

(`ellipsoid` and `convex_hull` only.)

## Options shared by all regions

| Option | Meaning |
|---|---|
| `padding` | cells of padding added around each level, in that level's cells |
| `align_top` | align every level with the coarsest grid cells |
| `no_shift` / `force_shift` | disable / force the recentring shift |
