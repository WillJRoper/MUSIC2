# Choosing parameters

How to pick the main settings. The [wiki](https://github.com/cosmo-sims/music/wiki)
lists every option; [`ics_example.conf`](../ics_example.conf) and
[`data/ExampleConfigs/swift_zoom.conf`](../data/ExampleConfigs/swift_zoom.conf)
are working starting points.

## Resolution: `levelmin`, `levelmax`, `levelmin_TF`

A level `L` is a grid of `2^L` cells per side across the whole box.

- `levelmin`: the base grid, covering the whole box. Its particles are the
  low-resolution background. Use the parent run's level for zooms.
- `levelmax`: the finest level, inside the refinement region. The
  high-resolution particle spacing is `boxlength / 2^levelmax`.
  `levelmin = levelmax` gives uniform-resolution ICs.
- `levelmin_TF`: the level on which the density field is first computed
  before being restricted to `levelmin`. Set it 1-2 levels above `levelmin` for
  accurate large-scale power in zooms; it cannot be below `levelmin` (MUSIC
  raises it with a warning).

Particle masses follow from the box: one high-resolution particle (DM + gas)
has mass `Omega_m × ρ_crit × (boxlength / 2^levelmax)^3`. With the
[masked load](swift-output.md#masked-particle-load) the gas and DM spacings
differ; MUSIC reports the actual masses in `/ICs_parameters/Derived` of the
SWIFT file.

Memory and time grow with the size of the largest grid, i.e. the padded
refinement region at `levelmax` (and the full `2^levelmin_TF` grid). A
`Panph6` descriptor adds about 200 bytes per cell of the largest padded level.

## Zoom geometry: `padding`, `overlap`, region

- `padding` (required for zooms): extra cells of each coarser level around the
  level above, so the transition between resolutions sits away from the
  object. 8 is typical.
- `overlap` (default 4): extra cells used when solving each refined patch,
  then discarded. Leave at 4.
- `align_top = no` lets refinement levels sit anywhere on the coarser grid.
- The region itself: [Refinement regions](refinement-regions.md).

## Starting redshift and LPT order

- `zstart`: start where the largest displacements are still small compared
  with the finest particle spacing. 2LPT allows a later start than 1LPT:
  z ≈ 50-127 is common for galaxy-scale zooms with 2LPT, higher for very small
  boxes or 1LPT.
- `use_2LPT = yes` adds second-order corrections; use it unless you need 1LPT.
- `use_LLA = yes` derives the gas density from the local Lagrangian
  approximation instead of the linear field.
- With `transfer = CLASS`, `zstart` must be below recombination (about 1088);
  CLASS fails above it.

## Baryons

`baryons = yes` writes gas as well as dark matter, using separate CDM and
baryon transfer functions. For SWIFT, gas is placed half a cell from the DM by
default (`baryon_staggering = 0.5`), or by the
[masked load](swift-output.md#masked-particle-load). The SWIFT plugin switches
on `do_SPH` itself, so the gas displacement field is computed at the gas
positions.

## Cosmology

Set the parameters directly in `[cosmology]`:

```ini
[cosmology]
Omega_m = 0.3099
Omega_b = 0.048891054
Omega_L = 0.6901
H0      = 67.742
sigma_8 = 0.8102      # or A_s = 2.1064e-9
n_s     = 0.96822
```

Optional: `w_0`, `w_a`, `Tcmb`, `YHe`, `N_ur`, `m_nu1..3`, `k_p`. With
`transfer = CLASS`, `A_s` takes precedence if both `sigma_8` and `A_s` are given.

`ParameterSet = Planck2018EE+BAO+SN` (and the other Planck 2018 sets) fills
these in for the cosmology code, but other parts of MUSIC (grid statistics, the
SWIFT writer) still read `Omega_m`, `Omega_b`, `Omega_L` and `H0` directly and
stop if they are missing. Give those explicitly in either case.

## Transfer functions

| `transfer` | Notes |
|---|---|
| `CLASS` | computed at run time; needs `ENABLE_CLASS`. Normalised to `sigma_8` (or `A_s`) at `ztarget` (default 0). The CLASS parameters are written to `input_class_parameters.ini` and to the SWIFT file. |
| `camb_file` | a CAMB transfer table given by `transfer_file`; it must be for `zstart`. Example: `data/EAGLE/Eagle_camb_file.txt`. |
| `eisenstein` | Eisenstein & Hu fit; no separate baryon field, fine for tests. |
| `eisenstein_wdm`, `bbks`, `music`, `linger++` | see the wiki. |

## White noise

- `generator = MUSIC` with `seed[L]` per level: MUSIC's own noise. Keep the
  seeds when changing `levelmax` so structure is preserved.
- `generator = PANPHASIA` with a `descriptor`: noise fixed independently of
  the grid, so any code or resolution reproduces the same structure. See
  [PANPHASIA](panphasia.md).

## Files MUSIC writes

All in the directory MUSIC is run from:

| File | Contents |
|---|---|
| `output/filename` | the initial conditions |
| `<config>_log.txt` | full log, including the grid structure, the domain shift and errors |
| `input_powerspec.txt` | the input power spectrum by component (also stored in the SWIFT file) |
| `input_class_parameters.ini` | CLASS input, for `transfer = CLASS` |
| `dump_transfer.txt` | the transfer table as read, for `transfer = camb_file` |
| `griddump.dat` | the level map, for `region = multibox` |
