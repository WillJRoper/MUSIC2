# Choosing parameters

How to pick the main settings. The [wiki](https://github.com/cosmo-sims/music/wiki)
lists every option; [`ics_example.conf`](../ics_example.conf) and
[`data/ExampleConfigs/swift_zoom.conf`](../data/ExampleConfigs/swift_zoom.conf)
are working starting points.

## Resolution

```ini
[setup]
boxlength   = 100   # Mpc/h
levelmin    = 7
levelmin_TF = 8
levelmax    = 10
```

A level `L` is a grid of `2^L` cells per side across the whole box.

- `boxlength`: side of the periodic box in Mpc/h.
- `levelmin`: the base grid covering the whole box; its particles are the
  low-resolution background. For a zoom, use the parent run's level.
- `levelmin_TF`: the level the density field is first computed on, before
  being restricted to `levelmin`. 1-2 levels above `levelmin` gives accurate
  large-scale power in zooms. Default `levelmin`; MUSIC raises it to
  `levelmin` with a warning if set lower.
- `levelmax`: the finest level, inside the refinement region. The
  high-resolution particle spacing is `boxlength / 2^levelmax`.
  `levelmax = levelmin` gives uniform-resolution ICs.

One high-resolution particle (DM + gas) has mass
`Omega_m × ρ_crit × (boxlength / 2^levelmax)^3`. With the
[masked load](swift-output.md#masked-particle-load) gas and DM spacings
differ; the SWIFT file records the actual masses in `/ICs_parameters/Derived`.

Memory and time scale with the largest grid: the padded refinement region at
`levelmax`, or the full `2^levelmin_TF` grid. A `Panph6` descriptor adds
about 200 bytes per cell of the largest padded level.

## Zoom geometry

```ini
[setup]
region    = box
padding   = 8
overlap   = 4
align_top = no
```

- `region`: the shape of the refinement region; see
  [Refinement regions](refinement-regions.md).
- `padding`: cells of each coarser level added around the level above, so
  resolution transitions sit away from the object. Required for zooms; 8 is
  typical.
- `overlap`: extra cells used when solving each refined patch, then
  discarded. Default 4; leave it.
- `align_top`: `yes` snaps every level to the coarsest grid's cells; `no`
  (default) lets levels sit anywhere on the grid above.

## Starting redshift and LPT order

```ini
[setup]
zstart   = 127
use_2LPT = yes
use_LLA  = no
```

- `zstart`: start while the largest displacements are still small compared
  with the finest particle spacing. With 2LPT, z ≈ 50-127 is common for
  galaxy-scale zooms; go higher for very small boxes or 1LPT. With
  `transfer = CLASS` it must be below recombination (about 1088).
- `use_2LPT`: add second-order corrections. Default `no`; use `yes` unless
  you need 1LPT.
- `use_LLA`: derive the gas density from the local Lagrangian approximation
  instead of the linear field. Default `no`.

## Baryons

```ini
[setup]
baryons           = yes
baryon_staggering = 0.5
```

- `baryons`: write gas as well as dark matter, using separate CDM and baryon
  transfer functions. Default `no`.
- `baryon_staggering`: offset of the gas from the DM, in cells. Default 0.5.
  Ignored with the SWIFT [masked load](swift-output.md#masked-particle-load),
  which places gas itself.

Gas-writing output plugins, including `swift`, switch on `do_SPH`
themselves, so the gas displacement field is computed at the gas positions.

## Cosmology

Give the values individually:

```ini
[cosmology]
Omega_m = 0.3099
Omega_b = 0.048891054
Omega_L = 0.6901
H0      = 67.742
sigma_8 = 0.8102
n_s     = 0.96822
A_s     = 2.1064e-9
```

- `Omega_m`, `Omega_b`, `Omega_L`: matter, baryon and dark-energy densities.
- `H0`: Hubble constant in km/s/Mpc.
- `sigma_8`: normalisation. Required by every transfer function except CLASS.
- `n_s`: spectral index.
- `A_s`: primordial amplitude, used only by CLASS; it overrides `sigma_8`
  there when both are given.

Optional: `w_0`, `w_a`, `Tcmb`, `YHe`, `N_ur`, `m_nu1..3`, `k_p`.

Or use a predefined set:

```ini
[cosmology]
ParameterSet = Planck2018EE+BAO+SN
transfer     = CLASS
sigma_8      = 0.8102
```

- `ParameterSet`: one of `Planck2018EE`, `Planck2018EE+BAO`,
  `Planck2018EE+SN`, `Planck2018EE+BAO+SN`. It replaces every individual value
  above.
- `sigma_8`: sets give only `A_s`, so add `sigma_8` for any transfer function
  other than CLASS.

## Transfer functions

```ini
[cosmology]
transfer      = CLASS
ztarget       = 0
transfer_file = data/EAGLE/Eagle_camb_file.txt   # camb_file only
```

`transfer` is one of:

- `CLASS`: computed at run time; needs `ENABLE_CLASS`. Normalised at
  `ztarget` (default 0). Its input is written to `input_class_parameters.ini`
  and to the SWIFT file.
- `camb_file`: the CAMB table in `transfer_file`, which must be for `zstart`.
- `eisenstein`: Eisenstein & Hu fit; no separate baryon field, fine for
  tests.
- `eisenstein_wdm`, `bbks`, `music`, `linger++`: see the wiki.

## White noise

```ini
[random]
generator  = PANPHASIA
descriptor = [Panph6,L18,(56034,71400,250000),S1,KK1025,CH3774755196,Flamingo_Gpc1]
```

or

```ini
[random]
generator = MUSIC
seed[7]   = 12345
seed[8]   = 23456
```

- `generator = PANPHASIA`: noise set by `descriptor`, independent of the
  grid, so any code or resolution reproduces the same structure. See
  [PANPHASIA](panphasia.md).
- `generator = MUSIC`: MUSIC's own noise, with one `seed[L]` per level. Keep
  the seeds when changing `levelmax` to preserve structure.

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
