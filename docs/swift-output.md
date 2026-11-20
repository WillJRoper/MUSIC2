# SWIFT output

The `swift` output plugin writes one SWIFT-compatible HDF5 initial-conditions
file. Every particle type is sorted into a grid of SWIFT top-level cells, and
the `/Cells` metadata lets SWIFT and
[swiftsimio](https://swiftsimio.readthedocs.io/) load spatial subsets.

## Options

```ini
[output]
format           = swift
filename         = ics_swift.hdf5
top_level_cells  = 32
swift_longids    = yes
swift_coarsetype = 2
swift_unitlength = 3.08567758149e24
swift_unitmass   = 1.98841e43
swift_unitvel    = 1e5
```

- `top_level_cells`: cells along each axis, 1-1290. Default 32. More cells
  give finer spatial selection, but the metadata grows as the cube.
- `swift_longids`: 64-bit particle IDs. Default `yes`.
- `swift_coarsetype`: particle type for the coarse background particles.
  Default 2.
- `swift_unitlength`, `swift_unitmass`, `swift_unitvel`: code units in cgs.
  Defaults are Mpc, 10^10 M☉ and km/s.

Fields are always double precision and always in one file, so
`swift_doubleprec` and `swift_num_files` are ignored.

## MPI

Build with `-DENABLE_MPI=ON` (the default). The grids are computed on every
rank; the finished particles are then split among ranks for writing. With
parallel HDF5, ranks write their parts at the same time; with serial HDF5
they take turns. Output is always one file.

Multi-rank runs need `format = swift`; other output plugins do not support
MPI.

## File contents

- `/Header`: corrected scale factor, three-component box size, counts for
  seven particle types, `Music_levelcounts`, and a zero `MassTable`. Every
  particle has its own mass, including background DM.
- `/Cells/{Counts,OffsetsInFile,Files,MinPositions,MaxPositions}` for every
  particle type. The min/max positions are the cell bounds, not particle
  extremes, so masked reads need no padding. Cells cover the whole box for
  every type.
- SWIFT unit exponents on every particle dataset.
- `/ICs_parameters/Config`: every MUSIC option used.
- `/ICs_parameters/Cosmology`, `/Build`, `/Derived`: provenance and derived
  quantities such as particle masses.
- `/ICs_parameters/InputPowerSpectrum`: MUSIC's input power spectrum.
- `/ICs_parameters/CLASS` and the target redshift for CLASS runs, and the
  PANPHASIA descriptor when used.

Check a file with

```
python3 tests/scripts/check_swift_cells.py ics_swift.hdf5
```

which checks cell counts, offsets, ordering and metadata, and, if
`swiftsimio` is installed, that it loads. Needs `h5py` and `numpy`.

## Masked particle load

```ini
[setup]
baryons          = yes
ParticleLoad     = masked
ParticleMaskType = 3
```

By default the high-resolution region has one DM and one gas particle per
finest cell, the gas offset by half a cell. The masked load, taken from
monofonIC, instead fills each 2x2x2 block of finest cells with one gas
particle at the block centre and DM in some of the other cells.

- `baryons`: must be `yes`; otherwise `masked` is ignored with a warning.
- `ParticleLoad`: `sc` (default) or `masked`. `masked` works only with
  `format = swift` and without a glass.
- `ParticleMaskType`: which cells of each block get DM. Default 3.

| `ParticleMaskType` | DM cells per block | DM:gas |
|---|---|---|
| 0 | corner | 1:1 |
| 1 | corner + edges | 4:1 |
| 2 | corner + faces | 4:1 |
| 3 | corner + edges + faces | 7:1 |

Masses are scaled so each species keeps its total mass. The gas smoothing
length is the 2-cell gas spacing, and the suggested softening uses the DM
spacing. Coarse levels still hold combined DM + gas particles.

Check a masked file with

```
python3 tests/scripts/check_swift_masked.py ics_swift.hdf5 <type> [unmasked.hdf5]
```

which checks particle counts, gas positions, smoothing lengths and, given an
unmasked run with the same seeds, total masses. It assumes small
displacements, so generate the file at high `zstart` (e.g. 1000; CLASS fails
above about 1088).
