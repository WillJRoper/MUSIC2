# SWIFT output

The `swift` output plugin writes one SWIFT-compatible HDF5 initial-conditions
file. Particle floating-point fields are always double precision. Every
particle type is stably sorted into a regular grid of SWIFT top-level cells;
the corresponding `/Cells` metadata allows SWIFT and
[swiftsimio](https://swiftsimio.readthedocs.io/) to load spatial subsets.

```ini
[output]
format = swift
filename = ics_swift.hdf5
swift_longids = yes
top_level_cells = 32
```

## Options

| Option | Default | Meaning |
|---|---|---|
| `top_level_cells` | 32 | cells along each axis, 1-1290. Larger values improve spatial selectivity but metadata grows as its cube. |
| `swift_longids` | `yes` | 64-bit particle IDs |
| `swift_coarsetype` | 2 | particle type for coarse (background) particles |
| `swift_unitlength`, `swift_unitmass`, `swift_unitvel` | | output code units |
| `swift_num_files` | | ignored; a single file is required for one global cell tree |
| `swift_doubleprec` | | ignored; fields are always double precision |

## MPI

MUSIC2's grid computation is replicated across MPI ranks; SWIFT output splits
the completed particle fields among ranks during finalization. Build with
`-DENABLE_MPI=ON`. When the HDF5 library has parallel support, ranks write
collectively to disjoint hyperslabs; with serial HDF5, ranks take turns
writing their ranges. Output is always one shared file.

Multi-rank runs currently require `output/format = swift`; the other output
plugins do not support this MPI writer path.

## File contents

- `/Cells/{Counts,OffsetsInFile,Files,MinPositions,MaxPositions}` for every
  particle type. `MinPositions`/`MaxPositions` are geometric cell bounds, not
  particle extrema, so masked reads do not require padding. Cell bounds cover
  the whole periodic box for every particle type, including zoom
  high-resolution and background particles.
- SWIFT unit exponents on every particle dataset.
- `/ICs_parameters/Config`, containing every stored MUSIC option.
- `/ICs_parameters/Cosmology`, `/Build`, and `/Derived` provenance.
- `/ICs_parameters/InputPowerSpectrum`, containing MUSIC's sampled input table.
- `/ICs_parameters/CLASS` and target redshift for CLASS transfer runs, plus the
  PANPHASIA descriptor when that RNG is selected.
- In `/Header`: corrected scale factor, three-component box size, seven
  particle-type counts, `Music_levelcounts`, and zero `MassTable` values. All
  masses are written per particle, including single-level background dark
  matter.

Check a file with

```
python3 tests/scripts/check_swift_cells.py ics_swift.hdf5
```

which verifies cell counts, offsets, ordering, metadata, and (when installed)
swiftsimio loading. It needs `h5py` and `numpy`; `swiftsimio` is optional.

## Masked particle load

By default the high-resolution region holds one dark matter and one gas
particle per finest-level cell, with the gas staggered by half a cell. With
baryons enabled, the SWIFT plugin can instead use monofonIC's masked load:

```ini
[setup]
baryons          = yes
ParticleLoad     = masked
ParticleMaskType = 3
```

Each 2x2x2 block of finest-level cells then holds one gas particle in the
block-centre cell and dark matter in the corner cell plus the edge-centre
cells (`ParticleMaskType` bit 0) and/or face-centre cells (bit 1):

| `ParticleMaskType` | DM cells per block | DM:gas |
|---|---|---|
| 0 | corner | 1:1 |
| 1 | corner + edges | 4:1 |
| 2 | corner + faces | 4:1 |
| 3 | corner + edges + faces | 7:1 |

Gas is not staggered, particle masses are scaled so each species keeps its
total mass, and the gas smoothing length is the 2-cell gas spacing. Coarse
levels still hold total-matter particles. The suggested high-resolution
softening is scaled with the masked DM spacing.

Restrictions: only for `swift` output and non-glass loads. `ParticleLoad`
accepts only `sc` (the default) and `masked`; `masked` is ignored with a
warning when baryons are disabled.

Check a masked file, generated at a high starting redshift, with

```
python3 tests/scripts/check_swift_masked.py ics_swift.hdf5 <type> [unmasked.hdf5]
```

which checks particle counts, gas lattice positions, smoothing lengths, and,
given an unmasked run with the same seeds, total masses.
