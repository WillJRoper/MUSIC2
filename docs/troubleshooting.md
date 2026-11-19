# Troubleshooting

MUSIC stops with a message on screen and in `<config>_log.txt`. Common ones:

## Build

| Message | Fix |
|---|---|
| `Could NOT find OpenMP` | Apple Clang has no OpenMP. Use Homebrew GCC: `CC=gcc-N CXX=g++-N cmake ...` |
| `Could NOT find FFTW3` / `GSL` | install them or pass `FFTW3_ROOT=<path>` / `GSL_ROOT_DIR=<path>` |
| no Fortran compiler | install `gfortran`, or `-DENABLE_PANPHASIA=OFF` |
| CLASS download fails | needs network at configure time; or `-DENABLE_CLASS=OFF` |
| `Could NOT find TIRPC` | harmless; only needed for Tipsy output |

## Configuration

| Message | Cause and fix |
|---|---|
| `... plug-in '<name>' not found` | misspelt `format`, `transfer`, `generator` or `region`, or that plugin was not built (e.g. PANPHASIA or CLASS switched off, or no HDF5 for `swift`) |
| `'cosmology/Omega_m' not found` | give `Omega_m`, `Omega_b`, `Omega_L`, `H0` explicitly, even with `ParameterSet` |
| `Found levelmin!=levelmax but neither ref_offset nor ref_center was specified` | a `box` zoom needs a position and a size, or set `region` |
| `ref_dims must be even` | use even `ref_dims` |
| `ParticleLoad must be 'sc' or 'masked'` | only those two loads exist for SWIFT |
| `... only supported for output format 'swift'` | the masked load is SWIFT-only |
| `ParticleMaskType must be 0, 1, 2 or 3` | see [Masked particle load](swift-output.md#masked-particle-load) |

## PANPHASIA

| Message | Cause and fix |
|---|---|
| `invalid descriptor` | the check digit (`CH...`) does not match; copy the descriptor exactly |
| `power-of-two base size` | `Panph6` needs `S` = 1, 2, 4, ... |
| `levelmin is too small for the descriptor's KK` | raise `levelmin` until 2^(levelmin-1) > √KK (e.g. 7 for KK1025) |
| `fix_mode_amplitude cannot be combined with a KK descriptor` | drop `fix_mode_amplitude`; `KK` already fixes those modes |

## CLASS

| Message | Cause and fix |
|---|---|
| `Shooting failed` / `You requested a very high z_pk` followed by `std::out_of_range` | `zstart` is above recombination; keep it below about 1088 |

## Results

- **Low-resolution particles in the zoomed halo**: the region is too small.
  Grow the radius used for the point file (2-3 × R200 is typical), or
  `padding`.
- **Zoom halo in the wrong place**: zoom positions are shifted relative to the
  parent; see the end of the [Quickstart](quickstart.md#3-zoom-on-a-halo).
- **`check_swift_masked.py` reports a large residual**: it assumes small
  displacements. Generate the file to check at high `zstart` (e.g. 1000).
