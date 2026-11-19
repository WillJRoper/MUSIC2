# Quickstart

From a fresh clone to SWIFT zoom initial conditions for a chosen halo.

## 1. Build

Dependencies:

| Needed for | Packages |
|---|---|
| always | CMake ≥ 3.11, a C++17 compiler with OpenMP, FFTW3 (with its OpenMP/threads libraries), GSL |
| SWIFT and other HDF5 outputs | HDF5 |
| PANPHASIA (on by default) | a Fortran compiler |
| MPI SWIFT writing (on by default if found) | MPI; parallel HDF5 for collective writes |
| CLASS transfer function (on by default) | network access at configure time: CMake downloads CLASS |
| tests and checkers | Python 3 with `numpy`, `h5py`; optionally `swiftsimio` |

On macOS with Homebrew: `brew install cmake gcc fftw gsl hdf5 open-mpi`.

```
git clone https://github.com/WillJRoper/MUSIC2.git
cd MUSIC2
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
```

On macOS use Homebrew's GCC instead of Apple Clang, which ships without OpenMP:

```
CC=gcc-16 CXX=g++-16 cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

(use whichever `gcc-N` Homebrew installed). If CMake cannot find FFTW or HDF5,
point it at them with `FFTW3_ROOT=<path> HDF5_ROOT=<path> cmake ...`.

CMake options, all `ON` by default:

| Option | Effect when `OFF` |
|---|---|
| `ENABLE_MPI` | no MPI; SWIFT output is written by one process |
| `ENABLE_PANPHASIA` | no `generator = PANPHASIA`; no Fortran compiler needed |
| `ENABLE_CLASS` | no `transfer = CLASS`; nothing is downloaded |

The configure log reports whether HDF5 is parallel (`Serial HDF5: MPI ranks will
take turns writing SWIFT output` otherwise). Check the build with
`cd build && ctest -L fast` (about a minute, see [Testing](testing.md)).

## 2. First run

```
cd build
./MUSIC ../ics_example.conf
```

MUSIC writes its output in the current directory, plus a log
(`<config>_log.txt`) and some side files (see
[Parameters](parameters.md#files-music-writes)).

[`data/ExampleConfigs/swift_zoom.conf`](../data/ExampleConfigs/swift_zoom.conf)
is a complete fork-style zoom: SWIFT output, gas with the masked particle load,
`Panph6` PANPHASIA noise, a CLASS transfer function, and a box-shaped zoom
region. It takes about two minutes on a laptop:

```
./MUSIC ../data/ExampleConfigs/swift_zoom.conf
python3 ../tests/scripts/check_swift_cells.py ics_swift_zoom.hdf5
```

At the end MUSIC prints the SWIFT compile and runtime settings the file needs
(`LONGIDS`, `INPUT_IN_DOUBLEPRECISION`, suggested softening and PM grid, ...).
The suggested softening and PM grid are also stored in `/Header`
(`suggested_highressoft`, `suggested_pmgrid`).

## 3. Zoom on a halo

A zoom resimulates one object of a uniform-resolution parent run at higher
resolution, using the same white noise so the large-scale structure is
unchanged.

1. **Parent ICs.** Make uniform ICs (`levelmin = levelmax`, e.g. 7 or 8) with
   the `[random]` and `[cosmology]` sections you will use for the zoom. With
   PANPHASIA the noise is fixed by the descriptor; with MUSIC's own generator,
   keep the same `seed[n]` values.
2. **Parent run.** Run SWIFT (or any code) to the target redshift and find
   the halo, e.g. with a halo finder. Note its centre and a radius, usually
   2-3 × R200, in the snapshot's length units.
3. **Lagrangian region.** Map the halo's particles back to their positions in
   the parent ICs:

   ```
   python3 tools/lagrangian_points.py parent_ics.hdf5 snapshot_0010.hdf5 X Y Z RADIUS -o points.txt
   ```

   This selects the particles inside the sphere, matches their IDs in the
   ICs, and writes their initial positions in box units. Both files must be
   SWIFT-style HDF5 (MUSIC's `swift` output and SWIFT snapshots are). For
   other codes, write `x y z` in box units per line yourself; see
   [Point files](refinement-regions.md#point-files).
4. **Zoom ICs.** Copy the parent config, raise `levelmax`, and replace the
   box region with the point file:

   ```ini
   [setup]
   levelmin          = 7    # same as the parent
   levelmax          = 10
   region            = ellipsoid      # or convex_hull
   region_point_file = points.txt
   ```

   `ellipsoid` is the usual choice. `convex_hull` fits elongated regions more
   tightly. Both handle regions that wrap the periodic boundary. See
   [Refinement regions](refinement-regions.md).
5. **Run the zoom** in SWIFT. Check that no low-resolution particles end up
   in the halo; if they do, grow the radius in step 3 or the `padding`.

MUSIC recentres the zoom region in the box and logs
`Domain will be shifted by (i, j, k)`. Zoom positions are offset from the
parent's by `(i, j, k) × boxlength / 2^levelmin`, so subtract that to compare
the two runs. The shift is only recorded in the log.
