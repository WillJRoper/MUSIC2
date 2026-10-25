[![CMake on multiple platforms](https://github.com/cosmo-sims/MUSIC/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/cosmo-sims/MUSIC/actions/workflows/cmake-multi-platform.yml)

MUSIC v2 - multi-scale cosmological initial conditions
======================================================

MUSIC is a computer program to generate nested grid initial conditions for high-resolution "zoom" cosmological simulations. A detailed description of some of the algorithms can be found in [Hahn & Abel (2011)](http://arxiv.org/abs/1103.6031), but version 2.0 is much improved and has a much lower memory footprint. 

*Be aware that v1 and v2 do not produce entirely consistent zooms, due to the changes in the way the noise hierarchy is composed. So do not upgrade if you need to exactly reproduce a zoom run generated with v1, there is currently no backwards compatibility mode.*

To learn how to use MUSIC, please [read the Wiki](https://github.com/cosmo-sims/music/wiki). 

Current MUSIC key features include:

- Supports output for RAMSES, ENZO, Arepo, Swift, Gadget-2/3, ART, Pkdgrav/Gasoline, NyX, and GAMER-2 via plugins. New codes can be added.
- Support for first (1LPT) and second order (2LPT) Lagrangian perturbation theory.
- Full integration of the [PANPHASIA](https://arxiv.org/abs/1306.5968) hierarchical Gaussian random fields
- Pluggable transfer functions, currently CAMB, Eisenstein&Hu, BBKS, Warm Dark Matter variants. Distinct baryon+CDM fields.
- Minimum bounding ellipsoid and convex hull shaped high-res regions supported with most codes, supports refinement mask generation for RAMSES.
- Parallelized with OpenMP
- Requires FFTW3, GSL (and HDF5 for output for some codes)

## Building MUSIC
While we still supply the old Makefile, using CMake is now the preferred way of building. CMake use out-of-source build, i.e. you create a build directory, and then configure the code using CMake. Inside the `music` directory, do
```
  mkdir build
  cd build
  ccmake ..
  make -j
```
to configure the code (you will se a menu), and then start a parallel compilation. If CMake has trouble finding your FFTW or HDF5 installation, you can add hints as follows
```
  FFTW3_ROOT=<path> HDF5_ROOT=<path> ccmake ..
```
If you want to build on macOS, then it is strongly recommended to use GNU (or Intel) compilers instead of Apple's Clang. Install them e.g. via homebrew and then configure cmake to use them instead of the macOS default compiler via
```
  CC=gcc-13 CXX=g++-13 ccmake ..
```
This is necessary since Apple's compilers haven't supported OpenMP for years.


## Running
There is an example parameter file 'example.conf' in the main directory. Possible options are explained in it, it can be run as a simple argument, e.g. from within the build directory:
```
  ./MUSIC ../ics_example.conf
```

## SWIFT output

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

`top_level_cells` is number of cells along each axis (default 32, valid range
1-1290). Larger values improve spatial selectivity but increase metadata as its
cube. Output is always one shared file. When MUSIC is built with MPI and
parallel HDF5, ranks write collectively to disjoint hyperslabs; with serial
HDF5, ranks take turns writing their ranges. MUSIC's current SWIFT tree uses
whole-box geometry for every particle type, including zoom high-resolution and
background particles.

MUSIC2's grid computation is replicated across MPI ranks; SWIFT output splits
the completed particle fields among ranks during finalization. The legacy
`swift_num_files` option is ignored because a single file is required for one
global cell tree. Particle fields are always double precision; `swift_doubleprec`
does not change output precision. `swift_longids` defaults to `yes`, and
`swift_coarsetype` defaults to particle type 2. Optional `swift_unitlength`,
`swift_unitmass`, and `swift_unitvel` set output code units.

Multi-rank runs currently require `output/format = swift`; the other legacy
output plugins do not support this MPI writer path.

The file also records:

- `/Cells/{Counts,OffsetsInFile,Files,MinPositions,MaxPositions}` for every particle type;
- SWIFT unit exponents on every particle dataset;
- `/ICs_parameters/Config`, containing every stored MUSIC option;
- `/ICs_parameters/Cosmology`, `/Build`, and `/Derived` provenance;
- `/ICs_parameters/InputPowerSpectrum`, containing MUSIC's sampled input table;
- `/ICs_parameters/CLASS` and target redshift for CLASS transfer runs, plus the
  PANPHASIA descriptor when that RNG is selected;
- corrected scale factor, three-component box size, seven particle-type counts, and zero `MassTable` values in `/Header`.

All masses are written per particle, including single-level background dark
matter. `/Cells/MinPositions` and `MaxPositions` contain geometric cell bounds,
not particle extrema, so masked reads do not require padding.

Enable MPI for the parallel writer with `-DENABLE_MPI=ON`. Collective writes
are enabled automatically only when detected HDF5 library has parallel
support. Use `python3 tests/scripts/check_swift_cells.py ics_swift.hdf5` to
check cell counts, offsets, ordering, metadata, and (when installed) swiftsimio
loading. The checker requires `h5py` and `numpy`; `swiftsimio` is optional.


## Disclaimer
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. By downloading and using MUSIC, you agree to the LICENSE, distributed with the source code in a text file of the same name.
