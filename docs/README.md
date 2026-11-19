# MUSIC2 fork documentation

The [upstream wiki](https://github.com/cosmo-sims/music/wiki) remains the
general reference for MUSIC's parameters. These pages cover what this fork
adds or changes, plus parts of MUSIC that the wiki leaves thin.

- [Refinement regions](refinement-regions.md): the shapes the high-resolution
  region can take (`box`, `ellipsoid`, `convex_hull`, `multibox`) and their
  options.
- [SWIFT output](swift-output.md): the `swift` output plugin, its cell
  metadata, MPI writing, and the masked particle load.
- [PANPHASIA descriptors](panphasia.md): `Panph1` and high-order `Panph6`
  (e.g. FLAMINGO) descriptors.
- [Testing](testing.md): the end-to-end tests and the SWIFT checker scripts.

Example parameter files: [`ics_example.conf`](../ics_example.conf) (box
region, all output formats) and
[`data/ExampleConfigs/AURIGA6_ics.conf`](../data/ExampleConfigs/AURIGA6_ics.conf)
(ellipsoid region, PANPHASIA, CAMB transfer file).
