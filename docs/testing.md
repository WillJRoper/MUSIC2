# Testing

## End-to-end tests

The tests in `tests/e2e` run MUSIC on small configurations and check the
results:

| Test | Checks |
|---|---|
| `test_power_spectrum.py` | measured power spectrum against the input spectrum |
| `test_lpt_velocity.py` | 1LPT velocities against displacements |
| `test_zoom_consistency.py` | zooms against uniform runs |
| `test_panphasia_kk.py` | PANPHASIA `KK` mode fixing |
| `test_config_errors.py` | invalid configurations are rejected |
| `test_swift_output.py` | SWIFT output |
| `test_lagrangian_points.py` | `tools/lagrangian_points.py` output and an ellipsoid zoom built from it |
| `test_monofonic.py` | MUSIC against monofonIC for a `Panph6` descriptor |

They need `python3` with `numpy` and `h5py`. From the build directory:

```
ctest -L fast --output-on-failure
```

takes about a minute and also runs in CI. `ctest -L monofonic` runs the
monofonIC comparison; it is skipped unless `MONOFONIC_EXECUTABLE` points to a
monofonIC built with PANPHASIA.

Each test is a script that can also be run directly:

```
python3 tests/e2e/test_x.py path/to/MUSIC workdir
```

Shared helpers live in `tests/e2e/musictest.py`.

## Checking SWIFT files

`tests/scripts` holds standalone checkers for SWIFT initial conditions:

- `check_swift_cells.py ics.hdf5`: cell metadata, ordering, swiftsimio loading.
- `check_swift_masked.py ics.hdf5 <type> [unmasked.hdf5]`: masked particle
  load.

See [SWIFT output](swift-output.md) for details.
