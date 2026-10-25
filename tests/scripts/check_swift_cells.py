#!/usr/bin/env python3
"""Validate SWIFT cell ordering and optional swiftsimio compatibility."""

import sys

import h5py
import numpy as np


def check(path):
    with h5py.File(path, "r") as f:
        dimensions = f["Cells/Meta-data"].attrs["dimension"]
        box = f["Header"].attrs["BoxSize"]
        assert len(dimensions) == 3 and dimensions[0] == dimensions[1] == dimensions[2]
        assert np.allclose(box, box[0])
        cdim = int(dimensions[0])

        for key in f:
            if not key.startswith("PartType") or "Coordinates" not in f[key]:
                continue
            coordinates = f[f"{key}/Coordinates"][:]
            cell_counts = f[f"Cells/Counts/{key}"][:]
            offsets = f[f"Cells/OffsetsInFile/{key}"][:]
            assert cell_counts.sum() == len(coordinates), key
            for field in ("Velocities", "Masses", "ParticleIDs"):
                assert len(f[f"{key}/{field}"]) == len(coordinates), f"{key}/{field}"
            ids_dtype = f[f"{key}/ParticleIDs"].dtype
            assert ids_dtype.kind == "u" and ids_dtype.itemsize in (4, 8), key
            assert np.array_equal(offsets, np.r_[0, np.cumsum(cell_counts[:-1])]), key

            cells = np.floor(coordinates * cdim / box[0]).astype(np.int64)
            cells = np.minimum(cells, cdim - 1)
            ids = cells[:, 2] + cdim * (cells[:, 1] + cdim * cells[:, 0])
            expected = np.repeat(np.arange(len(cell_counts)), cell_counts)
            assert np.array_equal(ids, expected), key
            lower = f[f"Cells/MinPositions/{key}"][:].reshape(-1, 3)
            upper = f[f"Cells/MaxPositions/{key}"][:].reshape(-1, 3)
            assert np.all(coordinates >= lower[ids]) and np.all(coordinates < upper[ids]), key
            assert f[f"Cells/Files/{key}"].shape == cell_counts.shape
            assert np.all(f[f"Cells/Files/{key}"][:] == 0)

    try:
        import swiftsimio
    except ImportError:
        print("Cell tree validated (swiftsimio unavailable)")
    else:
        data = swiftsimio.load(path)
        mask = swiftsimio.mask(path)
        boxsize = float(data.metadata.boxsize[0].value)
        mask.constrain_spatial([[0.0, boxsize / 2]] * 3)
        subset = swiftsimio.load(path, mask)
        assert 0 < subset.metadata.n_dark_matter <= data.metadata.n_dark_matter
        print(f"Cell tree and swiftsimio validated: {data.metadata.n_dark_matter} DM particles")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {sys.argv[0]} FILE.hdf5")
    check(sys.argv[1])
