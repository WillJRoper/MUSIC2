#!/usr/bin/env python3
"""Validate a SWIFT IC file written with 'ParticleLoad = masked'.

Usage: check_swift_masked.py ics.hdf5 MASK_TYPE [UNMASKED_REFERENCE.hdf5]

Checks the DM:gas ratio on the finest level, that every gas particle sits
near the centre cell of its 2x2x2 block (generate at high zstart so
displacements are small), and, given an unmasked run with the same seeds,
that the total high-res mass per species and the coarse particles agree.
"""

import sys

import h5py
import numpy as np

DM_PER_BLOCK = {0: 1, 1: 4, 2: 4, 3: 7}


def check(path, mask_type, reference=None):
    with h5py.File(path, "r") as f:
        params = f["ICs_parameters/Config/setup"].attrs
        levelmax = int(params["levelmax"])
        box = f["Header"].attrs["BoxSize"][0]
        dx = box / 2**levelmax

        gas = f["PartType0/Coordinates"][:]
        n_gas, n_dm = len(gas), len(f["PartType1/Coordinates"])
        assert n_dm == DM_PER_BLOCK[mask_type] * n_gas, (n_dm, n_gas)

        # gas lives on cell centres of the odd (i,j,k) cell in each block,
        # i.e. at 1.5 dx modulo 2 dx
        residual = np.mod(gas - 1.5 * dx + dx, 2 * dx) - dx
        assert np.abs(residual).max() < 0.25 * dx, np.abs(residual).max()

        # the high-res gas smoothing length is the gas interparticle spacing
        assert np.allclose(f["PartType0/SmoothingLength"][:], 2 * dx)

        masses = {k: f[f"{k}/Masses"][:].sum() for k in f if k.startswith("PartType")}

    if reference:
        with h5py.File(reference, "r") as f:
            ref = {k: f[f"{k}/Masses"][:].sum() for k in f if k.startswith("PartType")}
        assert masses.keys() == ref.keys(), (masses.keys(), ref.keys())
        for k in masses:
            assert np.isclose(masses[k], ref[k], rtol=1e-10), (k, masses[k], ref[k])

    print(f"OK {path}: {n_dm} DM, {n_gas} gas")


if __name__ == "__main__":
    check(sys.argv[1], int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
