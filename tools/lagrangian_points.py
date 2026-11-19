#!/usr/bin/env python3
"""Write a MUSIC region point file for a halo in a uniform-resolution parent run.

Usage:
  lagrangian_points.py PARENT_ICS SNAPSHOT X Y Z RADIUS [-o points.txt] [--parttype 1]

Selects the particles of PARTTYPE within RADIUS of (X, Y, Z) in SNAPSHOT (all in
the snapshot's own length units, periodic), finds the same particle IDs in
PARENT_ICS, and writes their initial-condition positions in box units (the unit
cube), one "x y z" line per particle. Use the result as `region_point_file` with
`region = ellipsoid`, `convex_hull` or `multibox`.

Both files must be SWIFT-style HDF5 (PartTypeN/Coordinates, PartTypeN/ParticleIDs,
Header/BoxSize), e.g. MUSIC's swift output and a SWIFT snapshot. The parent run
must be uniform (levelmin = levelmax), so its ICs are not shifted.
"""

import argparse

import h5py
import numpy as np


def box_size(f):
    return float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])


def lagrangian_points(ics, snapshot, centre, radius, parttype=1):
    group = f"PartType{parttype}"
    with h5py.File(snapshot, "r") as f:
        box = box_size(f)
        d = f[group]["Coordinates"][:] - centre
        d -= box * np.round(d / box)  # periodic separation
        ids = f[group]["ParticleIDs"][:][np.einsum("ij,ij->i", d, d) < radius**2]
    if len(ids) == 0:
        raise SystemExit("no particles within the radius")

    # ponytail: reads every ICs ID into memory; chunk over the file for very large parents
    with h5py.File(ics, "r") as f:
        box = box_size(f)
        ics_ids = f[group]["ParticleIDs"][:]
        order = np.argsort(ics_ids)
        idx = np.searchsorted(ics_ids, ids, sorter=order)
        idx = order[np.minimum(idx, len(order) - 1)]
        if not np.array_equal(ics_ids[idx], ids):
            raise SystemExit("some snapshot IDs are missing from the ICs")
        pos = f[group]["Coordinates"][:][np.sort(idx)]
    return np.mod(pos / box, 1.0), len(ids)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("ics")
    p.add_argument("snapshot")
    p.add_argument("centre", type=float, nargs=3, metavar="X/Y/Z")
    p.add_argument("radius", type=float)
    p.add_argument("-o", "--output", default="points.txt")
    p.add_argument("--parttype", type=int, default=1)
    a = p.parse_args()
    pts, n = lagrangian_points(a.ics, a.snapshot, np.array(a.centre), a.radius, a.parttype)
    np.savetxt(a.output, pts, fmt="%.8f")
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    print(f"wrote {n} points to {a.output}; extent {hi - lo} of the box "
          "(a value near 1 means the region wraps the periodic boundary, which is fine)")


if __name__ == "__main__":
    main()
