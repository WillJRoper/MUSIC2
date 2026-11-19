# ctest-labels: fast
"""tools/lagrangian_points.py: points of a sphere straddling the periodic boundary in uniform SWIFT ICs
come out in box units, and a zoom with an ellipsoid fitted to them runs."""

import os
import subprocess
import sys

import h5py
import numpy as np

from musictest import WORKDIR, check, finish, run_music

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools", "lagrangian_points.py")
SETUP = {"boxlength": 100, "levelmin": 6, "levelmin_TF": 6, "levelmax": 6}

parent = run_music("lp_parent", setup=SETUP, output={"format": "swift"})
if check("parent ICs written", parent.ok()):
    with h5py.File(parent.output, "r") as f:
        box = f["Header"].attrs["BoxSize"][0]
        pos = f["PartType1/Coordinates"][:]
    centre, radius = np.array([0.02, 0.5, 0.98]) * box, 0.08 * box
    points = os.path.join(WORKDIR, "points.txt")
    # the ICs stand in for the snapshot: the selected points are then the sphere itself
    proc = subprocess.run([sys.executable, TOOL, parent.output, parent.output, *map(str, centre), str(radius),
                           "-o", points], capture_output=True, text=True)
    if check("tool runs", proc.returncode == 0, proc.stderr[-500:]):
        pts = np.loadtxt(points)
        d = pts - centre / box
        d -= np.round(d)
        n_in = np.sum(np.sum(((pos - centre + box / 2) % box - box / 2) ** 2, axis=1) < radius**2)
        check("every particle in the sphere is written", len(pts) == n_in, f"{len(pts)} vs {n_in}")
        check("points are in box units inside the sphere",
              np.all((pts >= 0) & (pts < 1)) and np.all(np.sum(d**2, axis=1) < (radius / box) ** 2))

        zoom = run_music("lp_zoom", setup={**SETUP, "levelmax": 7, "ref_center": None, "ref_extent": None,
                                           "region": "ellipsoid", "region_point_file": points},
                         output={"format": "swift"})
        check("zoom with the fitted ellipsoid runs", zoom.ok())

finish()
