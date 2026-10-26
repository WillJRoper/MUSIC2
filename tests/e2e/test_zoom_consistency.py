# ctest-labels: fast
"""The finest level of a zoom reproduces a uniform run at that resolution inside the refined region.

PANPHASIA defines the field independently of the grid, so zoom and uniform runs must agree on all
scales. MUSIC's own generator only shares the large-scale modes between the two."""

import numpy as np

from musictest import EAGLE, FLAMINGO, check, cross_spectrum, cut, finish, run_music

BINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9]
# random section, coarse level of the zoom (FLAMINGO's KK needs 7), and the fraction of Nyquist up
# to which the modes must agree
GENERATORS = {
    "panph1": ({"generator": "PANPHASIA", "descriptor": EAGLE}, 6, 0.9),
    "panph6": ({"generator": "PANPHASIA", "descriptor": FLAMINGO}, 7, 0.9),
    "music": ({"seed[6]": 54321}, 6, 0.3),
}

for name, (random, lmin, kmax) in GENERATORS.items():
    fine = lmin + 1
    uniform = run_music(f"uni_{name}", setup={"levelmin": fine, "levelmin_TF": fine, "levelmax": fine}, random=random)
    extent = 0.2 * 2 ** (7 - lmin)  # the same number of fine cells for every coarse level
    zoom = run_music(f"zoom_{name}", setup={"levelmin": lmin, "levelmin_TF": lmin, "levelmax": fine,
                                            "ref_extent": f"{extent}, {extent}, {extent}"}, random=random)
    if not check(f"{name}: MUSIC runs", uniform.ok() and zoom.ok()):
        continue
    full, _ = uniform.level(fine, "dx")
    patch, offset = zoom.level(fine, "dx")
    m = 4  # skip the cells next to the patch boundary
    patch = patch[m:-m, m:-m, m:-m]
    for lo, hi, r, amp in cross_spectrum(patch, cut(full, offset + m, patch.shape), BINS):
        rmin = 0.999 if lo == 0 else 0.9999  # the lowest bin has few modes, all near the patch size
        if hi <= kmax:
            check(f"{name}: zoom matches uniform at {lo}-{hi} kNyq", r > rmin and abs(amp - 1) < 0.01,
                  f"r = {r:.6f}, amplitude ratio {amp:.4f}")

finish()
