# ctest-labels: fast
"""In 1LPT the peculiar velocity is the displacement times a H(a) f(a), on every level of a zoom."""

import numpy as np

from musictest import check, finish, run_music


def expected_ratio(config):
    """100 a E(a) f(a) for flat LCDM, the factor MUSIC's generic output uses between v and dx."""
    om, ol = float(config["cosmology"]["Omega_m"]), float(config["cosmology"]["Omega_L"])
    a = 1.0 / (1.0 + float(config["setup"]["zstart"]))
    e = lambda x: np.sqrt(om / x**3 + ol)

    def growth(x):
        y = np.linspace(0, x, 100001)[1:]
        f = 1.0 / (y * e(y)) ** 3
        return e(x) * np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(y))

    dlnd = (np.log(growth(a * 1.0001)) - np.log(growth(a / 1.0001))) / (2 * np.log(1.0001))
    return 100 * a * e(a) * dlnd


for zstart in (100, 10):
    run = run_music(f"vel_z{zstart}", setup={"zstart": zstart, "levelmax": 8})
    if not check(f"z={zstart}: MUSIC runs", run.ok()):
        continue
    want = expected_ratio(run.config)
    for level in (7, 8):
        for c in "xyz":
            d, _ = run.level(level, f"d{c}")
            v, _ = run.level(level, f"v{c}")
            ratio = np.sum(v * d) / np.sum(d * d)
            scatter = np.std(v - ratio * d) / np.std(v)
            check(f"z={zstart} level {level} v{c} = {want:.3f} d{c}", abs(ratio / want - 1) < 1e-3 and scatter < 1e-6,
                  f"ratio {ratio:.4f}, scatter {scatter:.1e}")

finish()
