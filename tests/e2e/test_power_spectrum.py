# ctest-labels: fast
"""The density power spectrum of uniform runs matches MUSIC's input_powerspec.txt, which is
smaller than the measured P(k) by 8 pi^3, for each random number generator."""

import numpy as np

from musictest import EAGLE, FLAMINGO, check, finish, run_music, wavenumbers

GENERATORS = {
    "music": {},
    "panph1": {"generator": "PANPHASIA", "descriptor": EAGLE},
    "panph6": {"generator": "PANPHASIA", "descriptor": FLAMINGO},
}

for name, random in GENERATORS.items():
    run = run_music(f"pk_{name}", random=random)
    if not check(f"{name}: MUSIC runs", run.ok()):
        continue
    delta, _ = run.level(7, "rho")
    n, box = delta.shape[0], float(run.config["setup"]["boxlength"])
    pk = np.abs(np.fft.rfftn(delta)) ** 2 * box**3 / n**6
    knyq = np.pi * n / box
    k = wavenumbers(delta.shape) * knyq
    table = run.powerspec  # k, P_cdm, P_theta_cdm, P_total
    for lo, hi in [(0.05, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.7), (0.7, 0.95)]:
        s = (k > lo * knyq) & (k <= hi * knyq)
        ratio = np.mean(pk[s] / np.interp(k[s], table[:, 0], table[:, 3])) / (8 * np.pi**3)
        tol = 3 / np.sqrt(s.sum()) + 0.01  # cosmic variance plus 1%
        check(f"{name}: P(k)/P_input at {lo}-{hi} kNyq", abs(ratio - 1) < tol, f"{ratio:.4f}, tolerance {tol:.4f}")

finish()
