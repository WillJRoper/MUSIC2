# ctest-labels: fast
"""A Panph6 descriptor's KK limit fixes every mode with |k|^2 <= KK (box units) to the mean power,
and the descriptor may be written with spaces between its fields."""

import numpy as np

from musictest import FLAMINGO, check, finish, run_music

run = run_music("kk", random={"generator": "PANPHASIA", "descriptor": FLAMINGO.replace(",", ", ")})
check("descriptor with spaces is accepted", run.ok())
finish()

delta, _ = run.level(7, "rho")
n, box = delta.shape[0], float(run.config["setup"]["boxlength"])
pk = np.abs(np.fft.rfftn(delta)) ** 2 * box**3 / n**6
k1, kz = np.fft.fftfreq(n) * n, np.fft.rfftfreq(n) * n
k2 = k1[:, None, None] ** 2 + k1[None, :, None] ** 2 + kz[None, None, :] ** 2
kphys = np.sqrt(k2) * 2 * np.pi / box
ratio = pk / np.interp(kphys, run.powerspec[:, 0], run.powerspec[:, 3]) / (8 * np.pi**3)

# fixed modes all have the same power at a given |k|, and the mean power up to the interpolation
# of the tabulated spectrum
fixed = (k2 > 0) & (k2 <= 1025)
shells = [ratio[k2 == s] for s in np.unique(k2[fixed])]
spread = max(np.ptp(r) / r.mean() for r in shells)
check("KK modes have the same power at each |k|", spread < 1e-6, f"max spread {spread:.1e}")
check("KK modes have the mean power", np.allclose(ratio[fixed], 1, atol=0.02), f"max deviation {np.abs(ratio[fixed] - 1).max():.1e}")
free = (k2 > 1025) & (k2 <= 4 * 1025)
check("modes beyond KK are not fixed", np.std(ratio[free]) > 0.5, f"scatter {np.std(ratio[free]):.2f}")
finish()
