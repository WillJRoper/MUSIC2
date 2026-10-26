# ctest-labels: monofonic
"""With the same Panph6 descriptor, MUSIC reproduces monofonIC's displacement field, both in a
uniform run and on the finest level of a zoom.

monofonIC places particles at cell corners and MUSIC at cell centres, so monofonIC's field is
shifted by half a cell before comparing. The two codes normalise the amplitude slightly
differently, so only the shape of the amplitude ratio is checked.

Needs MONOFONIC_EXECUTABLE pointing at a monofonIC built with PANPHASIA (which needs MPI)."""

import os
import subprocess

import h5py
import numpy as np

from musictest import FLAMINGO, WORKDIR, check, cross_spectrum, cut, finish, run_music, skip

MONOFONIC = os.environ.get("MONOFONIC_EXECUTABLE")
if not MONOFONIC:
    skip("MONOFONIC_EXECUTABLE is not set")

PANPH6 = {"generator": "PANPHASIA", "descriptor": FLAMINGO}
uniform = run_music("music_uniform", setup={"levelmin": 8, "levelmin_TF": 8, "levelmax": 8}, random=PANPH6)
zoom = run_music("music_zoom", setup={"levelmax": 8}, random=PANPH6)
c = uniform.config
with open(os.path.join(WORKDIR, "monofonic.conf"), "w") as f:
    f.write(f"""[setup]
GridRes = 256
BoxLength = {c["setup"]["boxlength"]}
zstart = {c["setup"]["zstart"]}
LPTorder = 1
DoBaryons = no
DoFixing = no
ParticleLoad = sc

[cosmology]
ParameterSet = none
Omega_m = {c["cosmology"]["Omega_m"]}
Omega_b = {c["cosmology"]["Omega_b"]}
Omega_L = {c["cosmology"]["Omega_L"]}
H0 = {c["cosmology"]["H0"]}
n_s = {c["cosmology"]["n_s"]}
sigma_8 = {c["cosmology"]["sigma_8"]}
m_nu1 = 0
m_nu2 = 0
m_nu3 = 0
transfer = eisenstein

[random]
generator = PANPHASIA
descriptor = {FLAMINGO}

[output]
format = generic
filename = monofonic.hdf5
""")
mono = subprocess.run([MONOFONIC, "monofonic.conf"], cwd=WORKDIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
if not check("monofonIC runs", mono.returncode == 0, mono.stdout[-2000:] if mono.returncode else ""):
    finish()
if not check("MUSIC runs", uniform.ok() and zoom.ok()):
    finish()

with h5py.File(os.path.join(WORKDIR, "monofonic.hdf5"), "r") as f:
    reference = f["DM_dx"][()]
n = reference.shape[0]
k = [np.fft.fftfreq(n) * n] * 2 + [np.fft.rfftfreq(n) * n]
shift = np.exp(1j * np.pi * (k[0][:, None, None] + k[1][None, :, None] + k[2][None, None, :]) / n)
shift[n // 2], shift[:, n // 2], shift[:, :, n // 2] = 0, 0, 0
reference = np.fft.irfftn(np.fft.rfftn(reference) * shift, s=reference.shape, axes=(0, 1, 2))

bins = [0.0, 0.1, 0.2, 0.4, 0.6, 0.9]
full, _ = uniform.level(8, "dx")
for lo, hi, r, amp in cross_spectrum(full, reference, bins, window=False):
    check(f"uniform matches monofonIC at {lo}-{hi} kNyq", r > 0.999999, f"r = {r:.7f}, amplitude ratio {amp:.4f}")
patch, offset = zoom.level(8, "dx")
patch = patch[4:-4, 4:-4, 4:-4]
spectra = cross_spectrum(patch, cut(reference, offset + 4, patch.shape), bins)
for lo, hi, r, amp in spectra:
    check(f"zoom matches monofonIC at {lo}-{hi} kNyq", r > (0.999 if lo == 0 else 0.9999), f"r = {r:.6f}, amplitude ratio {amp:.4f}")
amps = [amp for _, _, _, amp in spectra]
check("amplitude ratio is the same at all k", np.ptp(amps) / np.mean(amps) < 0.01, f"{min(amps):.4f}-{max(amps):.4f}")

finish()
