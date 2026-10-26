"""Helpers for MUSIC end-to-end tests.

Each test script is run by ctest as `python3 test_x.py MUSIC_EXECUTABLE WORKDIR`. It runs MUSIC on
small configurations inside WORKDIR, checks the output and exits non-zero if a check fails, or with
SKIP (77) if something it needs is not available.
"""

import copy
import os
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np

SKIP = 77

MUSIC = os.path.abspath(sys.argv[1])
WORKDIR = os.path.abspath(sys.argv[2])
shutil.rmtree(WORKDIR, ignore_errors=True)
os.makedirs(WORKDIR)

# descriptors with valid check digits: EAGLE (original PANPHASIA) and FLAMINGO (high order, with KK)
EAGLE = "[Panph1,L16,(31250,23438,39063),S12,CH1050187043,EAGLE_L0100_VOL1]"
FLAMINGO = "[Panph6,L18,(56034,71400,250000),S1,KK1025,CH3774755196,Flamingo_Gpc1]"

BASE = {
    "setup": {
        "boxlength": 1000,
        "zstart": 100,
        "levelmin": 7,
        "levelmin_TF": 7,
        "levelmax": 7,
        "padding": 8,
        "overlap": 4,
        "ref_center": "0.3, 0.6, 0.45",
        "ref_extent": "0.2, 0.2, 0.2",
        "align_top": "no",
        "baryons": "no",
        "use_2LPT": "no",
    },
    "cosmology": {
        "Omega_m": 0.305,
        "Omega_L": 0.695,
        "Omega_b": 0.045,
        "H0": 67.77,
        "sigma_8": 0.811,
        "n_s": 0.961,
        "transfer": "eisenstein",
    },
    "random": {"generator": "MUSIC", "seed[7]": 12345, "seed[8]": 23456, "seed[9]": 34567},
    "output": {"format": "generic"},
    "poisson": {"fft_fine": "yes", "accuracy": 1e-5, "laplace_order": 6, "grad_order": 6},
}

failures = []


def check(what, ok, detail=""):
    """Record and print one check."""
    print(f"{'PASS' if ok else 'FAIL'}: {what}" + (f" ({detail})" if detail else ""))
    if not ok:
        failures.append(what)
    return ok


def finish():
    if failures:
        print(f"{len(failures)} check(s) failed")
        sys.exit(1)


def skip(reason):
    finish()  # report failures found before the skip
    print(f"SKIP: {reason}")
    sys.exit(SKIP)


class Run:
    """One MUSIC run: its log, output file and grid hierarchy."""

    def __init__(self, name, config):
        self.name = name
        self.config = config
        conf = os.path.join(WORKDIR, f"{name}.conf")
        with open(conf, "w") as f:
            for section, values in config.items():
                f.write(f"[{section}]\n")
                f.writelines(f"{k} = {v}\n" for k, v in values.items())
                f.write("\n")
        proc = subprocess.run([MUSIC, conf], cwd=WORKDIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self.returncode = proc.returncode
        self.log = proc.stdout
        self.output = os.path.join(WORKDIR, config["output"]["filename"])
        # MUSIC overwrites input_powerspec.txt on every run, so keep this run's copy
        self.powerspec = np.loadtxt(os.path.join(WORKDIR, "input_powerspec.txt")) if self.returncode == 0 else None
        if "plug-in 'PANPHASIA' not found" in self.log:
            skip("MUSIC was built without PANPHASIA")
        m = re.search(r"shifted by \(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)", self.log)
        self.shift = np.array([int(x) for x in m.groups()]) if m else np.zeros(3, int)

    def ok(self):
        if self.returncode != 0:
            print(self.log[-3000:])
        return self.returncode == 0

    def level(self, level, field):
        """A generic-output field on one level without ghost cells, and the absolute index of its
        first cell in the unshifted box."""
        with h5py.File(self.output, "r") as f:
            data = f[f"level_{level:03d}_DM_{field}"][()][4:-4, 4:-4, 4:-4]
            lmin = int(f["header"].attrs["levelmin"])
            offset = np.zeros(3, int)
            for l in range(lmin + 1, level + 1):  # offsets are in cells of the parent level
                offset = 2 * (offset + np.array([f[f"header/grid_off_{a}"][()][l - lmin] for a in "xyz"]))
        return data, offset - self.shift * 2 ** (level - lmin)


def run_music(name, **sections):
    """Run MUSIC with BASE updated by the given sections; a value of None removes a key."""
    config = copy.deepcopy(BASE)
    config["output"]["filename"] = f"{name}.hdf5"
    for section, values in sections.items():
        config.setdefault(section, {}).update(values)
        config[section] = {k: v for k, v in config[section].items() if v is not None}
    return Run(name, config)


def cut(full, offset, shape):
    """Periodic sub-box of a full-box field."""
    n = full.shape[0]
    return full[np.ix_(*[(offset[k] + np.arange(shape[k])) % n for k in range(3)])]


def wavenumbers(shape):
    """|k| in units of the Nyquist wavenumber, for an rfftn of the given shape."""
    k = [np.fft.fftfreq(n) * 2 for n in shape[:2]] + [np.fft.rfftfreq(shape[2]) * 2]
    return np.sqrt(k[0][:, None, None] ** 2 + k[1][None, :, None] ** 2 + k[2][None, None, :] ** 2)


def cross_spectrum(a, b, bins, window=True):
    """Per |k| bin (units of Nyquist): correlation coefficient and amplitude ratio of a against b.
    Use the window for non-periodic sub-boxes."""
    if window:
        w = [np.hanning(n) for n in a.shape]
        w = w[0][:, None, None] * w[1][None, :, None] * w[2][None, None, :]
        a, b = a * w, b * w
    fa, fb = np.fft.rfftn(a - a.mean()), np.fft.rfftn(b - b.mean())
    k = wavenumbers(a.shape)
    out = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = (k > lo) & (k <= hi)
        x, y = fa[s], fb[s]
        r = np.real(np.sum(x * np.conj(y))) / np.sqrt(np.sum(abs(x) ** 2) * np.sum(abs(y) ** 2))
        out.append((lo, hi, r, np.sqrt(np.sum(abs(x) ** 2) / np.sum(abs(y) ** 2))))
    return out
