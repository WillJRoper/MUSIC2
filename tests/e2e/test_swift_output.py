# ctest-labels: fast
"""SWIFT initial conditions: total mass and baryon fraction, particles inside the box, a valid cell
structure, and masked particle loads that conserve the mass of every particle type."""

import os
import sys

import h5py
import numpy as np

from musictest import check, finish, run_music

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import check_swift_cells  # noqa: E402
import check_swift_masked  # noqa: E402

# small boxes at high redshift so that the masked-load checker can find each block's gas cell
SETUP = {"boxlength": 100, "zstart": 1000, "levelmin": 5, "levelmin_TF": 5, "levelmax": 6, "baryons": "yes",
         "ref_center": "0.5, 0.5, 0.5", "ref_extent": "0.3, 0.3, 0.3"}
SWIFT = {"format": "swift"}


def passes(what, function, *args):
    try:
        function(*args)
    except AssertionError as err:
        return check(what, False, f"assertion failed: {err}")
    return check(what, True)


def check_masses(run, path):
    with h5py.File(path, "r") as f:
        units, cosmo = dict(f["Units"].attrs), dict(f["Cosmology"].attrs)
        box = f["Header"].attrs["BoxSize"][0]
        mass = {k: f[f"{k}/Masses"][:].sum() for k in f if k.startswith("PartType")}
        inside = all(np.all((f[f"{k}/Coordinates"][:] >= 0) & (f[f"{k}/Coordinates"][:] < box)) for k in mass)
        levels = len(f["Header"].attrs["Music_levelcounts"])
    h0 = cosmo["H0 [internal units]"][0] * 1e5 / 3.08567758e24  # km/s/Mpc in 1/s
    rho_crit = 3 * h0**2 / (8 * np.pi * 6.6743e-8) * units["Unit length in cgs (U_L)"][0] ** 3 / units["Unit mass in cgs (U_M)"][0]
    want = cosmo["Omega_m"][0] * rho_crit * box**3
    total = sum(mass.values())
    check(f"{run.name}: total mass is Omega_m rho_crit V", abs(total / want - 1) < 1e-3, f"{total:.6g} vs {want:.6g}")
    if levels == 1:  # baryon fraction only holds everywhere without coarse total-matter particles
        fb = mass["PartType0"] / total
        want_fb = cosmo["Omega_b"][0] / cosmo["Omega_m"][0]
        check(f"{run.name}: gas fraction is Omega_b/Omega_m", abs(fb / want_fb - 1) < 1e-6, f"{fb:.6f}")
    check(f"{run.name}: all particles inside the box", inside)


uniform = run_music("swift_uniform", setup={**SETUP, "levelmin": 6, "levelmin_TF": 6}, output=SWIFT)
zoom = run_music("swift_zoom", setup=SETUP, output=SWIFT)
for run in (uniform, zoom):
    if check(f"{run.name}: MUSIC runs", run.ok()):
        check_masses(run, run.output)
        passes(f"{run.name}: cell structure", check_swift_cells.check, run.output)

# a ParameterSet gives the cosmology without Omega_m etc. in the config file
pset = run_music("swift_pset", setup={**SETUP, "levelmin": 6, "levelmin_TF": 6}, output=SWIFT,
                 cosmology={"ParameterSet": "Planck2018EE+BAO+SN", "Omega_m": None, "Omega_L": None,
                            "Omega_b": None, "H0": None, "n_s": None})
if check(f"{pset.name}: MUSIC runs", pset.ok()):
    check_masses(pset, pset.output)
    with h5py.File(pset.output, "r") as f:
        om, h = f["Header"].attrs["Omega0"], f["Header"].attrs["HubbleParam"]
    check(f"{pset.name}: header has the set's cosmology", np.isclose(om, 0.3099) and np.isclose(h, 0.67742), f"{om}, {h}")

for mask_type in range(4):
    run = run_music(f"swift_masked{mask_type}", setup={**SETUP, "ParticleLoad": "masked", "ParticleMaskType": mask_type},
                    output=SWIFT)
    if check(f"{run.name}: MUSIC runs", run.ok()):
        check_masses(run, run.output)
        passes(f"{run.name}: cell structure", check_swift_cells.check, run.output)
        if zoom.returncode == 0:
            passes(f"{run.name}: masked load", check_swift_masked.check, run.output, mask_type, zoom.output)

finish()
