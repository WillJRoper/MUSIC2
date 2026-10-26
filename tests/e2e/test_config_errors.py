# ctest-labels: fast
"""Invalid configurations stop MUSIC with a clear message instead of producing wrong output."""

from musictest import FLAMINGO, check, finish, run_music

PANPH6 = {"generator": "PANPHASIA", "descriptor": FLAMINGO}
CASES = {
    "bcc_swift": ({"setup": {"ParticleLoad": "bcc", "baryons": "yes"}, "output": {"format": "swift"}},
                  "ParticleLoad must be 'sc' or 'masked'"),
    "masked_generic": ({"setup": {"ParticleLoad": "masked", "baryons": "yes"}},
                       "only supported for output format 'swift'"),
    "mask_type": ({"setup": {"ParticleLoad": "masked", "ParticleMaskType": 7, "baryons": "yes"}, "output": {"format": "swift"}},
                  "ParticleMaskType must be 0, 1, 2 or 3"),
    "odd_ref_dims": ({"setup": {"levelmax": 8, "ref_extent": None, "ref_dims": "51, 52, 52"}},
                     "ref_dims must be even"),
    "kk_levelmin": ({"setup": {"levelmin": 6, "levelmin_TF": 6, "levelmax": 6}, "random": PANPH6},
                    "levelmin is too small for the descriptor's KK"),
    "kk_fixing": ({"setup": {"fix_mode_amplitude": "yes"}, "random": PANPH6},
                  "fix_mode_amplitude cannot be combined with a KK descriptor"),
    "base_size": ({"random": {"generator": "PANPHASIA", "descriptor": "[Panph6,L18,(56034,71400,250000),S3,CH-999,test]"}},
                  "power-of-two base size"),
    "check_digit": ({"random": {"generator": "PANPHASIA", "descriptor": FLAMINGO.replace("CH3774755196", "CH3774755197")}},
                    "invalid descriptor"),
}

for name, (sections, message) in CASES.items():
    run = run_music(f"err_{name}", **sections)
    check(f"{name} is rejected", run.returncode != 0 and message in run.log,
          f"exit {run.returncode}" + ("" if message in run.log else f", no '{message}' in the log"))

finish()
