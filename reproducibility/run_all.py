"""Run every step in order and finish with the comparison against the paper.
Usage:  python run_all.py            (about 1-2 minutes; needs ../data)
        python run_all.py --paper    (also copies the regenerated figures to ../fig)
"""
import subprocess, sys, time
from pathlib import Path

here = Path(__file__).resolve().parent
steps = ["01_saturation_metrics.py", "02_energy_and_heating.py", "03_coarse_scoring_and_windows.py",
         "04_weight_sensitivity.py", "05_noise_monte_carlo.py", "06_refined_sweeps_and_mesh.py",
         "07_switching_time_and_run_checks.py", "08_make_figures.py", "09_check_paper_numbers.py"]
extra = ["--paper"] if "--paper" in sys.argv else []
for s in steps:
    t0 = time.time()
    print(f"\n{'=' * 78}\n{s}\n{'=' * 78}", flush=True)
    r = subprocess.run([sys.executable, str(here / s)] + (extra if s.startswith("08") else []), cwd=here)
    print(f"[{s}: exit {r.returncode}, {time.time() - t0:.0f} s]")
    if r.returncode != 0 and not s.startswith("09"):
        sys.exit(r.returncode)
