# Reproducibility package

Regenerates every table value, score, thickness window, weight-sensitivity and Monte Carlo result,
mesh comparison and computed figure of the paper from the exported COMSOL time series in `../data`.
COMSOL is **not** needed. Runtime about one minute.

```
pip install -r requirements.txt          # numpy, pandas, scipy, matplotlib (Python 3.10+)
python run_all.py                        # runs steps 01-09, ends with a PASS/FAIL list against the paper
python run_all.py --paper                # additionally copies the regenerated figures to ../fig
```
The data folder is found automatically (`../data`, the package root, or `$MEM_DATA`).
Random seeds are fixed, so the Monte Carlo numbers are identical on every run.

## Steps

| script | what it computes | where it appears in the paper |
|---|---|---|
| `memlib.py` | shared functions: saturation criterion, switching energy, automated t_sw rule, min-max scoring, Pareto test, windows | Sec. 4 |
| `01_saturation_metrics.py` | R_LRS, saturated R_HRS, Lambda, t_sat (99.5 % and 99 %), t50, t90, E90, automated t_sw for all 1 V runs | Sec. 4.1 |
| `02_energy_and_heating.py` | E_sw from raw R(t); mean/peak ec.Qh; total power; slab temperature estimate; P = V^2/R check; electrode resistance; Peclet numbers | Sec. 2.2, 4.2, 4.3, 5.3 |
| `03_coarse_scoring_and_windows.py` | score variants C, D, N, T, U; Pareto test; correlations; windows at 1 V and (t_sw, E_sw, heating) common to all 11 voltages | Sec. 6.1, 6.4 |
| `04_weight_sensitivity.py` | 20 000 Dirichlet(1) weight vectors; leave-one-objective-out margins | Sec. 6.2 |
| `05_noise_monte_carlo.py` | propagation of run-to-run scatter (20 000 draws) | Sec. 6.2 |
| `06_refined_sweeps_and_mesh.py` | Extra Fine vs Extremely Fine, saturation statistics, quadratic trends, ranking under noise | Sec. 5.4, 6.3 |
| `07_switching_time_and_run_checks.py` | initial dip of R(t), global vs selected max of dR/dt, automated t_sw rule, coarse-vs-refined run-to-run differences | Sec. 3.3, 4.1, 4.3 |
| `08_make_figures.py` | the computed figures (not the schematic or the COMSOL vacancy images) | figures |
| `09_check_paper_numbers.py` | compares the regenerated numbers with those printed in the paper | - |

Outputs go to `results/` (CSV and text summaries) and `figures/`.

## Inputs that are not recomputed

* **Switching times t_sw** were selected manually from diagnostic plots of dR/dt (Sec. 4.1). They are read from
  `coarse_sweep_1V_values_used_in_paper.csv` (`t_sw_selected_s`), `refined_sweep_80_110nm_values_used_in_paper.csv`
  (`switching_time_s`) and `resistance_all_thicknesses_summary.csv` (`switching_time_s`, voltage sweep). The automated
  rule in `memlib.tsw_automated` reproduces them for all 62 refined runs within 0.01 s, but not for 130 and 190 nm in the
  coarse sweep (step 07 lists the cases).
* **Tabulated E_sw** (used in the scores) comes from the same CSV files. Step 02 recomputes it from the raw R(t):
  the agreement is within 0.4 % for all 9 coarse and all 62 refined runs.
* **Noise levels of the Monte Carlo** (Lambda 0.5 %, t_sw 0.015 s, t_sat 1 s, E_sw 3 %, heating 0.5 %) are input
  assumptions taken from the run-to-run spread and the refined-sweep residuals (paper Sec. 6.2); they are set at the top of
  `05_noise_monte_carlo.py`. They describe run-to-run scatter only, not the mesh dependence (step 06).
* **Ag conductivity** used for the electrode-resistance estimate in step 02 is the bulk value 6.3e7 S/m (COMSOL's
  built-in silver is temperature dependent and gives the same order of magnitude).

## What cannot be reproduced from these files

* The COMSOL simulations themselves (they need the .mph model and a COMSOL licence; the model files are available from the authors on request).
* The statement that the t_sat scatter of the Extra Fine sweep is caused by step-like changes of 0.3-0.5 % of the
  plateau (inspected visually); the dependence on the smoothing window is tested in step 06.
