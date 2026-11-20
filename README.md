# Data and code: switching-layer thickness study of Ag/WO3/FTO valence-change memristors

Data and analysis code for the pre-print

> **A Physics-Based Multi-Objective Framework for Selecting the Switching-Layer Thickness of Ag/WO3/FTO Valence-Change Memristors: A Coupled Electrical-Ionic-Joule-Heating Simulation Study**
> Ayush Gupta (Department of Electronics and Communication Engineering, National Institute of Engineering, Mysuru) and
> Tribikram Gupta (Department of Physics, R. V. College of Engineering, Bengaluru).

The study is computational. A two-dimensional COMSOL Multiphysics 6.2 model (Electric Currents, Transport of Diluted
Species for oxygen vacancies, Heat Transfer) of an Ag/WO3/FTO valence-change memristor is used to map the saturated
HRS/LRS ratio, switching time, saturation time, switching energy and volumetric heat generation (`ec.Qh`) over WO3
thicknesses of 50-210 nm and voltages of 0.2-1.2 V, and to derive thickness windows for different design goals.
The results are simulation results and have not been validated experimentally.

## Contents

```
README.md                  this file
LICENSE.md                 licences (CC BY 4.0 for data and documentation, MIT for code)
data/                      exported COMSOL time series and per-run tables (about 51 MB)
  README.md                column definitions, run list, definitions of HRS/LRS, t_sat and t_sw
reproducibility/           Python scripts that regenerate every number and computed figure of the paper
  README.md                what each script does and which inputs are not recomputed
  memlib.py                shared functions
  01_ ... 09_*.py          analysis steps; 09 compares the regenerated numbers with the paper
  run_all.py               runs all steps
  requirements.txt         numpy, pandas, scipy, matplotlib
  results/                 reference outputs (CSV and text) of the scripts
  figures/                 computed figures of the paper
```

## Reproducing the results

COMSOL is not needed. Python 3.10 or later is required.

```
cd reproducibility
pip install -r requirements.txt
python run_all.py
```

The run takes about one minute, regenerates the tables, scores, thickness windows, weight-sensitivity and Monte Carlo
results, mesh comparison and computed figures, and ends with a PASS/FAIL list against the numbers printed in the paper.
Random seeds are fixed, so the Monte Carlo results are identical on every run. The reference outputs of a complete run
are in `reproducibility/results/`.

Two inputs are not recomputed: the switching times `t_sw`, which were selected manually from diagnostic plots of dR/dt
(an automated rule is provided as a cross-check), and the noise levels used in the Monte Carlo propagation (see
`reproducibility/README.md`).

## Simulated runs

| Mesh | Sweep | Thicknesses (nm) | Voltages (V) | Run length |
|---|---|---|---|---|
| Extra Fine | voltage sweep | 50, 70, 90, 110, 130, 150, 170, 190, 210 | 0.2-1.2 in steps of 0.1 | 15 s |
| Extra Fine | refined sweep | 80-110 in steps of 1 | 1 | 12 s (80-85 nm), 15 s (86-110 nm) |
| Extremely Fine | refined sweep | 80-110 in steps of 1 | 1 | 12 s |
| Extra Fine | long runs (resistance only) | 110, 130, 150, 170, 190, 210 | 1 | 20 s; 40 s for 150, 170, 190, 210 |

Units: time in s, resistance in ohm, `ec.Qh` in W/m^3, thickness in nm, voltage in V. Output interval 0.01 s.

## COMSOL model files

The COMSOL model files (.mph) are not included in this record. They are available from the authors on request.
The model is described in Section 2 and Section 3 of the paper (geometry, materials, governing equations, boundary
conditions, mesh presets and study settings), and the exported data in `data/` are the inputs of all analyses.

## Citation

If you use these data or code, please cite the paper and this record (DOI assigned by Zenodo on publication).

## Licence

Data, figures, result tables and documentation: CC BY 4.0. Source code: MIT. See `LICENSE.md`.
