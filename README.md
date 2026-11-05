# Memristors_Thickness_study

Data and scripts for the paper

> **A Physics-Based Multi-Objective Framework for Selecting the Switching-Layer Thickness of Ag/WO3/FTO Valence-Change Memristors: A Coupled Electrothermal Simulation Study**
> Ayush Gupta (Department of Electronics and Communication Engineering, National Institute of Engineering, Mysuru) and Tribikram Gupta (Department of Physics, R. V. College of Engineering, Bangalore).

The study uses a two-dimensional COMSOL Multiphysics model of an Ag/WO3/FTO valence-change memristor (electric currents + oxygen-vacancy drift-diffusion + Joule heating) to map the HRS/LRS ratio, switching time, switching energy and volumetric heat generation (`ec.Qh`) over WO3 thicknesses of 50-210 nm and voltages of 0.2-1.2 V, and to derive thickness windows for different design goals. The device architecture follows the fabricated Ag/WO3/FTO devices of K. Hadiyal, N. Raghavan and R. M. Thamankar, *Adv. Electron. Mater.* **12**, e00658 (2026), doi:10.1002/aelm.202500658. The results are simulation results and have not been validated experimentally.

## Repository contents

| Item | What it is |
|---|---|
| `resistance_all_thicknesses_timeseries.csv` | Terminal resistance versus time for every simulated run |
| `resistance_all_thicknesses_summary.csv` | One row per run: initial/end/minimum resistance, HRS/LRS, switching time, saturation time (where read) |
| `ecQh_all_thicknesses_timeseries.csv` | `ec.Qh` (surface average over Ag + WO3 + FTO) versus time for every run |
| `ecQh_all_thicknesses_summary.csv` | One row per run: mean over the full run and over 0-12 s, peak value |
| `vapp.py` | Plots resistance versus time for each thickness and voltage |
| `vapp_2_find_st_.py` | Extracts switching time, HRS/LRS, switching energy and average power from resistance data |
| `vapp_3_switching_energy_.py` | Recomputes switching energy using manually corrected switching times |
| `wkb_tunneling_comsol.py` | WKB transmission and Tsu-Esaki tunnelling current from 1-D COMSOL potential profiles |
| `run_ag_schottky_reset.py` | Wrapper that applies the WKB tool to the RESET polarity (Ag at +V, FTO grounded) |

The COMSOL model file itself is not part of this repository unless it is listed in the table above.

## Data

All times are in seconds, resistance in ohm, `ec.Qh` in W/m^3, thickness in nm, voltage in V. Output was written every 0.01 s.

### Simulated runs

| Mesh | Sweep | Thicknesses (nm) | Voltages (V) | Run length |
|---|---|---|---|---|
| Extra Fine | voltage sweep | 50, 70, 90, 110, 130, 150, 170, 190, 210 | 0.2-1.2 in steps of 0.1 | 15 s |
| Extra Fine | refined sweep | 80-85 | 1 | 12 s |
| Extra Fine | refined sweep | 86-95 | 1 | 15 s |
| Extra Fine | refined sweep | 96-110 | 1 | 15 s |
| Extremely Fine | refined sweep | 80-85 | 1 | 12 s |
| Extremely Fine | refined sweep | 86-95 | 1 | 12 s |
| Extremely Fine | refined sweep | 96-110 | 1 | 12 s |
| Extra Fine | long runs (resistance only) | 110, 130, 150, 170, 190, 210 | 1 (see note) | 20 s; 40 s for 150, 170, 190, 210 |

All Extremely Fine runs come from raw exports (the original export columns repeat their last value after 12 s; those repeated values are omitted).

### Columns of the time-series files

`mesh`, `sweep`, `thickness_nm`, `voltage_V`, `time_s`, then `resistance_ohm` (resistance files) or `ecQh_W_per_m3` and `terminal_current_A_per_m2` (ec.Qh files; the current is empty where it was not exported), and `source_file` (original export). The resistance file also has `voltage_note`.

### Notes

- **Averaging domain.** `ec.Qh` is the COMSOL surface average over all three domains (Ag, FTO and WO3). Multiply by (d + 40 nm) / d to estimate the WO3-only mean if electrode heating is negligible. The 110-210 nm voltage-sweep values are a re-export over all three domains; the earlier export covered WO3 only and is not included.
- **Window.** The full-run mean of `ec.Qh` is a 15 s mean for the voltage sweep and the Extra Fine 96-110 nm runs, and a 12 s mean for the other runs. Use `ecQh_mean_0_12s` in the summary file for a common 12 s window.
- **HRS/LRS** in the summary file is `R_end / R_initial`, so it depends on the run length; thick devices are still rising at the end of a 15 s run.
- **Long runs.** The original long-run exports have no voltage column; they are listed at 1 V because their first 15 s agree with the 1 V runs to a median of 0.4-0.8 %. The 210 nm device has both a 20 s and a 40 s run.
- **Switching times** are the maximum of dR/dt of the main early peak, read from diagnostic plots. `switching_time_source` in the summary file says whether a value was `author hand-picked` or came from the `automatic rule` (median of nine variants of "maximum of dR/dt before R has completed 10, 20 or 30 % of its rise"). Solver noise limits them to about +/-0.01-0.02 s. The literal global maximum of dR/dt is often a late solver jump, which is why plain `np.argmax` is not used for the final values.
- **Saturation times** (`saturation_time_author_s`) are read from the long runs as the time at which the resistance levels off, not from a numerical threshold.
- 110 nm appears in more than one sweep; the `sweep` column tells them apart.

## Scripts

Dependencies: Python 3.9+, `numpy`, `pandas`, `matplotlib` (the WKB scripts use only the standard library). The scripts have hard-coded input file names; edit them at the top of each script.

**NumPy 2.x:** `np.trapz` was removed. Replace it with `np.trapezoid`, or use NumPy < 2.

### Resistance post-processing

COMSOL export format: five header lines, then four columns `Thickness, Voltage, Time, Resistance`.

1. `vapp.py` - plots R(t) for each thickness, one curve per voltage (reads `ultra_fine_86nm_sweep_resis.csv` as written; change the file name).
2. `vapp_2_find_st_.py` - for every thickness-voltage pair: switching time = time of maximum dR/dt (`np.gradient`), LRS = R(first sample), HRS = R(last sample), HRS/LRS, switching energy = integral of V^2/R dt up to the switching point (trapezoidal rule), average power = energy / switching time. Shows a diagnostic plot of R(t) and dR/dt with the marked point and writes `SwitchingResults.csv`.
3. `vapp_3_switching_energy_.py` - reads the resistance CSV and the (edited) `SwitchingResults.csv`, and recomputes the switching energy with those switching times. Writes `SwitchingEnergyResults.csv`.

Workflow used in the study: run step 2, inspect each diagnostic plot, correct the `Switching Time (s)` column of `SwitchingResults.csv` by hand where the global maximum is a late solver jump, then run step 3. The energy integral stops at the last sample at or before the chosen switching time.

### Tunnelling estimate (one-way post-processing)

A preliminary estimate of tunnelling current through the Ag/WO3 barrier from the classical COMSOL potential. It does not feed back into the drift-diffusion-heating model, so it is not a self-consistent quantum-transport calculation.

- `wkb_tunneling_comsol.py` reads a COMSOL text export of 1-D potential cuts (arc length in nm, potential in V; the arc length restarts at 0 for each stored time step), computes the WKB transmission and the 3-D Tsu-Esaki current density, and writes `wkb_tunnelling_current.csv` with columns `solution_index, time_s, points, thickness_nm, bias_v, wkb_integral_EF_sqrt_eV_nm, wkb_exponent_EF, transmission_EF, current_density_A_per_m2, current_A`. Default settings: barrier height 0.56 eV, tunnelling effective mass 0.45 m0, transverse supply mass 0.5 m0, T = 293 K, device area 2.5e-13 m^2 (500 nm x 500 nm), time step 0.01 s, 801 energy points. Set `INPUT_FILE` to your export.
- `run_ag_schottky_reset.py` imports `wkb_tunneling_comsol.py` and replaces the supply function for the RESET polarity (FTO at 0 V, Ag at +V): the barrier is U(x, E) = Phi_Ag/WO3 - V(x) - E (eV), referenced to the FTO Fermi level. The COMSOL cut must contain only WO3 and be used in its native order (x = 0 at the WO3/FTO interface, x = L at the WO3/Ag interface); do not reverse it. Run `python run_ag_schottky_reset.py` with both files in the same folder.

## Citation

Please cite the paper (reference to be added when available) and, for the device architecture, Hadiyal et al., *Adv. Electron. Mater.* **12**, e00658 (2026).

## Licence

No licence file has been added yet. Add one before reuse (for example CC BY 4.0 for the data and the MIT licence for the code).
