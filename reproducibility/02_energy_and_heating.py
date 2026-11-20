"""Step 2. Switching energy, heat-generation metrics and the consistency checks of Sec. 4.2-4.3.

  * E_sw = integral_0^t_sw V^2/R dt from the raw R(t) with the manually selected t_sw
    (compared with the energies tabulated in the paper);
  * mean and peak ec.Qh (full-stack surface average), mean total Joule power
    <P> = <Qh>*(d+40 nm)*500 nm*500 nm, slab temperature estimate dT = Q_ox d^2/(8k);
  * check that P from ec.Qh equals V^2/R (all voltages, all times);
  * electrode vertical resistance and Peclet-number estimates quoted in Sec. 2.2 / 4.3.

Outputs: results/energy_heating_coarse_1V.csv, results/energy_heating_refined.csv,
         results/power_consistency_check.csv
"""
import numpy as np, pandas as pd
import memlib as m

A_M2 = (m.AREA_NM * 1e-9) ** 2
W_M = m.AREA_NM * 1e-9

QH_SWEEP = {  # sweep names of the ec.Qh file for the coarse nodes at 1 V
    50: "voltage sweep (50-90 nm)", 70: "voltage sweep (50-90 nm)", 90: "voltage sweep (50-90 nm)"}
def qh_sweep_coarse(th):
    return QH_SWEEP.get(th, "voltage sweep (110-210 nm, full-stack re-export)")

def stack_volume_m3(th_nm):
    return (th_nm + m.STACK_EXTRA_NM) * 1e-9 * W_M * W_M

# ---------------------------------------------------------------- coarse sweep, 1 V
C = m.load_coarse_table()
rows = []
for th in m.COARSE_NODES:
    t, r = m.coarse_run(th)
    t15, r15 = m.get_run(m.resistance(), "Extra Fine", m.VS, th)
    tq, q = m.get_run(m.ecqh(), "Extra Fine", qh_sweep_coarse(th), th, col="ecQh_W_per_m3")
    tsw = C.t_sw_selected_s[th]
    e_raw = m.switching_energy(t15, r15, tsw)
    t90 = C.t90_s[th]
    qmean, qpeak = q.mean(), q.max()
    p_mean = qmean * stack_volume_m3(th)
    q_ox = qpeak * (th + m.STACK_EXTRA_NM) / th
    dT = q_ox * (th * 1e-9) ** 2 / (8 * m.K_WO3)
    rows.append(dict(thickness_nm=th, t_sw_selected_s=tsw, E_sw_raw_J=e_raw, E_sw_table_J=C.E_sw_J[th],
                     E_sw_rel_diff=e_raw / C.E_sw_J[th] - 1,
                     E90_J=m.switching_energy(t15, r15, t90) if t90 <= t15[-1] else np.nan,
                     mean_ecQh=qmean, mean_ecQh_table=C.mean_ecQh_full_stack_W_per_m3[th],
                     peak_ecQh=qpeak, peak_first_sample_frac=q[0] / qpeak,
                     mean_total_power_uW=p_mean * 1e6, dT_max_K=dT))
out = pd.DataFrame(rows).set_index("thickness_nm")
out.to_csv(m.RESULTS / "energy_heating_coarse_1V.csv")
print(out[["E_sw_raw_J", "E_sw_table_J", "E_sw_rel_diff", "mean_ecQh", "mean_ecQh_table",
           "mean_total_power_uW", "dT_max_K"]].to_string(float_format=lambda x: f"{x:.4g}"))
print("power ratio 50->210 nm: %.2f ; volumetric mean ratio: %.2f" %
      (out.mean_total_power_uW[50] / out.mean_total_power_uW[210], out.mean_ecQh[50] / out.mean_ecQh[210]))

# ---------------------------------------------------------------- refined sweeps, 1 V
R = m.load_refined_table()
rows = []
for _, p in R.iterrows():
    mesh, th = p.mesh, int(p.thickness_nm)
    t, r = m.refined_run(mesh, th)
    e_raw = m.switching_energy(t, r, p.switching_time_s)
    qs = {"Extra Fine": {True: "refined sweep 80-85 nm", False: None}}.get(mesh)
    sweep = ("refined sweep 80-85 nm" if th < 86 else
             ("refined sweep 86-95 nm (15 s)" if mesh == "Extra Fine" else "refined sweep 86-95 nm (12 s)")
             if th < 96 else "refined sweep 96-110 nm")
    if mesh == "Extremely Fine" and th < 86:
        sweep = "refined sweep 80-85 nm (12 s)"
    tq, q = m.get_run(m.ecqh(), mesh, sweep, th, col="ecQh_W_per_m3")
    q12 = q[tq <= 12 + 1e-9].mean()            # common 12 s basis
    rows.append(dict(mesh=mesh, thickness_nm=th, E_sw_raw_J=e_raw, E_sw_table_J=p.switching_energy_J,
                     E_sw_rel_diff=e_raw / p.switching_energy_J - 1,
                     mean_ecQh_12s=q12, mean_ecQh_table=p.mean_ecQh_W_per_m3,
                     mean_ecQh_run=q.mean()))
ref = pd.DataFrame(rows)
ref.to_csv(m.RESULTS / "energy_heating_refined.csv", index=False)
print("\nrefined E_sw: |raw/table-1| <=0.4%% for %d of %d runs, <=1%% for %d, max %.2f%%" % (
    (ref.E_sw_rel_diff.abs() <= 0.004).sum(), len(ref), (ref.E_sw_rel_diff.abs() <= 0.01).sum(),
    100 * ref.E_sw_rel_diff.abs().max()))
print("refined mean Qh (12 s) vs table: max rel diff %.3f%%" % (100 * (ref.mean_ecQh_12s / ref.mean_ecQh_table - 1).abs().max()))

# ---------------------------------------------------------------- P = V^2/R consistency
res, qh = m.resistance(), m.ecqh()
pairs = [  # (mesh, qh sweep, resistance sweep) with matching thickness sets
    ("Extra Fine", "voltage sweep (50-90 nm)", m.VS, (50, 70, 90)),
    ("Extra Fine", "voltage sweep (110-210 nm, full-stack re-export)", m.VS, (110, 130, 150, 170, 190, 210)),
]
recs, allr = [], []
for mesh, qs, rs, ths in pairs:
    for th in ths:
        for v in sorted(res[(res.sweep == rs) & (res.thickness_nm == th)].voltage_V.unique()):
            t1, r1 = m.get_run(res, mesh, rs, th, v)
            t2, q2 = m.get_run(qh, mesh, qs, th, v, col="ecQh_W_per_m3")
            n = min(len(t1), len(t2))
            ratio = q2[:n] * stack_volume_m3(th) / (v ** 2 / r1[:n])
            recs.append(dict(thickness_nm=th, voltage_V=v, median_ratio=np.median(ratio), ratio_t0=ratio[0],
                             ratio_min=ratio.min(), ratio_max=ratio.max()))
            allr.append(pd.DataFrame(dict(v=v, th=th, ratio=ratio)))
pc = pd.DataFrame(recs)
pc.to_csv(m.RESULTS / "power_consistency_check.csv", index=False)
print("\nP(ec.Qh)/(V^2/R) over %d runs: median of medians %.4f, range of medians %.4f-%.4f; t=0 within %.2f%%; "
      "extremes %+.1f%% to %+.1f%%" % (len(pc), pc.median_ratio.median(), pc.median_ratio.min(), pc.median_ratio.max(),
                                        100 * (pc.ratio_t0 - 1).abs().max(), 100 * (pc.ratio_min.min() - 1),
                                        100 * (pc.ratio_max.max() - 1)))

allr = pd.concat(allr)
print("median ratio over all samples of all runs: %.4f ; per voltage (min..max of 11): %.4f..%.4f ; "
      "per thickness (9): %.4f..%.4f ; 1st-99th percentile of all samples: %.3f..%.3f" % (
          allr.ratio.median(), allr.groupby("v").ratio.median().min(), allr.groupby("v").ratio.median().max(),
          allr.groupby("th").ratio.median().min(), allr.groupby("th").ratio.median().max(),
          allr.ratio.quantile(.01), allr.ratio.quantile(.99)))

# ---------------------------------------------------------------- electrodes and Peclet number
sig_fto, sig_ag = 5e5, 6.3e7        # S/m (FTO: Table 1; Ag: bulk value near 293 K, COMSOL built-in is similar)
t_el = 20e-9
R_fto, R_ag = t_el / (sig_fto * A_M2), t_el / (sig_ag * A_M2)
Rlrs_min = min(m.coarse_run(th)[1][0] for th in m.COARSE_NODES)
print("\nelectrode vertical resistance: FTO %.3f ohm, Ag %.2f mohm; R_LRS,min = %.0f ohm; "
      "electrode share of Joule power <= %.1e" % (R_fto, 1e3 * R_ag, Rlrs_min, (R_fto + R_ag) / Rlrs_min))
kT = 8.617e-5 * 293.15                # eV
pe_ox = 2 * 1.0 / kT                  # z=2, 1 V across the oxide: Pe = z e V / (k T)
d, h = 100.0, 1.0
print("Peclet: oxide (1 V, z=2) Pe = %.0f ; mesh-cell Pe_h = z e V h /(kT d) = h/(%.2f nm) at d = 100 nm"
      % (pe_ox, kT * d / 2.0))
