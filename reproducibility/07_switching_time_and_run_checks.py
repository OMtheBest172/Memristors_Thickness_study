"""Step 7. Checks quoted in Sec. 4.1 and 3.3: behaviour of R(t) (initial dip), the selection of the switching
time (global maximum of dR/dt vs the selected early peak vs the automated rule), run-to-run differences at the thicknesses that appear in two sweeps, the
position of the heating maximum.  Output: results/run_checks.txt
"""
import io
import numpy as np, pandas as pd
import memlib as m

buf = io.StringIO()
def say(*a):
    print(*a); print(*a, file=buf)

res, qh = m.resistance(), m.ecqh()
summ = pd.read_csv(m.DATA / "resistance_all_thicknesses_summary.csv")
vs = summ[summ.sweep == m.VS]

# ---- initial dip of R(t) at 1 V (all runs at 1 V in the per-run summary file)
one = summ[np.isclose(summ.voltage_V, 1.0)]
dip = 100 * (1 - one.R_min_ohm / one.R_initial_ohm)
say(f"R(t) dip below R(0) at 1 V over {len(one)} runs: up to {dip.max():.1f} %; the minimum occurs at "
    f"t = {one.time_of_R_min_s.min():.2f}-{one.time_of_R_min_s.max():.2f} s")

# ---- selection of t_sw in the 99 voltage-sweep runs
rows = []
for _, p in vs.iterrows():
    t, r = m.get_run(res, "Extra Fine", m.VS, int(p.thickness_nm), p.voltage_V)
    dy = np.gradient(r, t)
    s = m.smooth(r); ref = np.median(s[t >= t[-1] - 1]); r0 = r[0]
    idx = np.where((s - r0) / (ref - r0) >= 0.2)[0]
    end = max(idx[0] if len(idx) else len(t) - 1, 2)
    rows.append(dict(th=int(p.thickness_nm), v=p.voltage_V, sel=p.switching_time_s,
                     glob=t[np.argmax(dy)], lim20=t[np.argmax(dy[:end + 1])],
                     auto=m.tsw_automated(t, r)))
V = pd.DataFrame(rows)
V = V[V.sel.notna()]
n = len(V)
same = lambda a: np.abs(V[a] - V.sel) < 0.005
say(f"voltage-sweep runs with a selected t_sw: {n}")
say(f"  global argmax(dR/dt) differs from the selected value in {int((~same('glob')).sum())} runs")
say(f"  argmax restricted to the first 20 % of the rise reproduces it in {int(same('lim20').sum())} runs; "
    f"the other {int((~same('lim20')).sum())} by thickness: {V[~same('lim20')].th.value_counts().sort_index().to_dict()}")
say(f"  automated rule within one output step (0.01 s): {int((np.abs(V.auto - V.sel) <= 0.0101).sum())} of {n}")
V.to_csv(m.RESULTS / "tsw_selection_voltage_sweep.csv", index=False)

# ---- automated rule vs selected value, refined sweeps and coarse 1 V
S1 = pd.read_csv(m.RESULTS / "saturation_metrics_1V.csv")
R = m.load_refined_table()
refd = S1[S1.set == "refined"].merge(R[["mesh", "thickness_nm", "switching_time_s"]], on=["mesh", "thickness_nm"])
dd = np.abs(refd.tsw_auto - refd.switching_time_s)
say(f"refined sweeps: automated rule within 0.01 s for {int((dd <= 0.0101).sum())} of {len(dd)} runs, identical for {100*(dd < 1e-3).mean():.0f} %")
C = m.load_coarse_table()
dc = (S1[S1.set == "coarse"].set_index("thickness_nm").tsw_auto - C.t_sw_selected_s)
say(f"coarse 1 V: automated rule within 0.02 s for {int((dc.abs() <= 0.0201).sum())} of 9; larger differences: "
    f"{ {int(k): (round(float(S1[S1.set=='coarse'].set_index('thickness_nm').tsw_auto[k]),2), float(C.t_sw_selected_s[k])) for k in dc.index[dc.abs() > 0.0201]} } (auto, selected)")

cs = S1[S1.set == "coarse"].set_index("thickness_nm")

# ---- run-to-run differences between coarse and refined runs at 90 and 110 nm
EF = R[R.mesh == "Extra Fine"].set_index("thickness_nm")
ef1 = S1[(S1.set == "refined") & (S1.mesh == "Extra Fine")].set_index("thickness_nm")
for th in (90, 110):
    say(f"{th} nm coarse vs refined (EF): t_sw {C.t_sw_selected_s[th]:.2f} vs {EF.switching_time_s[th]:.2f} s; "
        f"E_sw {100*(EF.switching_energy_J[th]/C.E_sw_J[th]-1):+.1f} %; Lambda {100*(ef1.lam[th]/cs.lam[th]-1):+.2f} %; "
        f"t_sat {cs.t_sat_995[th]:.2f} vs {ef1.t_sat_995[th]:.2f} s")

# ---- heating: position of the maximum (voltage-sweep runs)
first = []
for qs, ths in (("voltage sweep (50-90 nm)", (50, 70, 90)), ("voltage sweep (110-210 nm, full-stack re-export)", (110, 130, 150, 170, 190, 210))):
    for th in ths:
        for v in sorted(qh[(qh.sweep == qs) & (qh.thickness_nm == th)].voltage_V.unique()):
            t, q = m.get_run(qh, "Extra Fine", qs, th, v, col="ecQh_W_per_m3")
            first.append(100 * q[0] / q.max())
allq = qh.groupby(["mesh", "sweep", "thickness_nm", "voltage_V"]).ecQh_W_per_m3.agg(first="first", peak="max")
say(f"\nec.Qh(0) as a fraction of the peak over all {len(allq)} runs: {100*(allq['first']/allq.peak).min():.0f}-{100*(allq['first']/allq.peak).max():.0f} %")

(m.RESULTS / "run_checks.txt").write_text(buf.getvalue())
