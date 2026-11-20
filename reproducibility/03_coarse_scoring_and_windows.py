"""Step 3. Equal-weight scoring of the coarse sweep (50-170 nm, 1 V), variants C, D, N, T, U,
Pareto test, objective correlations, and the design windows of Sec. 6.4 / Table 7 / Fig. 10.

Variants (all with the saturated Lambda maximised):
  C  t_sw, E_sw, heating minimised, t_sat maximised (primary)
  D  as C with t_sat minimised
  N  as C without t_sat
  T  t_90 and E_90 in place of t_sw and E_sw, t_sat maximised
  U  as T without t_sat
Inputs: data/*values_used_in_paper.csv (manually selected t_sw, tabulated E_sw), results of step 1.
Outputs: results/scores_coarse.csv, results/windows_coarse_1V.csv, results/windows_all_voltages.csv
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import memlib as m

C = m.load_coarse_table()
S1 = pd.read_csv(m.RESULTS / "saturation_metrics_1V.csv")
S1 = S1[S1.set == "coarse"].set_index("thickness_nm")
H = pd.read_csv(m.RESULTS / "energy_heating_coarse_1V.csv", index_col=0)

D = pd.DataFrame({
    "lam": S1.lam, "reached": S1.reached, "tsw": C.t_sw_selected_s, "tsat": S1.t_sat_995,
    "t90": S1.t90, "E": C.E_sw_J, "E90": S1.E90_J, "Q": C.mean_ecQh_full_stack_W_per_m3,
    "P": H.mean_total_power_uW, "dT": H.dT_max_K})
print("E90 (script) vs paper table E_to_t90_J: max rel diff %.2e" % (D.E90 / C.E_to_t90_J - 1).abs().max())

VARIANTS = {
    "C": {"lam": 1, "tsw": 0, "tsat": 1, "E": 0, "Q": 0},
    "D": {"lam": 1, "tsw": 0, "tsat": 0, "E": 0, "Q": 0},
    "N": {"lam": 1, "tsw": 0, "E": 0, "Q": 0},
    "T": {"lam": 1, "t90": 0, "tsat": 1, "E90": 0, "Q": 0},
    "U": {"lam": 1, "t90": 0, "E90": 0, "Q": 0},
}
cc = D.loc[m.SCORING_NODES]
tab, comp = {}, {}
for k, obj in VARIANTS.items():
    Sc, res = m.score_table(cc, obj)
    tab[k] = res
    comp[k] = Sc
    print(f"\nvariant {k}: " + ", ".join(f"{i}: {r.S:.3f} ({r['rank']})" for i, r in res.iterrows()),
          "| Pareto-optimal:", [int(i) for i in res.index[res.pareto]])
scores = pd.concat({k: v[["S", "rank", "pareto"]] for k, v in tab.items()}, axis=1)
scores.to_csv(m.RESULTS / "scores_coarse.csv")

print("\ncomponent scores (C) for 90 and 110 nm:\n", comp["C"].loc[[90, 110]].round(3).to_string())
dC = tab["C"].S[90] - tab["C"].S[110]
print("C: gap 90-110 = %.3f ; per-objective score differences:" % dC, (comp["C"].loc[90] - comp["C"].loc[110]).round(3).to_dict())
print("\nSpearman: tsw-tsat %.2f, tsw-Q %.2f, E-lam %.2f" % (spearmanr(cc.tsw, cc.tsat)[0], spearmanr(cc.tsw, cc.Q)[0],
                                                              spearmanr(cc.E, cc.lam)[0]))
ssum = comp["C"].tsw + comp["C"].tsat
print("score(tsw)+score(tsat) range: %.2f-%.2f" % (ssum.min(), ssum.max()))
cw = cc.assign(win=cc.tsat - cc.tsw)
_, rw = m.score_table(cw, {"lam": 1, "tsw": 0, "win": 1, "E": 0, "Q": 0})
print("window (tsat - tsw) instead of tsat:", rw.S.round(3).to_dict())
print("tsw/tsat max over scored nodes: %.1f %%" % (100 * (cc.tsw / cc.tsat).max()))

# ---------------------------------------------------------------- windows at 1 V (9 nodes)
rows = []
for tol in (0.10, 0.25):
    w = {"Large saturated Lambda": m.window(D.lam, True, tol),
         "Fast t_sw": m.window(D.tsw, False, tol),
         "Fast t_90": m.window(D.t90.iloc[:7], False, tol),         # 190/210: t90 not defined (no saturation)
         "Low E_sw": m.window(D.E, False, tol),
         "Low E_90": m.window(D.E90.iloc[:7], False, tol),
         "Low heating <Qh>": m.window(D.Q, False, tol),
         "Low power <P>": m.window(D.P, False, tol),
         "Low dT_max": m.window(D.dT, False, tol)}
    # t_sat: the maximum is censored (190 and 210 nm exceed 40 s); a window is only bracketed
    longest = D.tsat.dropna().max()
    w["Long t_sat (190,210 > 40 s)"] = [190, 210] if tol == 0.10 else [170, 190, 210]
    for k, v in w.items():
        rows.append(dict(tolerance=tol, purpose=k, window_nm=v))
    for nm in ("C", "D", "N", "T", "U"):
        r = tab[nm]
        rows.append(dict(tolerance=tol, purpose=f"balanced ({nm})",
                         window_nm=[int(i) for i in r.index[r.S >= r.S.max() * (1 - tol)]]))
W = pd.DataFrame(rows)
W.to_csv(m.RESULTS / "windows_coarse_1V.csv", index=False)
print("\nwindows at 1 V:")
print(W.to_string(index=False))
print("t_sat(170 nm) = %.2f s: inside the 25%% window only if t_sat(210 nm) <= %.1f s" % (longest, longest / 0.75))

# ---------------------------------------------------------------- windows common to all eleven voltages
res_all, qh_all = m.resistance(), m.ecqh()
summ = pd.read_csv(m.DATA / "resistance_all_thicknesses_summary.csv")
recs = []
for v in np.round(np.arange(0.2, 1.21, 0.1), 1):
    for th in m.COARSE_NODES:
        t, r = m.get_run(res_all, "Extra Fine", m.VS, th, v)
        sw = summ[(summ.sweep == m.VS) & (summ.thickness_nm == th) & (np.isclose(summ.voltage_V, v))]
        tsw = float(sw.switching_time_s.iloc[0])
        qs = "voltage sweep (50-90 nm)" if th <= 90 else "voltage sweep (110-210 nm, full-stack re-export)"
        tq, q = m.get_run(qh_all, "Extra Fine", qs, th, v, col="ecQh_W_per_m3")
        recs.append(dict(voltage_V=v, thickness_nm=th, tsw=tsw,
                         E=m.switching_energy(t, r, tsw, v), Q=q.mean()))
A = pd.DataFrame(recs)
A.to_csv(m.RESULTS / "voltage_sweep_metrics.csv", index=False)
out = []
for tol in (0.10, 0.25):
    for name, col, mx in [("t_sw", "tsw", False), ("E_sw", "E", False), ("heating", "Q", False)]:
        sets = [set(m.window(g.set_index("thickness_nm")[col], mx, tol)) for _, g in A.groupby("voltage_V")]
        common = sorted(set.intersection(*sets))
        out.append(dict(tolerance=tol, metric=name, common_to_all_voltages=common))
        print(f"tol {tol:.2f} {name:10s} common to all 11 voltages: {common}")
pd.DataFrame(out).to_csv(m.RESULTS / "windows_all_voltages.csv", index=False)
