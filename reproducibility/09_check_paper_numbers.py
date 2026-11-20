"""Step 9. Compare the regenerated results with the numbers printed in the paper.
Every line is PASS/FAIL; tolerances are the rounding of the number as printed. Exit status 1 if any FAIL.
"""
import ast, sys
import numpy as np, pandas as pd
import memlib as m

R = m.RESULTS
fails = 0
def check(name, got, want, tol):
    global fails
    ok = abs(got - want) <= tol
    fails += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got:.4g}, paper {want:.4g} (tol {tol:g})")

sat = pd.read_csv(R / "saturation_metrics_1V.csv")
co = sat[sat.set == "coarse"].set_index("thickness_nm")
for th, lam, ts in [(50, 26.74, 2.53), (70, 29.64, 5.25), (90, 30.30, 10.48), (110, 29.98, 13.73), (130, 29.53, 18.74),
                    (150, 29.31, 28.30), (170, 29.10, 33.91)]:
    check(f"coarse {th} nm Lambda", co.lam[th], lam, 0.006)
    check(f"coarse {th} nm t_sat", co.t_sat_995[th], ts, 0.006)
check("190 nm Lambda lower bound", co.lam[190], 29.19, 0.006)
check("210 nm Lambda lower bound", co.lam[210], 28.52, 0.006)
check("R_LRS 50 nm (ohm)", co.R_LRS_ohm[50], 5861.57, 0.01)
check("R_LRS 210 nm / 50 nm", co.R_LRS_ohm[210] / co.R_LRS_ohm[50], 9.4, 0.05)

sc = pd.read_csv(R / "scores_coarse.csv", header=[0, 1], index_col=0)
paper = {"C": {50: .338, 70: .653, 90: .781, 110: .756, 130: .626, 150: .547, 170: .533},
         "D": {50: .538, 70: .818, 90: .880, 110: .813, 130: .619, 150: .419, 170: .333},
         "N": {50: .423, 70: .794, 90: .913, 110: .855, 130: .653, 150: .479, 170: .416},
         "T": {50: .400, 70: .636, 90: .695, 110: .636, 130: .595, 150: .578, 170: .533},
         "U": {50: .500, 70: .773, 90: .805, 110: .706, 130: .615, 150: .517, 170: .416}}
for v, d in paper.items():
    for th, s in d.items():
        check(f"score {v} {th} nm", sc[(v, "S")][th], s, 0.0006)

ws = pd.read_csv(R / "weight_sensitivity.csv")
g = lambda t, k, v="C": float(ws[(ws.variant == v) & (ws.test == t) & (ws.key.astype(str) == str(k))].value.iloc[0])
check("Dirichlet C: 90 nm first (%)", g("dirichlet_first_place_pct", 90), 71, 0.6)
check("Dirichlet C: 170 nm first (%)", g("dirichlet_first_place_pct", 170), 17, 0.6)
check("Dirichlet C: 110 nm first (%)", g("dirichlet_first_place_pct", 110), 12, 0.6)
check("Dirichlet D: 90 nm first (%)", g("dirichlet_first_place_pct", 90, "D"), 78, 0.6)
for k, v in [("tsat", .058), ("Q", .057), ("lam", .010), ("E", .002), ("tsw", .002)]:
    check(f"leave-one-out margin without {k}", g("leave_one_out_margin", k), v, 0.0006)

mc = pd.read_csv(R / "noise_monte_carlo.csv")
f = lambda v, k: float(mc[(mc.variant == v) & (mc.key == k)].value.iloc[0])
check("MC C: 90 nm first (%)", f("C", "first_90nm_pct"), 79, 0.6)
check("MC D: 90 nm first (%)", f("D", "first_90nm_pct"), 97, 0.6)
check("MC N: 90 nm first (%)", f("N", "first_90nm_pct"), 94, 0.6)
check("MC D: 90 outranks 110 (%)", f("D", "90_outranks_110_pct"), 99, 0.6)

W = pd.read_csv(R / "windows_coarse_1V.csv"); W["window_nm"] = W.window_nm.apply(ast.literal_eval)
def win(tol, name): return list(W[(np.isclose(W.tolerance, tol)) & (W.purpose == name)].window_nm.iloc[0])
expected = {(0.1, "Large saturated Lambda"): list(range(70, 211, 20)), (0.1, "Fast t_sw"): [50],
            (0.1, "Low E_sw"): [70, 90, 110], (0.1, "Low heating <Qh>"): [190, 210], (0.1, "Low power <P>"): list(range(110, 211, 20)),
            (0.1, "Low dT_max"): [190, 210], (0.25, "Low heating <Qh>"): [170, 190, 210], (0.25, "Low power <P>"): list(range(90, 211, 20)),
            (0.25, "Low dT_max"): [150, 170, 190, 210], (0.1, "balanced (C)"): [90, 110], (0.1, "balanced (D)"): [70, 90, 110],
            (0.25, "balanced (C)"): [70, 90, 110, 130], (0.25, "balanced (D)"): [70, 90, 110], (0.25, "balanced (N)"): [70, 90, 110]}
for (tol, name), want in expected.items():
    ok = win(tol, name) == want; fails += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}  window {name} @{tol:.0%}: {win(tol, name)}")

tr = pd.read_csv(R / "refined_trend_fits.csv").set_index(["mesh", "metric"])
for (mesh, k), res, chg in [(("EF", "tsw"), 5.2, 43), (("EF", "E"), 3.0, -0.2), (("EF", "lam"), 0.28, -0.5), (("EF", "Q"), 0.55, -31),
                            (("XF", "tsw"), 3.0, 86), (("XF", "E"), 1.9, 19), (("XF", "lam"), 0.17, -6.8), (("XF", "Q"), 0.25, -29)]:
    check(f"{mesh} {k} residual %", tr.loc[(mesh, k), "residual_pct"], res, 0.06)
    check(f"{mesh} {k} trend change %", tr.loc[(mesh, k), "trend_change_pct"], chg, 0.6)

txt = (R / "refined_mesh_summary.txt").read_text()
for s in ("t_sw      mean 0.57", "mean -38.6%", "<Qh>      mean -9.1%", "29/31 reach saturation; t_sat 6.92-13.12",
          "19/31 reach saturation; t_sat 8.39-10.98", "mean +9.6%", "Spearman(t_sw, d): EF 0.84, XF 0.99",
          "raw top5 91 (0.820), 84 (0.717), 90 (0.709), 104 (0.628), 107 (0.620)", "raw top5 86 (0.770)",
          "Pareto-optimal 15/31", "Pareto-optimal 26/31"):
    ok = s in txt; fails += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}  refined/mesh summary contains: {s}")
print(f"\n{fails} check(s) failed")
sys.exit(1 if fails else 0)
