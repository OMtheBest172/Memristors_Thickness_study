"""Step 6. Refined sweeps 80-110 nm on the Extra Fine (EF) and Extremely Fine (XF) meshes:
mesh sensitivity (Sec. 5.4), saturation statistics, quadratic trends and residuals, ranking under
noise (Sec. 6.3).  Inputs: data/refined_sweep_80_110nm_values_used_in_paper.csv and the raw R(t).
Output: results/refined_mesh_summary.txt (everything printed) and results/refined_trend_fits.csv
"""
import io, sys, contextlib
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import memlib as m

buf = io.StringIO()
def say(*a):
    print(*a); print(*a, file=buf)

R = m.load_refined_table()
S1 = pd.read_csv(m.RESULTS / "saturation_metrics_1V.csv")
S1 = S1[S1.set == "refined"]
ef = R[R.mesh == "Extra Fine"].set_index("thickness_nm")
xf = R[R.mesh == "Extremely Fine"].set_index("thickness_nm")
sef = S1[S1.mesh == "Extra Fine"].set_index("thickness_nm")
sxf = S1[S1.mesh == "Extremely Fine"].set_index("thickness_nm")
x = np.arange(80, 111)

# ------------------------------------------------ mesh sensitivity (XF relative to EF)
def stat(name, v, fmt):
    say(f"  {name:9s} mean {fmt(v.mean())}  range {fmt(v.min())} .. {fmt(v.max())}")
say("== XF relative to EF over 31 thicknesses")
stat("t_sw", xf.switching_time_s / ef.switching_time_s, lambda a: f"{a:.2f}")
stat("E_sw", 100 * (xf.switching_energy_J / ef.switching_energy_J - 1), lambda a: f"{a:+.1f}%")
stat("<Qh>", 100 * (xf.mean_ecQh_W_per_m3 / ef.mean_ecQh_W_per_m3 - 1), lambda a: f"{a:+.1f}%")
stat("R_LRS", 100 * (xf.R_initial_ohm / ef.R_initial_ohm - 1), lambda a: f"{a:+.3f}%")
both = sef.reached & sxf.reached
say(f"  Lambda over the {int(both.sum())} thicknesses where both runs saturate:")
stat("Lambda", 100 * (sxf.lam[both] / sef.lam[both] - 1), lambda a: f"{a:+.1f}%")
say("  Spearman(t_sw, d): EF %.2f, XF %.2f" % (spearmanr(x, ef.switching_time_s)[0], spearmanr(x, xf.switching_time_s)[0]))

# ------------------------------------------------ saturation statistics
say("\n== saturation statistics (99.5 % threshold)")
for nm, s in (("EF", sef), ("XF", sxf)):
    ok = s.reached
    fit = np.polyfit(s.index[ok], s.t_sat_995[ok], 1)
    resid = s.t_sat_995[ok] - np.polyval(fit, s.index[ok])
    say(f"  {nm}: {int(ok.sum())}/31 reach saturation; t_sat {s.t_sat_995.min():.2f}-{s.t_sat_995.max():.2f} s; "
        f"slope {fit[0]:.2f} s/nm; scatter about linear fit {resid.std():.2f} s; not reached: {list(s.index[~ok])}")
    ok99 = s.t_sat_99.notna()
    say(f"      99 % threshold: {int(ok99.sum())}/31 reach; t_sat {s.t_sat_99.min():.2f}-{s.t_sat_99.max():.2f} s")
# relative slope of R over the last 3 s of the XF runs
sl = {}
for th in x:
    t, r = m.refined_run("Extremely Fine", th)
    mk = t >= t[-1] - 3
    sl[th] = 100 * np.polyfit(t[mk], r[mk], 1)[0] / r[mk].mean()   # relative to the mean over the window
sl = pd.Series(sl)
say("  XF relative slope over the last 3 s (%%/s) for d >= 93 nm: %.1f .. %.1f" % (sl[sl.index >= 93].min(), sl[sl.index >= 93].max()))
# smoothing-window test of the t_sat scatter (EF)
for win in (11, 21, 51, 101):
    ts = {}
    for th in x:
        t, r = m.refined_run("Extra Fine", th)
        ts[th] = m.saturation(t, r, win=win)[0]
    ts = pd.Series(ts).dropna()
    fit = np.polyfit(ts.index, ts.values, 1)
    say(f"  EF t_sat scatter about linear fit with a {win/100:.2f} s median window: {(ts - np.polyval(fit, ts.index)).std():.2f} s (n={len(ts)})")
# 12 s vs 15 s heating mean, EF 96-110 nm
q = []
for th in range(96, 111):
    tq, qq = m.get_run(m.ecqh(), "Extra Fine", "refined sweep 96-110 nm", th, col="ecQh_W_per_m3")
    q.append(100 * (qq[tq <= 12 + 1e-9].mean() / qq.mean() - 1))
say("  EF 96-110 nm: 12 s mean heating exceeds the 15 s mean by %.1f-%.1f %%" % (min(q), max(q)))

# ------------------------------------------------ quadratic trends and residuals
def quad(y):
    p = np.polyfit(x, y, 2); f = np.polyval(p, x)
    return f, (y - f).std(ddof=3)

OBJ = {"lam": 1, "tsw": 0, "E": 0, "Q": 0}
cols = {"lam": "lam", "tsw": "switching_time_s", "E": "switching_energy_J", "Q": "mean_ecQh_W_per_m3"}
fits, out = {}, []
say("\n== quadratic trends over 31 thicknesses (residual std with 3 fitted parameters)")
for nm, d, s in (("EF", ef, sef), ("XF", xf, sxf)):
    df = pd.DataFrame({k: (s.lam.values if k == "lam" else d[c].values) for k, c in cols.items()}, index=x)
    tr, sg = {}, {}
    for k in OBJ:
        tr[k], sg[k] = quad(df[k].values)
    fits[nm] = (df, pd.DataFrame(tr, index=x), sg)
    say(f"  {nm}: residual % " + ", ".join(f"{k} {100*sg[k]/np.mean(tr[k]):.2f}" for k in OBJ))
    say(f"      trend change over 80-110 nm: " + ", ".join(
        f"{k} {100*(tr[k][-1]/tr[k][0]-1):+.1f}% ({abs(tr[k][-1]-tr[k][0])/sg[k]:.1f} sigma)" for k in OBJ))
    for k in OBJ:
        out.append(dict(mesh=nm, metric=k, residual_pct=100 * sg[k] / np.mean(tr[k]),
                        trend_change_pct=100 * (tr[k][-1] / tr[k][0] - 1), sigma_multiple=abs(tr[k][-1] - tr[k][0]) / sg[k]))
pd.DataFrame(out).to_csv(m.RESULTS / "refined_trend_fits.csv", index=False)

# ------------------------------------------------ ranking (four objectives) and noise draws
NDRAW, rng = 4000, np.random.default_rng(2)
say("\n== ranking over 80-110 nm, objectives Lambda max; t_sw, E_sw, <Qh> min, equal weights")
for nm in ("EF", "XF"):
    df, T, sg = fits[nm]
    _, r = m.score_table(df, OBJ)
    top = r.S.sort_values(ascending=False).head(5)
    say(f"  {nm}: raw top5 " + ", ".join(f"{i} ({v:.3f})" for i, v in top.items()) + f"; Pareto-optimal {int(r.pareto.sum())}/31")
    _, rt = m.score_table(T, OBJ)
    say(f"      trend scores: 80 nm {rt.S[80]:.3f}, 100 nm {rt.S[100]:.3f}, 110 nm {rt.S[110]:.3f}; minimum at {int(rt.S.idxmin())} nm")
    cnt = pd.Series(0.0, index=x)
    for _ in range(NDRAW):
        N = pd.DataFrame({k: T[k].values + rng.normal(0, sg[k], len(x)) for k in OBJ}, index=x)
        _, rn = m.score_table(N, OBJ)
        cnt[rn.S.idxmax()] += 1
    say("      first place under noise (%d draws): " % NDRAW + ", ".join(f"{int(i)} nm {100*v/NDRAW:.1f}%" for i, v in cnt.sort_values(ascending=False).head(4).items()))
    if nm == "EF":
        ts = sef.t_sat_995.fillna(sef.run_length_s)
        _, r5 = m.score_table(df.assign(tsat=ts.values), {**OBJ, "tsat": 1})
        say("      EF with t_sat added (15 s for the two unsaturated runs): top5 " + ", ".join(
            f"{i} ({v:.3f})" for i, v in r5.S.sort_values(ascending=False).head(5).items()))
(m.RESULTS / "refined_mesh_summary.txt").write_text(buf.getvalue())
