"""Step 8. Regenerate the computed figures of the paper from the CSV files in ../data and the results of
steps 1-6.  Figures that are not produced here: the schematic and the COMSOL vacancy-profile images.
Output: figures/*.png   (use  --paper  to also copy them to ../fig, where the LaTeX source looks for them)
"""
import ast, shutil, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import memlib as m

plt.rcParams.update({"font.size": 8, "axes.titlesize": 7.5, "axes.labelsize": 8})
m.FIG.mkdir(exist_ok=True)
BLUE, ORANGE = "#1f77b4", "#ff7f0e"
cmap = plt.get_cmap("viridis")

S1 = pd.read_csv(m.RESULTS / "saturation_metrics_1V.csv")
C = m.load_coarse_table()
H = pd.read_csv(m.RESULTS / "energy_heating_coarse_1V.csv", index_col=0)
sc = S1[S1.set == "coarse"].set_index("thickness_nm")
x = sc.index.values
ok = sc.reached.values

# ------------------------------------------------------------------ coarse sweep at 1 V
c = pd.DataFrame({"lam": sc.lam, "tsw": C.t_sw_selected_s, "tsat": sc.t_sat_995, "E": C.E_sw_J,
                  "Q": C.mean_ecQh_full_stack_W_per_m3})
fig, ax = plt.subplots(1, 5, figsize=(8.6, 2.4))
def panel(a, y, title):
    a.plot(x, y, "-", color=BLUE, lw=1.2); a.plot(x[ok], np.asarray(y)[ok], "o", color=BLUE, ms=4)
    a.plot(x[~ok], np.asarray(y)[~ok], "o", mfc="white", color=BLUE, ms=4)
    a.set_title(title); a.set_xlabel("$d_{ox}$ (nm)"); a.grid(alpha=.3)
panel(ax[0], c.lam, "HRS/LRS (saturated)")
panel(ax[1], c.tsw, "$t_{sw}$ (s)"); ax[1].plot(x[~ok], c.tsw[~ok], "o", color=BLUE, ms=4)
ts = c.tsat.copy(); ts[~ok] = 40
ax[2].plot(x, ts, "-", color=BLUE, lw=1.2, alpha=.5); ax[2].plot(x[ok], ts[ok], "o", color=BLUE, ms=4)
ax[2].plot(x[~ok], ts[~ok], "^", mfc="white", color=BLUE, ms=5)
ax[2].set_title("$t_{sat}$ (s); open = $>$40 s"); ax[2].set_xlabel("$d_{ox}$ (nm)"); ax[2].grid(alpha=.3)
for a, y, t in [(ax[3], c.E, "$E_{sw}$ (J)"), (ax[4], c.Q, "mean ec.Qh (W/m$^3$)")]:
    a.plot(x, y, "o-", color=BLUE, ms=4, lw=1.2); a.set_title(t); a.set_xlabel("$d_{ox}$ (nm)"); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(m.FIG / "coarse1V.png", dpi=300); plt.close()

# ------------------------------------------------------------------ refined sweeps
R = m.load_refined_table()
rs = S1[S1.set == "refined"]
ef = R[R.mesh == "Extra Fine"].set_index("thickness_nm"); xf = R[R.mesh == "Extremely Fine"].set_index("thickness_nm")
sef = rs[rs.mesh == "Extra Fine"].set_index("thickness_nm"); sxf = rs[rs.mesh == "Extremely Fine"].set_index("thickness_nm")
xs = np.arange(80, 111)
fig, ax = plt.subplots(2, 3, figsize=(10.5, 5.4)); ax = ax.ravel()
def rpanel(a, getter, ylabel, hollow=False):
    for d, s, col, mk, nm in [(ef, sef, ORANGE, "s", "Extra Fine"), (xf, sxf, BLUE, "o", "Extremely Fine")]:
        y = getter(d, s)
        a.plot(xs, np.polyval(np.polyfit(xs, y, 2), xs), "-", color=col, lw=1)
        if hollow:
            r = s.reached.values; a.plot(xs[r], y[r], mk, color=col, ms=3.5, label=nm); a.plot(xs[~r], y[~r], mk, mfc="white", color=col, ms=3.5)
        else:
            a.plot(xs, y, mk, color=col, ms=3.5, label=nm)
    a.set_xlabel("$d_{ox}$ (nm)"); a.set_ylabel(ylabel); a.grid(alpha=.3)
rpanel(ax[0], lambda d, s: s.lam.values, "HRS/LRS (saturated)", hollow=True); ax[0].legend(fontsize=6, loc="lower left")
for s, col, mk in [(sef, ORANGE, "s"), (sxf, BLUE, "o")]:
    r = s.reached.values
    ax[1].plot(xs[r], s.t_sat_995.values[r], mk, color=col, ms=3.5)
    ax[1].plot(xs[~r], np.full((~r).sum(), s.run_length_s.iloc[0]), "^", mfc="white", color=col, ms=4)
ax[1].set_xlabel("$d_{ox}$ (nm)"); ax[1].set_ylabel("$t_{sat}$ (s)"); ax[1].grid(alpha=.3)
rpanel(ax[2], lambda d, s: d.switching_time_s.values, "$t_{sw}$ (s)")
rpanel(ax[3], lambda d, s: d.switching_energy_J.values, "$E_{sw}$ (J)")
rpanel(ax[4], lambda d, s: d.mean_ecQh_W_per_m3.values, "mean $Q_h$ (W/m$^3$)")
ax[5].axis("off")
ax[5].text(0, .7, "Open symbols: HRS not saturated within the\nrun (lower bound for HRS/LRS); open triangles in\nthe $t_{sat}$ panel: 99.5% criterion not met within\nthe run (plotted at the run length). Curves:\nquadratic fits to all 31 points (not drawn for\n$t_{sat}$).", fontsize=7, va="top")
plt.tight_layout(); plt.savefig(m.FIG / "refined_80_110.png", dpi=250); plt.close()

# ------------------------------------------------------------------ scores (variants C, D, N)
sco = pd.read_csv(m.RESULTS / "scores_coarse.csv", header=[0, 1], index_col=0)
fig, a = plt.subplots(figsize=(3.6, 2.6)); i = np.arange(len(sco)); w = .27
a.bar(i - w, sco[("C", "S")], w, color=BLUE, label="C: $t_{sat}$ maximized (primary)")
a.bar(i, sco[("D", "S")], w, color=ORANGE, label="D: $t_{sat}$ minimized")
a.bar(i + w, sco[("N", "S")], w, color="#7f7f7f", label="N: no $t_{sat}$")
a.set_xticks(i); a.set_xticklabels(sco.index); a.set_xlabel("$d_{ox}$ (nm)"); a.set_ylabel("Equal-weight score $S$")
a.legend(fontsize=5.5, loc="upper right"); a.grid(alpha=.3, axis="y"); a.set_ylim(0, 1.15)
plt.tight_layout(); plt.savefig(m.FIG / "scores.png", dpi=250); plt.close()

# ------------------------------------------------------------------ windows
Wd = pd.read_csv(m.RESULTS / "windows_coarse_1V.csv")
Wd["window_nm"] = Wd.window_nm.apply(ast.literal_eval)
def get(tol, name):
    return list(Wd[(np.isclose(Wd.tolerance, tol)) & (Wd.purpose == name)].window_nm.iloc[0])
ROWS = [("Large resistance window", "Large saturated Lambda", [190, 210]),
        ("Fast switching ($t_{sw}$)", "Fast t_sw", []), ("Fast completion ($t_{90}$)", "Fast t_90", []),
        ("Low energy to peak rate ($E_{sw}$)", "Low E_sw", []), ("Low energy to 90% ($E_{90}$)", "Low E_90", []),
        ("Low volumetric heating ($\\bar{Q}_h$)", "Low heating <Qh>", []), ("Low total power ($P$)", "Low power <P>", []),
        ("Low temperature rise ($\\Delta T$)", "Low dT_max", []),
        ("Long analog window ($t_{sat}$)", "Long t_sat (190,210 > 40 s)", [170, 190, 210]),
        ("Balanced (score within 10% / 25% of best)", "balanced (C)", [])]
fig, a = plt.subplots(figsize=(7.4, 4.0))
for k, (label, key, prov) in enumerate(ROWS):
    y = len(ROWS) - k
    for tol, col, dy, al in [(0.25, ORANGE, .12, .55), (0.10, BLUE, -.12, .8)]:
        wnd = get(tol, key)
        if not wnd: continue
        a.plot([min(wnd), max(wnd)], [y + dy] * 2, "-", color=col, lw=5, alpha=al, solid_capstyle="butt")
        for n in wnd: a.plot(n, y + dy, "o", color=col, ms=4.5, mfc="white" if n in prov else col)
a.set_yticks(range(len(ROWS), 0, -1)); a.set_yticklabels([r[0] for r in ROWS]); a.set_xticks(m.COARSE_NODES)
a.set_xlim(40, 222); a.grid(alpha=.3, axis="x")
a.set_xlabel("WO$_3$ thickness $d_{ox}$ (nm); markers are simulated nodes (1 V); open markers: provisional")
a.legend([Line2D([0], [0], color=BLUE, lw=5, alpha=.8), Line2D([0], [0], color=ORANGE, lw=5, alpha=.55)],
         ["within 10% of best", "within 25% of best"], fontsize=6.5, loc="lower right")
plt.tight_layout(); plt.savefig(m.FIG / "windows.png", dpi=250); plt.close()

# ------------------------------------------------------------------ voltage sweep figures (tsw, energy, heating)
A = pd.read_csv(m.RESULTS / "voltage_sweep_metrics.csv")
ths = m.COARSE_NODES; vs_ = sorted(A.voltage_V.unique())
def two_panel(col, ylabel, fname, log=False, sci=False):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
    for k, th in enumerate(ths):
        d = A[A.thickness_nm == th]
        ax[0].plot(d.voltage_V, d[col], "o-", ms=3.5, lw=1.2, color=cmap(k / (len(ths) - 1)), label=f"{th} nm")
    for k, v in enumerate((0.2, 0.6, 1.0, 1.2)):
        d = A[np.isclose(A.voltage_V, v)]
        ax[1].plot(d.thickness_nm, d[col], "o-", ms=3.5, lw=1.2, label=f"{v:.1f} V")
    ax[0].set_xlabel("Voltage (V)"); ax[1].set_xlabel("WO$_3$ thickness (nm)")
    for a in ax:
        a.set_ylabel(ylabel); a.grid(alpha=.3)
        if log: a.set_yscale("log")
    ax[0].legend(fontsize=5.5, ncol=2); ax[1].legend(fontsize=5.5)
    plt.tight_layout(); plt.savefig(m.FIG / fname, dpi=200); plt.close()
two_panel("tsw", "Switching time (s)", "tsw.png", log=True)
two_panel("E", "Switching energy (J)", "energy.png")
fig, a = plt.subplots(figsize=(3.6, 2.8))
for k, v in enumerate(vs_):
    d = A[np.isclose(A.voltage_V, v)]
    a.plot(d.thickness_nm, d.Q, "o-", ms=2.5, lw=1, color=cmap(k / (len(vs_) - 1)), label=f"{v:.1f} V")
a.set_xlabel("$d_{ox}$ (nm)"); a.set_ylabel("15-s mean $Q_h$, Ag+WO$_3$+FTO (W/m$^3$)"); a.grid(alpha=.3)
a.legend(fontsize=5, ncol=2); plt.tight_layout(); plt.savefig(m.FIG / "qh_three_domain.png", dpi=200); plt.close()

# ------------------------------------------------------------------ switching-time example (50 nm, 0.2 V)
t, r = m.get_run(m.resistance(), "Extra Fine", m.VS, 50, 0.2)
summ = pd.read_csv(m.DATA / "resistance_all_thicknesses_summary.csv")
tsw = float(summ[(summ.sweep == m.VS) & (summ.thickness_nm == 50) & np.isclose(summ.voltage_V, 0.2)].switching_time_s.iloc[0])
dr = np.gradient(r, t); E = m.switching_energy(t, r, tsw, 0.2)
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.5))
ax[0].plot(t, r, color="royalblue"); ax[0].axvline(tsw, color="r", ls="--", lw=1); ax[0].plot(tsw, r[np.argmin(abs(t - tsw))], "ro")
ax[1].plot(t, dr, color="darkgreen"); ax[1].axvline(tsw, color="r", ls="--", lw=1); ax[1].plot(tsw, dr[np.argmin(abs(t - tsw))], "ro", label="selected maximum of dR/dt")
ax[0].set_xlabel("Time (s)"); ax[0].set_ylabel("Resistance ($\\Omega$)"); ax[1].set_xlabel("Time (s)"); ax[1].set_ylabel("dR/dt ($\\Omega$/s)")
ax[1].legend(fontsize=6)
fig.suptitle(f"50 nm, 0.2 V: $t_{{sw}}$ = {tsw:.2f} s, $E_{{sw}}$ = {E:.3e} J", fontsize=8)
plt.tight_layout(); plt.savefig(m.FIG / "switching_example.png", dpi=200); plt.close()
print("figures written to", m.FIG)
if "--paper" in sys.argv:
    for f in m.FIG.glob("*.png"):
        shutil.copy(f, m.ROOT.parent / "fig" / f.name)
    print("copied to ../fig")
