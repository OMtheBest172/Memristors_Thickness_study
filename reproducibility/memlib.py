"""Shared functions for the Ag/WO3/FTO thickness study (post-processing of COMSOL exports).

All scripts in this folder import this module. Inputs are the CSV files in ../data;
nothing here needs COMSOL. Tested with Python 3.10+, numpy>=1.24 (np.trapezoid is
used if available, otherwise np.trapz), pandas, scipy, matplotlib.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
def _find_data():
    """../data (this package layout), else the package root, else $MEM_DATA."""
    import os
    for cand in (os.environ.get("MEM_DATA"), ROOT.parent / "data", ROOT.parent, ROOT / "data"):
        if cand and (Path(cand) / "resistance_all_thicknesses_timeseries.csv").exists():
            return Path(cand)
    raise SystemExit("Cannot find resistance_all_thicknesses_timeseries.csv; set MEM_DATA to the data folder.")
DATA = _find_data()
RESULTS = ROOT / "results"
FIG = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)

trapz = getattr(np, "trapezoid", None) or np.trapz

# ---------------------------------------------------------------- constants (paper)
STACK_EXTRA_NM = 40.0          # two 20 nm electrodes
AREA_NM = 500.0                # lateral width and out-of-plane depth, nm
K_WO3 = 3.0                    # W/m/K, used only for the slab temperature estimate
COARSE_NODES = list(range(50, 211, 20))
SCORING_NODES = list(range(50, 171, 20))   # 190 and 210 nm did not saturate
VS = "voltage sweep (50-210 nm, 0.2-1.2 V, 15 s)"

_res = None
_qh = None


def resistance():
    """Resistance time series of every run (cached)."""
    global _res
    if _res is None:
        _res = pd.read_csv(DATA / "resistance_all_thicknesses_timeseries.csv", low_memory=False)
    return _res


def ecqh():
    """ec.Qh (surface average over Ag+WO3+FTO) time series of every run (cached)."""
    global _qh
    if _qh is None:
        _qh = pd.read_csv(DATA / "ecQh_all_thicknesses_timeseries.csv", low_memory=False)
    return _qh


def get_run(df, mesh, sweep, th, v=1.0, col="resistance_ohm"):
    d = df[(df.mesh == mesh) & (df.sweep == sweep) & (df.thickness_nm == th) & (df.voltage_V == v)]
    d = d.sort_values("time_s")
    return d.time_s.values, d[col].values


def coarse_sweep_name(th):
    """Which export gives the saturated metrics of the coarse nodes at 1 V (see paper Sec. 3.3/4.1)."""
    if th <= 90:
        return VS
    return "long run (20 s)" if th <= 130 else "long run (40 s)"


def coarse_run(th):
    return get_run(resistance(), "Extra Fine", coarse_sweep_name(th), th)


def refined_sweep_name(mesh, th):
    if th < 86:
        return "refined sweep 80-85 nm (12 s)"
    tail = "(15 s)" if mesh == "Extra Fine" else "(12 s)"
    return ("refined sweep 86-95 nm " if th < 96 else "refined sweep 96-110 nm ") + tail


def refined_run(mesh, th):
    return get_run(resistance(), mesh, refined_sweep_name(mesh, th), th)


# ---------------------------------------------------------------- resistance metrics
def smooth(r, win=11):
    """Centred rolling median (11 samples = 0.1 s at the 0.01 s output interval)."""
    return pd.Series(r).rolling(win, center=True, min_periods=1).median().values


def saturation(t, r, win=11, thr=0.995, hold=1.0, final=1.0):
    """Saturated HRS and saturation time.

    R_ref  = median of the smoothed R over the final `final` seconds of the run.
    t_sat  = earliest time at which smoothed R >= thr*R_ref and stays so for `hold` s more.
    Returns (t_sat or nan, R_ref). nan => criterion not met within the run (lower bound).
    """
    s = smooth(r, win)
    ref = np.median(s[t >= t[-1] - final + 1e-9])
    ok = s >= thr * ref
    dt = t[1] - t[0]
    n = int(round(hold / dt))
    bad = np.cumsum(~ok)
    for i in range(len(t) - n):
        if ok[i] and bad[i + n] - bad[i] == 0:
            return float(t[i]), float(ref)
    return float("nan"), float(ref)


def completion_times(t, r, fractions=(0.5, 0.9)):
    """Times at which the smoothed R has covered a given fraction of the rise R_LRS -> R_ref."""
    s = smooth(r)
    ref = np.median(s[t >= t[-1] - 1 + 1e-9])
    r0 = r[0]
    out = []
    for f in fractions:
        i = np.where(s >= r0 + f * (ref - r0))[0]
        out.append(float(t[i[0]]) if len(i) else float("nan"))
    return out


def switching_energy(t, r, t_end, v=1.0):
    """E = integral_0^t_end V^2/R dt (trapezoid on the exported samples with t <= t_end)."""
    m = t <= t_end + 1e-9
    return float(trapz(v ** 2 / r[m], t[m]))


def tsw_automated(t, r):
    """Automated switching-time rule: median of nine variants of argmax dR/dt restricted to
    the part of the run before R has completed 10/20/30 % of its rise, with no smoothing or a
    Savitzky-Golay filter of 5 or 9 samples (order 2). Used only as a cross-check of the
    manually selected t_sw."""
    from scipy.signal import savgol_filter
    s = smooth(r)
    ref = np.median(s[t >= t[-1] - 1])
    r0 = r[0]
    vals = []
    for sm in (None, 5, 9):
        y = r if sm is None else savgol_filter(r, sm, 2)
        dy = np.gradient(y, t)
        for f in (0.1, 0.2, 0.3):
            idx = np.where((s - r0) / (ref - r0) >= f)[0]
            end = idx[0] if len(idx) else len(t) - 1
            end = max(end, 2)
            vals.append(t[np.argmax(dy[:end + 1])])
    return float(np.median(vals))


# ---------------------------------------------------------------- multi-objective scoring
def minmax_scores(s, maximize):
    lo, hi = s.min(), s.max()
    return (s - lo) / (hi - lo) if maximize else (hi - s) / (hi - lo)


def score_table(df, objectives, weights=None):
    """Benefit scores in [0,1] (larger better), weighted sum and Pareto flag.

    objectives: dict {column: True if maximized else False}. Extrema are taken over the rows
    of `df`, i.e. the candidate set (scores of different candidate sets are not comparable).
    """
    S = pd.DataFrame({k: minmax_scores(df[k], m) for k, m in objectives.items()})
    w = np.ones(len(objectives)) / len(objectives) if weights is None else np.asarray(weights, float) / np.sum(weights)
    total = S.values @ w
    X = np.column_stack([df[k].values * (1 if m else -1) for k, m in objectives.items()])
    pareto = [not any(np.all(X[j] >= X[i]) and np.any(X[j] > X[i]) for j in range(len(df)) if j != i)
              for i in range(len(df))]
    out = pd.DataFrame({"S": total, "pareto": pareto}, index=df.index)
    out["rank"] = out.S.rank(ascending=False, method="min").astype(int)
    return S, out


def window(values, maximize, tol):
    """Thicknesses whose value lies within `tol` of the best value (best = max or min of `values`)."""
    best = values.max() if maximize else values.min()
    ok = values >= best * (1 - tol) if maximize else values <= best * (1 + tol)
    return [int(i) for i in values.index[ok]]


def load_coarse_table():
    """Paper table values of the coarse sweep at 1 V (manually selected t_sw, table E_sw)."""
    return pd.read_csv(DATA / "coarse_sweep_1V_values_used_in_paper.csv").set_index("thickness_nm")


def load_refined_table():
    return pd.read_csv(DATA / "refined_sweep_80_110nm_values_used_in_paper.csv")
