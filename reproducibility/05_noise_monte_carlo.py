"""Step 5. Propagation of run-to-run scatter through the coarse ranking (Sec. 6.2, 'Noise propagation').

Each of the seven scored nodes (50-170 nm) is perturbed with independent Gaussian errors,
Lambda +-0.5 %, t_sw +-0.015 s, t_sat +-1 s, E_sw +-3 %, <Qh> +-0.5 % (1 sigma; relative for Lambda,
E_sw and <Qh>, absolute for the times), the min-max scores are recomputed for every draw and the first
place is tallied. The errors are run-to-run scatter only (taken from the repeated 90 and 110 nm runs and
from the residuals of the refined sweeps); they do NOT include the mesh dependence (script 06).
Output: results/noise_monte_carlo.csv
"""
import numpy as np, pandas as pd
import memlib as m

SEED, NDRAW = 7, 20000
SIG = {"lam": ("rel", 0.005), "tsw": ("abs", 0.015), "tsat": ("abs", 1.0), "E": ("rel", 0.03), "Q": ("rel", 0.005)}
C = m.load_coarse_table()
S1 = pd.read_csv(m.RESULTS / "saturation_metrics_1V.csv")
S1 = S1[S1.set == "coarse"].set_index("thickness_nm")
D = pd.DataFrame({"lam": S1.lam, "tsw": C.t_sw_selected_s, "tsat": S1.t_sat_995, "E": C.E_sw_J,
                  "Q": C.mean_ecQh_full_stack_W_per_m3}).loc[m.SCORING_NODES]
VARS = {"C": {"lam": 1, "tsw": 0, "tsat": 1, "E": 0, "Q": 0},
        "D": {"lam": 1, "tsw": 0, "tsat": 0, "E": 0, "Q": 0},
        "N": {"lam": 1, "tsw": 0, "E": 0, "Q": 0}}
rng = np.random.default_rng(SEED)
nodes = np.array(m.SCORING_NODES)
noisy = {}
for k, (kind, s) in SIG.items():
    base = D[k].values[None, :]
    z = rng.normal(0, 1, (NDRAW, len(nodes)))
    noisy[k] = base * (1 + s * z) if kind == "rel" else base + s * z
rows = []
for nm, obj in VARS.items():
    tot = np.zeros((NDRAW, len(nodes)))
    for k, maxi in obj.items():
        x = noisy[k]
        lo, hi = x.min(1, keepdims=True), x.max(1, keepdims=True)
        tot += ((x - lo) if maxi else (hi - x)) / (hi - lo)
    tot /= len(obj)
    first = nodes[tot.argmax(1)]
    i90, i110 = list(nodes).index(90), list(nodes).index(110)
    p = pd.Series(first).value_counts(normalize=True)
    beat = float((tot[:, i90] > tot[:, i110]).mean())
    print(f"variant {nm}: first place {({int(k): round(100*v, 1) for k, v in p.items()})}; 90 nm outranks 110 nm in {100*beat:.1f} %")
    for k, v in p.items():
        rows.append(dict(variant=nm, key=f"first_{int(k)}nm_pct", value=100 * v))
    rows.append(dict(variant=nm, key="90_outranks_110_pct", value=100 * beat))
pd.DataFrame(rows).to_csv(m.RESULTS / "noise_monte_carlo.csv", index=False)
