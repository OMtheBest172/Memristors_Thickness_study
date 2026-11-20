"""Step 4. Weight sensitivity of the coarse ranking (Sec. 6.2): random Dirichlet(1) weights,
leave-one-objective-out ranking and its margins.   Outputs: results/weight_sensitivity.csv"""
import numpy as np, pandas as pd
import memlib as m

SEED, NDRAW = 1, 20000
C = m.load_coarse_table()
S1 = pd.read_csv(m.RESULTS / "saturation_metrics_1V.csv")
S1 = S1[S1.set == "coarse"].set_index("thickness_nm")
D = pd.DataFrame({"lam": S1.lam, "tsw": C.t_sw_selected_s, "tsat": S1.t_sat_995, "E": C.E_sw_J,
                  "Q": C.mean_ecQh_full_stack_W_per_m3}).loc[m.SCORING_NODES]
VARS = {"C": {"lam": 1, "tsw": 0, "tsat": 1, "E": 0, "Q": 0},
        "D": {"lam": 1, "tsw": 0, "tsat": 0, "E": 0, "Q": 0},
        "N": {"lam": 1, "tsw": 0, "E": 0, "Q": 0}}
rng = np.random.default_rng(SEED)
rows = []
for nm, obj in VARS.items():
    Sc, _ = m.score_table(D, obj)
    W = rng.dirichlet(np.ones(len(obj)), NDRAW)
    win = np.array(Sc.index)[(Sc.values @ W.T).argmax(0)]
    share = pd.Series(win).value_counts(normalize=True)
    print(f"variant {nm}: first place over {NDRAW} Dirichlet(1) weight vectors:",
          {int(k): round(100 * v, 1) for k, v in share.items()})
    for k, v in share.items():
        rows.append(dict(variant=nm, test="dirichlet_first_place_pct", key=int(k), value=100 * v))
    if nm == "C":
        for drop in obj:
            o2 = {a: b for a, b in obj.items() if a != drop}
            _, r = m.score_table(D, o2)
            top = r.S.sort_values(ascending=False)
            print(f"  without {drop:5s}: ranking {list(top.index[:3])}  margin 1st-2nd = {top.iloc[0]-top.iloc[1]:.4f}")
            rows.append(dict(variant="C", test="leave_one_out_margin", key=drop, value=top.iloc[0] - top.iloc[1]))
pd.DataFrame(rows).to_csv(m.RESULTS / "weight_sensitivity.csv", index=False)
