"""Step 1. R_LRS, saturated R_HRS, Lambda = R_HRS/R_LRS, t_sat (99.5 % and 99 %), t50, t90, E90
and the automated switching time for every 1 V run (coarse sweep 50-210 nm; refined sweeps
80-110 nm on the Extra Fine and Extremely Fine meshes).

Output: results/saturation_metrics_1V.csv   (compare with data/*values_used_in_paper.csv)
"""
import numpy as np, pandas as pd
import memlib as m

rows = []
for th in m.COARSE_NODES:
    t, r = m.coarse_run(th)
    t995, ref = m.saturation(t, r, thr=0.995)
    t99, _ = m.saturation(t, r, thr=0.99)
    t50, t90 = m.completion_times(t, r)
    e90 = m.trapz(1 / r[t <= t90], t[t <= t90]) if np.isfinite(t90) else np.nan   # per V^2
    # automated t_sw is evaluated on the 15 s voltage-sweep run (as in the paper table)
    t15, r15 = m.get_run(m.resistance(), "Extra Fine", m.VS, th)
    rows.append(dict(set="coarse", mesh="Extra Fine", thickness_nm=th, run_length_s=t[-1], R_LRS_ohm=r[0],
                     R_ref_ohm=ref, lam=ref / r[0], reached=np.isfinite(t995), t_sat_995=t995, t_sat_99=t99,
                     t50=t50, t90=t90, E90_J=e90, tsw_auto=m.tsw_automated(t15, r15)))
for mesh in ("Extra Fine", "Extremely Fine"):
    for th in range(80, 111):
        t, r = m.refined_run(mesh, th)
        t995, ref = m.saturation(t, r, thr=0.995)
        t99, _ = m.saturation(t, r, thr=0.99)
        rows.append(dict(set="refined", mesh=mesh, thickness_nm=th, run_length_s=t[-1], R_LRS_ohm=r[0],
                         R_ref_ohm=ref, lam=ref / r[0], reached=np.isfinite(t995), t_sat_995=t995, t_sat_99=t99,
                         tsw_auto=m.tsw_automated(t, r)))
out = pd.DataFrame(rows)
out.to_csv(m.RESULTS / "saturation_metrics_1V.csv", index=False)

# ---- check against the values tabulated in the paper
C = m.load_coarse_table()
c = out[out.set == "coarse"].set_index("thickness_nm")
print("coarse: max |diff| vs paper table")
print("  Lambda  ", float(np.abs(c.lam - C.HRS_over_LRS_saturated).max()))
print("  R_LRS   ", float(np.abs(c.R_LRS_ohm - C.R_LRS_ohm).max()))
print("  t_sat   ", float(np.abs(c.t_sat_995 - C["t_sat_99.5pct_s"]).max()), "(nan = not reached)")
print("  t50,t90 ", float(np.abs(c.t50 - C.t50_s).max()), float(np.abs(c.t90 - C.t90_s).max()))
print("  tsw_auto", float(np.abs(c.tsw_auto - C.t_sw_automated_rule_s).max()))
R = m.load_refined_table()
for mesh in ("Extra Fine", "Extremely Fine"):
    r = out[(out.set == "refined") & (out.mesh == mesh)].set_index("thickness_nm")
    P = R[R.mesh == mesh].set_index("thickness_nm")
    print(mesh, "refined: Lambda", float(np.abs(r.lam - P.HRS_over_LRS_saturated).max()),
          " t_sat", float(np.abs(r.t_sat_995 - P["t_sat_99.5pct_s"]).max()),
          " reached", int(r.reached.sum()), "/31",
          " tsw_auto", float(np.abs(r.tsw_auto - P.t_sw_automated_rule_s).max()))
