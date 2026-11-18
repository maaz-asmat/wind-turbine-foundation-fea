"""Step 3 - Parametric study and optimisation: the lightest foundation that passes every check.

Design variables: radius R, edge thickness h_e and thickness at the pedestal h_c (the embedment depth
follows the pedestal: D = h_c - 0.2 m). Every candidate is first screened with the fast analytical checks;
the survivors get the full FE soil-structure + EC2 design. Then the optimisation is repeated for four
soil classes to show which criterion governs the size.
"""
import sys, pathlib, csv, itertools, time
from dataclasses import replace
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from foundation import Geometry, Loads, Soil, Concrete, Criteria, analytical_checks
from design import design

OUT = ROOT / "results"
L, c, cr = Loads(), Concrete(), Criteria()
base = Geometry()


def geom(R, he, hc):
    return replace(base, R=R, h_e=he, h_c=hc, D=hc - 0.2)


def util(a, g, s):
    return {
        "no gapping (normal load)": a["e_normal_m"] / (g.R / 4),
        "contact >= 50 % (extreme)": cr.contact_extreme_min / max(a["contact_extreme"], 1e-9),
        "bearing pressure": a["q_effective_kPa"] / s.q_allow,
        "rotational stiffness": cr.K_rot_required / 1e6 / a["K_rot_GNm_rad"],
        "overturning": cr.fs_overturning / a["overturning_ratio"],
    }


# indicative unit prices (EUR), only used to compare designs with each other
PRICE_CONCRETE, PRICE_STEEL, PRICE_EXCAVATION = 160.0, 1600.0, 20.0     # per m3, per t, per m3 dug and backfilled


def cost(g, o):
    return (PRICE_CONCRETE * o["concrete_m3"] + PRICE_STEEL * o["steel_t"]
            + PRICE_EXCAVATION * np.pi * (g.R + 1.0) ** 2 * g.D)


grid_R = np.arange(8.0, 12.01, 0.5)
grid_he = np.round(np.arange(0.6, 1.61, 0.2), 2)
grid_hc = np.round(np.arange(2.4, 3.61, 0.3), 2)


def optimise(s: Soil, run_all_fe=False, label=""):
    cand = []
    for R, he, hc in itertools.product(grid_R, grid_he, grid_hc):
        g = geom(R, he, hc)
        a = analytical_checks(g, L, s, c, cr)
        u = util(a, g, s)
        cand.append((g, a, u, max(u.values()) <= 1.0))
    ok = sorted([x for x in cand if x[3]], key=lambda x: x[1]["concrete_m3"])
    results = []
    for g, a, u, _ in (ok if run_all_fe else ok):
        o, fe, prof = design(g, s=s)
        results.append((g, o, u))
        if not run_all_fe and o["all_ok"]:
            break
    return cand, results


# ---------------------------------------------------------------- 1. reference soil: full parametric study
t0 = time.time()
soil = Soil()
cand, res = optimise(soil, run_all_fe=True)
rows = []
for g, o, u in res:
    rows.append([g.R, g.h_e, g.h_c, round(o["concrete_m3"], 1), round(o["steel_t"], 1), round(o["N_total_kN"]),
                 round(o["contact_normal"], 3), round(o["fe_contact_extreme"], 3), round(o["fe_p_max_kPa"], 1),
                 round(o["K_rot_GNm_rad"], 1), round(o["m_Ed_face_kNm_per_m"]), o["shear_links_needed"], o["all_ok"],
                 max(u, key=u.get), round(max(u.values()), 3), round(cost(g, o) / 1e3)])
hdr = ["R_m", "h_edge_m", "h_pedestal_m", "concrete_m3", "steel_t", "N_total_kN", "contact_normal", "contact_extreme_FE",
       "p_max_FE_kPa", "K_rot_GNm_rad", "m_Ed_face_kNm_per_m", "shear_links", "all_checks_ok", "governing_check", "max_utilisation", "indicative_cost_kEUR"]
with open(OUT / "03_parametric_study.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
feasible = [r for r in rows if r[12]]
best = min(feasible, key=lambda r: r[3])
o_ref, _, _ = design(base, s=soil)
print(f"{len(cand)} candidates, {len(res)} passed the analytical screen, {len(feasible)} pass all checks ({time.time() - t0:.0f} s)")
print("reference:", round(o_ref["concrete_m3"], 1), "m3,", round(o_ref["steel_t"], 1), "t")
print("optimum  :", best)
best_cost = min(feasible, key=lambda r: r[15])
print("min cost :", best_cost)

with open(OUT / "03_optimum_vs_reference.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["design", "R_m", "h_edge_m", "h_pedestal_m", "concrete_m3", "steel_t", "governing_check"])
    w.writerow(["reference", base.R, base.h_e, base.h_c, round(o_ref["concrete_m3"], 1), round(o_ref["steel_t"], 1), "-"])
    w.writerow(["optimum", best[0], best[1], best[2], best[3], best[4], best[13]])
    w.writerow(["change_%", "", "", "", round(100 * (best[3] / o_ref["concrete_m3"] - 1), 1),
                round(100 * (best[4] / o_ref["steel_t"] - 1), 1), ""])
    w.writerow(["min_indicative_cost", best_cost[0], best_cost[1], best_cost[2], best_cost[3], best_cost[4], best_cost[13]])

# ---------------------------------------------------------------- 2. soil classes
soils = {"very soft clay": Soil(G=4e3, q_allow=120), "soft clay": Soil(G=8e3, q_allow=150), "firm clay / loose sand": Soil(G=20e3, q_allow=200),
         "medium-dense sand (reference)": Soil(G=50e3, q_allow=300), "dense sand / gravel": Soil(G=80e3, q_allow=400)}
srows = []
for name, s in soils.items():
    _, r_ = optimise(s)
    feas = [x for x in r_ if x[1]["all_ok"]]
    if not feas:
        srows.append([name, s.G / 1e3, s.q_allow, "", "", "", "no feasible design in grid", "", ""]); continue
    g, o, u = feas[0]
    srows.append([name, s.G / 1e3, s.q_allow, g.R, g.h_e, g.h_c, round(o["concrete_m3"], 1), round(o["steel_t"], 1),
                  max(u, key=u.get)])
    print(srows[-1])
with open(OUT / "03_soil_classes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["soil", "G_MPa", "q_allow_kPa", "R_m", "h_edge_m", "h_pedestal_m", "concrete_m3", "steel_t", "governing_check"])
    w.writerows(srows)

# ---------------------------------------------------------------- figures
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
F = np.array([[r[0], r[3], r[4]] for r in feasible])
NF = np.array([[r[0], r[3]] for r in rows if not r[12]]) if any(not r[12] for r in rows) else np.zeros((0, 2))
sc = ax[0].scatter(F[:, 0], F[:, 1], c=F[:, 2], cmap="viridis", s=28, label="passes all checks")
if len(NF):
    ax[0].scatter(NF[:, 0], NF[:, 1], marker="x", color="0.6", s=20, label="fails an FE/EC2 check")
ax[0].plot(base.R, o_ref["concrete_m3"], "s", color="crimson", ms=9, label="reference design")
ax[0].plot(best[0], best[3], "*", color="gold", mec="k", ms=16, label="optimum")
fig.colorbar(sc, ax=ax[0], label="Reinforcement [t]")
ax[0].set(xlabel="Foundation radius R [m]", ylabel="Concrete volume [m³]", title="(a) Designs passing the analytical screen")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
# (b) which criterion governs for the optimum geometry as the soil gets softer
g_opt = geom(best[0], best[1], best[2])
Gs = np.linspace(3e3, 100e3, 80)
U = {}
for G in Gs:
    s_ = Soil(G=G, q_allow=float(np.interp(G, [4e3, 8e3, 20e3, 50e3, 80e3], [120, 150, 200, 300, 400])))
    for k_, v_ in util(analytical_checks(g_opt, L, s_, c, cr), g_opt, s_).items():
        U.setdefault(k_, []).append(v_)
for (k_, v_), col in zip(U.items(), ["#1f5f8b", "#0e7c7b", "#d97706", "#b91c1c", "#7c3aed"]):
    ax[1].plot(Gs / 1e3, v_, color=col, lw=2, label=k_)
ax[1].axhline(1.0, color="k", ls="--", lw=1); ax[1].text(60, 1.03, "limit", fontsize=8)
ax[1].set(xscale="log", xlabel="Soil shear modulus G [MPa]", ylabel="Utilisation (demand / capacity)",
          title=f"(b) Optimum design (R = {best[0]:.1f} m): governing check vs soil", ylim=(0, 1.6))
ax[1].legend(fontsize=7, loc="upper right"); ax[1].grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig(OUT / "03_parametric_optimisation.png", dpi=150)
print(f"total {time.time() - t0:.0f} s")
