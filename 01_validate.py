"""Step 1 - Verification of the plate/soil FE model.

(a) Simply supported circular plate under uniform pressure vs the analytical Mindlin solution.
(b) A very stiff slab on tension-less springs vs the rigid-footing solution (contact area, edge pressure).
(c) Mesh convergence of the reference foundation.
"""
import sys, pathlib, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from foundation import Geometry, Loads, Soil, Concrete, rigid_contact, subgrade_modulus, vertical_load
from plate_fem import solve_plate, radial_profiles
from design import design

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
g, L, s, c = Geometry(), Loads(), Soil(), Concrete()
E = 22e3 * ((c.fck + 8) / 10) ** 0.3 * 1e3

# ------------------------------------------------------------ (a) simply supported plate
rows = []
q = 100.0
for h in (1.0, 3.0):
    D = E * h**3 / (12 * (1 - 0.2**2)); G = E / (2 * 1.2)
    w_exact = q * g.R**4 * (5 + 0.2) / (64 * D * (1 + 0.2)) + q * g.R**2 / (4 * 5 / 6 * G * h)
    for hf in (2.0, 1.0, 0.6):
        res = solve_plate(g, s, c, 0, 0, h_far=hf, q_uniform=q, ss_edge=True, const_h=h)
        w_fe = res.w.max()
        rows.append([h, hf, res.n_dofs, round(w_fe * 1000, 4), round(w_exact * 1000, 4), round(100 * (w_fe / w_exact - 1), 2)])
        print(rows[-1])
with open(OUT / "01a_plate_verification.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["thickness_m", "element_size_m", "dofs", "w_centre_FE_mm", "w_centre_Mindlin_mm", "error_%"]); w.writerows(rows)

# ------------------------------------------------------------ (b) rigid limit
N = vertical_load(g, L, s, c)
rows = []
for M in np.linspace(0, 0.42 * g.R * N, 8):
    # the self-weight acts as a uniform pressure in the FE; the tower moment is applied on the anchor ring
    pmax_r, cont_r, e = rigid_contact(N, M, g.R)
    res = solve_plate(g, s, c, L.N_turbine, M, E=E * 1000)
    rows.append([round(e / g.R, 3), round(cont_r, 3), round(res.contact_fraction, 3), round(pmax_r, 1), round(res.p_max, 1)])
    print(rows[-1])
with open(OUT / "01b_rigid_limit.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["e_over_R", "contact_rigid_theory", "contact_FE", "p_max_rigid_theory_kPa", "p_max_FE_kPa"]); w.writerows(rows)
rig = np.array(rows, float)

# ------------------------------------------------------------ (c) mesh convergence of the reference design
mrows = []
for hf in (1.2, 0.8, 0.6, 0.4):
    o, fe, prof = design(g, h_far=hf)
    mrows.append([hf, fe.n_dofs, round(fe.p_max, 1), round(fe.contact_fraction, 3), round(o["m_Ed_face_kNm_per_m"]),
                  round(o["v_Ed_kN_per_m"]), round(fe.runtime_s, 1)])
    print(mrows[-1])
with open(OUT / "01c_mesh_convergence.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["element_size_m", "dofs", "p_max_kPa", "contact_fraction", "m_Ed_pedestal_face_kNm_per_m", "v_Ed_kN_per_m", "runtime_s"])
    w.writerows(mrows)

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
ax[0].plot(rig[:, 0], 100 * rig[:, 1], "k-", label="rigid footing theory")
ax[0].plot(rig[:, 0], 100 * rig[:, 2], "o", color="#1f5f8b", label="FE, very stiff slab")
ax[0].axvline(0.25, color="0.5", ls=":"); ax[0].text(0.255, 85, "kern\ne = R/4", fontsize=8, color="0.4")
ax[0].set(xlabel="Eccentricity e / R", ylabel="Base in contact [%]", title="(a) Contact area (gapping)")
ax[1].plot(rig[:, 0], rig[:, 3], "k-", label="rigid footing theory")
ax[1].plot(rig[:, 0], rig[:, 4], "o", color="#1f5f8b", label="FE, very stiff slab")
ax[1].set(xlabel="Eccentricity e / R", ylabel="Max. soil pressure [kPa]", title="(b) Edge pressure")
for a in ax:
    a.grid(alpha=.3); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "01_rigid_limit_verification.png", dpi=150)
