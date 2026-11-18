"""Step 2 - Reference foundation: all design checks, soil pressure / gapping map, settlement and the
reinforcement demand along the radius."""
import sys, pathlib, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from foundation import Geometry, Loads, Soil, Concrete, Steel, ec2_bending_reinforcement
from design import design

OUT = ROOT / "results"
g, L, s, c = Geometry(), Loads(), Soil(), Concrete()
o, fe, prof = design(g)
with open(OUT / "02_reference_design.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["quantity", "value"])
    for k, v in o.items():
        w.writerow([k, round(v, 3) if isinstance(v, float) else v])
for k, v in o.items():
    print(k, round(v, 3) if isinstance(v, float) else v)

# ---------------------------------------------------------------- maps: soil pressure and settlement
m = fe.mesh
tri = mtri.Triangulation(m.p[0], m.p[1], m.t.T)
wv = fe.w[:m.p.shape[1]]
p = fe.k * np.clip(wv, 0, None)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
cs = ax[0].tricontourf(tri, p, levels=np.linspace(0, max(p.max(), 1), 21), cmap="viridis")
ax[0].tricontour(tri, wv, levels=[0], colors="r", linewidths=2)
fig.colorbar(cs, ax=ax[0], label="Soil pressure [kPa]")
cs = ax[1].tricontourf(tri, 1000 * wv, levels=21, cmap="coolwarm")
fig.colorbar(cs, ax=ax[1], label="Settlement [mm] (negative = uplift)")
for a, t in zip(ax[:2], [f"(a) Soil pressure, extreme load\n{100 * fe.contact_fraction:.0f} % in contact, red line = lift-off",
                         "(b) Settlement, extreme load"]):
    a.set_aspect("equal"); a.set_title(t); a.set_xlabel("x [m] (moment acts about y)"); a.set_ylabel("y [m]")
    for rr, ls in [(g.r_p, "-"), (g.r_a, ":")]:
        a.add_patch(plt.Circle((0, 0), rr, fill=False, color="k", ls=ls, lw=1))
# radial reinforcement demand
gF = 1.35
r = prof["r"]
As_b = [ec2_bending_reinforcement(gF * mm, h, c, Steel())[0] for mm, h in zip(prof["mrr_bottom"], g.thickness(r))]
As_t = [ec2_bending_reinforcement(gF * mm, h, c, Steel())[0] for mm, h in zip(prof["mrr_top"], g.thickness(r))]
ax[2].plot(r, gF * prof["mrr_bottom"], color="#1f5f8b", label="bottom tension, m_Ed [kNm/m]")
ax[2].plot(r, gF * prof["mrr_top"], color="#d97706", label="top tension, m_Ed [kNm/m]")
ax2 = ax[2].twinx()
ax2.plot(r, As_b, "--", color="#1f5f8b", label="A_s bottom [mm²/m]")
ax2.plot(r, As_t, "--", color="#d97706", label="A_s top [mm²/m]")
ax2.set_ylabel("Required radial reinforcement A_s [mm²/m] (dashed)")
ax[2].axvline(g.r_p, color="0.5", ls=":"); ax[2].text(g.r_p + 0.1, ax[2].get_ylim()[1] * 0.9, "pedestal face", fontsize=8)
ax[2].set(xlabel="Radius [m]", ylabel="Design radial moment (ULS) [kNm/m]", title="(c) Radial bending and EC2 reinforcement")
ax[2].legend(fontsize=7, loc="upper right"); ax[2].grid(alpha=.3)
fig.tight_layout(); fig.savefig(OUT / "02_reference_foundation.png", dpi=150)
