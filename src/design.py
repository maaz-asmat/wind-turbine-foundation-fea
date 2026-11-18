"""Complete foundation design check: analytical stability checks + FE soil-structure interaction +
Eurocode 2 reinforcement, shear and anchor-cage (foundation basket) bearing checks + material quantities."""
from dataclasses import dataclass, replace
import math
import numpy as np
from foundation import (Geometry, Loads, Soil, Concrete, Steel, Criteria, analytical_checks,
                        ec2_bending_reinforcement, ec2_shear_resistance)
from plate_fem import solve_plate, radial_profiles


def design(g: Geometry, L: Loads = Loads(), s: Soil = Soil(), c: Concrete = Concrete(), st: Steel = Steel(),
           cr: Criteria = Criteria(), h_far=0.6, flange_w=0.6, allowance=0.15):
    out = analytical_checks(g, L, s, c, cr)

    # ---- FE: characteristic extreme loads (soil pressure, gapping) and ULS (all loads x gamma_F) for RC design
    M_x = L.M_extreme + L.H_extreme * g.D
    fe = solve_plate(g, s, c, L.N_turbine, M_x, h_far=h_far)
    out.update(fe_contact_extreme=fe.contact_fraction, fe_p_max_kPa=fe.p_max, fe_dofs=fe.n_dofs,
               fe_runtime_s=fe.runtime_s, fe_tilt_mrad=_tilt(fe))
    # ULS: the FE problem is linear apart from the contact; scaling all loads by gamma_F keeps the same
    # contact zone, so the ULS stress resultants are gamma_F times the characteristic ones.
    prof = radial_profiles(fe, c)
    gF = cr.gamma_F
    r = prof["r"]
    h_r = g.thickness(r)

    # ---- reinforcement along the radius (bottom and top faces, radial and tangential bars)
    As = {}
    for face in ("bottom", "top"):
        for comp in ("mrr", "mtt"):
            As[f"{comp}_{face}"] = np.array([ec2_bending_reinforcement(gF * m, h, c, st)[0]
                                             for m, h in zip(prof[f"{comp}_{face}"], h_r)])   # mm2/m
    dr = r[1] - r[0]
    ring_area = 2 * math.pi * r * dr
    steel_m3 = sum(float(np.sum(a * 1e-6 * ring_area)) for a in As.values())
    # ---- shear along the radius: where v_Ed > v_Rd,c, vertical links are needed (EC2 6.2.3, cot(theta) = 2.5)
    links_m3 = 0.0
    shear_zone = []
    for i, ri in enumerate(r):
        if ri < g.r_p:
            continue
        h = float(h_r[i]); d = h - 0.107
        vE = gF * prof["qr_abs"][i]
        vR = ec2_shear_resistance(h, As["mrr_bottom"][i], c)
        if vE > vR:
            asw = vE / (0.9 * d * st.fyd * 1e3 * 2.5)              # m2 of links per m2 of slab
            links_m3 += asw * h * ring_area[i]
            shear_zone.append(ri)
    steel_m3 += links_m3
    steel_t = steel_m3 * st.density * (1 + allowance)

    # design values at the pedestal face (critical section for bending)
    i_face = int(np.argmin(np.abs(r - (g.r_p + 0.1))))
    m_face = gF * prof["mrr_bottom"][i_face]
    As_face, As_min_face, d_face = ec2_bending_reinforcement(m_face, g.thickness(g.r_p), c, st)

    # ---- shear: one-way shear at a distance d outside the pedestal face (EC2 6.2.2, no shear links)
    r_v = g.r_p + d_face
    i_v = int(np.argmin(np.abs(r - r_v)))
    v_Ed = gF * prof["qr_abs"][i_v]
    v_Rd = ec2_shear_resistance(float(g.thickness(r_v)), As["mrr_bottom"][i_v], c)

    # ---- anchor cage (foundation basket): concrete pressure under the bottom anchor flange (EC2 6.7)
    a = g.r_a
    q_line = gF * (L.N_turbine / (2 * math.pi * a) + M_x / (math.pi * a**2))      # kN/m on the compression side
    sigma_flange = q_line / flange_w / 1000                                       # MPa
    f_Rdu = min(c.fcd * math.sqrt(3.0), 3.0 * c.fcd)                              # A_c1/A_c0 ~ 3 (spread below flange)
    t_line = gF * M_x / (math.pi * a**2) - L.N_turbine / (2 * math.pi * a)        # kN/m tension side (bolt pull)
    # hanger (suspension) reinforcement that carries the anchor-cage uplift from the bottom flange to the top
    # of the slab: total tension on the windward half T = 2 M / (pi a), bars ~2 x pedestal height long
    T_up = 2 * gF * M_x / (math.pi * a) - gF * L.N_turbine / 2
    hanger_m3 = max(T_up, 0) / (st.fyd * 1e3) * 2 * g.h_c
    steel_m3 += hanger_m3
    steel_t = steel_m3 * st.density * (1 + allowance)

    out.update(
        m_Ed_face_kNm_per_m=m_face, As_face_mm2_per_m=As_face, As_min_mm2_per_m=As_min_face,
        rho_face_pct=100 * As_face * 1e-6 / d_face,
        m_top_max_kNm_per_m=gF * float(prof["mrr_top"].max()),
        v_Ed_kN_per_m=v_Ed, v_Rd_c_kN_per_m=v_Rd, flange_stress_MPa=sigma_flange, f_Rdu_MPa=f_Rdu,
        bolt_tension_kN_per_m=t_line, steel_t=steel_t, steel_kg_per_m3=1000 * steel_t / out["concrete_m3"],
        shear_links_t=links_m3 * st.density * (1 + allowance),
        hanger_t=hanger_m3 * st.density * (1 + allowance),
        shear_link_zone_m=f"{min(shear_zone):.1f}-{max(shear_zone):.1f}" if shear_zone else "none",
    )
    out["ok_fe_contact_extreme"] = bool(fe.contact_fraction >= cr.contact_extreme_min)
    out["ok_fe_bearing"] = bool(fe.p_max <= 1.0 * s.q_allow * 1.3)        # edge peak may exceed the mean allowable by 30 %
    out["ok_bending"] = bool(out["rho_face_pct"] <= 2.0)
    out["shear_links_needed"] = bool(v_Ed > v_Rd)
    out["ok_anchor_bearing"] = bool(sigma_flange <= f_Rdu)
    out["all_ok"] = all(v for k, v in out.items() if k.startswith("ok_"))
    return out, fe, prof


def _tilt(fe):
    """Rotation of the slab about the y-axis from the FE settlements (least-squares plane), mrad."""
    p = fe.mesh.p
    A = np.column_stack([np.ones(p.shape[1]), p[0], p[1]])
    w = fe.w[:p.shape[1]]
    coef, *_ = np.linalg.lstsq(A, w, rcond=None)
    return 1000 * abs(coef[1])
