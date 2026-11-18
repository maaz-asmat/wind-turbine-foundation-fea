"""Finite element model of the foundation slab: Reissner-Mindlin plate of variable thickness on a
tension-less Winkler soil (soil-structure interaction with gapping).

Unknowns: settlement w (positive downwards) and rotations (theta_x, theta_y), all quadratic (P2) triangles.
The soil only pushes, never pulls: springs are active where w > 0 (contact) and switched off where the
slab lifts off (gapping). The contact zone is found by iteration (active-set method).

Loads: self-weight of concrete and backfill, and the tower loads transferred by the anchor cage
("foundation basket") as a line load on the bolt circle: N / (2 pi a) + M cos(phi) / (pi a^2).
"""
from dataclasses import dataclass
import math
import time
import numpy as np
import triangle
from scipy.sparse import bmat
from skfem import MeshTri, Basis, ElementTriP2, BilinearForm, LinearForm, asm, solve
from skfem.helpers import grad, dot
from foundation import Geometry, Soil, Concrete, subgrade_modulus

KAPPA = 5.0 / 6.0


@dataclass
class PlateResult:
    mesh: object
    basis: object
    w: np.ndarray            # nodal settlement [m]
    tx: np.ndarray
    ty: np.ndarray
    k: float                 # subgrade modulus [kPa/m]
    contact_fraction: float
    p_max: float             # max soil pressure [kPa]
    iterations: int
    runtime_s: float
    n_dofs: int
    g: Geometry


def disc_mesh(g: Geometry, h_far=0.6, h_ring=0.2, ring_w=0.6):
    """Triangular mesh of the circular slab with circles at the anchor ring and the pedestal edge."""
    def circle(r, h, marker):
        n = max(24, int(2 * math.pi * r / h))
        t = np.linspace(0, 2 * math.pi, n, endpoint=False)
        pts = np.column_stack([r * np.cos(t), r * np.sin(t)])
        seg = np.column_stack([np.arange(n), (np.arange(n) + 1) % n])
        return pts, seg
    P, S, off = [], [], 0
    for r, h in [(g.R, h_far), (g.r_p, h_ring), (g.r_a - ring_w / 2, h_ring), (g.r_a + ring_w / 2, h_ring)]:
        p, s = circle(r, h, 0)
        P.append(p); S.append(s + off); off += len(p)
    t = triangle.triangulate({"vertices": np.vstack(P), "segments": np.vstack(S)}, f"pq30a{0.5 * h_far**2:.4f}")
    # refine inside the pedestal / anchor zone
    for _ in range(2):
        tri = t["triangles"]; v = t["vertices"]
        c = v[tri].mean(axis=1); rr = np.hypot(c[:, 0], c[:, 1])
        area = np.where(rr < g.r_p + 1.0, 0.5 * h_ring**2, 0.5 * h_far**2)
        t = triangle.triangulate({"vertices": v, "triangles": tri, "segments": t["segments"],
                                  "triangle_max_area": area}, "rpq30")
    return MeshTri(t["vertices"].T.copy(), t["triangles"].T.copy())


def solve_plate(g: Geometry, s: Soil, c: Concrete, N_turbine, M, h_far=0.6, ring_w=0.6,
                k=None, E=None, max_iter=50, mesh=None, q_uniform=None, ss_edge=False, const_h=None):
    """Solve the slab. Test options: q_uniform (only a uniform pressure, no self-weight/tower load),
    ss_edge (w = 0 on the outer edge, no springs) and const_h (constant thickness)."""
    if const_h is not None:
        from dataclasses import replace
        g = replace(g, h_c=const_h, h_e=const_h)
    t0 = time.time()
    mesh = mesh if mesh is not None else disc_mesh(g, h_far=h_far, ring_w=ring_w)
    b = Basis(mesh, ElementTriP2(), intorder=4)
    E = E if E is not None else 22e3 * ((c.fck + 8) / 10) ** 0.3 * 1e3     # Ecm (EC2) in kPa
    nu_c = 0.2
    k = k if k is not None else subgrade_modulus(g, s)

    def h_at(x):
        return g.thickness(np.hypot(x[0], x[1]))

    def Dm(x):
        return E * h_at(x) ** 3 / (12 * (1 - nu_c**2))

    def Gh(x):
        return KAPPA * E / (2 * (1 + nu_c)) * h_at(x)

    # stiffness blocks (u: w, tx, ty)
    @BilinearForm
    def kxx(u, v, w):   # bending, theta_x - theta_x
        return Dm(w.x) * (u.grad[0] * v.grad[0] + 0.5 * (1 - nu_c) * u.grad[1] * v.grad[1])

    @BilinearForm
    def kyy(u, v, w):
        return Dm(w.x) * (u.grad[1] * v.grad[1] + 0.5 * (1 - nu_c) * u.grad[0] * v.grad[0])

    @BilinearForm
    def kxy(u, v, w):   # u = theta_y (trial), v = theta_x (test)
        return Dm(w.x) * (nu_c * u.grad[1] * v.grad[0] + 0.5 * (1 - nu_c) * u.grad[0] * v.grad[1])

    @BilinearForm
    def s_ww(u, v, w):
        return Gh(w.x) * dot(grad(u), grad(v))

    @BilinearForm
    def s_wtx(u, v, w):  # trial theta_x, test w
        return -Gh(w.x) * u * v.grad[0]

    @BilinearForm
    def s_wty(u, v, w):
        return -Gh(w.x) * u * v.grad[1]

    @BilinearForm
    def s_tt(u, v, w):
        return Gh(w.x) * u * v

    @BilinearForm
    def spring(u, v, w):
        return k * w.chi * u * v

    Kxx, Kyy, Kxy = asm(kxx, b), asm(kyy, b), asm(kxy, b)
    Sww, Swx, Swy, Stt = asm(s_ww, b), asm(s_wtx, b), asm(s_wty, b), asm(s_tt, b)

    # loads [kPa]
    a = g.r_a

    @LinearForm
    def load(v, w):
        r = np.hypot(w.x[0], w.x[1])
        if q_uniform is not None:
            return q_uniform * v + 0 * r
        h = g.thickness(r)
        q = c.unit_weight * h + s.gamma * np.clip(g.D - h, 0, None)
        ring = (np.abs(r - a) <= ring_w / 2)
        cosphi = w.x[0] / np.maximum(r, 1e-9)
        q = q + ring * (N_turbine / (2 * math.pi * a) + M * cosphi / (math.pi * a**2)) / ring_w
        return q * v

    f_w = asm(load, b)
    n = b.N
    F = np.concatenate([f_w, np.zeros(n), np.zeros(n)])

    chi_nodes = np.ones(n)
    if ss_edge:                                            # simply supported test plate, no soil
        K = bmat([[Sww, Swx, Swy], [Swx.T, Kxx + Stt, Kxy], [Swy.T, Kxy.T, Kyy + Stt]], format="csr")
        from skfem import condense
        D_w = b.get_dofs().all()                          # boundary nodes of the w field
        U = solve(*condense(K, F, D=D_w))
        n_ = b.N
        return PlateResult(mesh, b, U[:n_], U[n_:2 * n_], U[2 * n_:], 0.0, 1.0, 0.0, 1, time.time() - t0, 3 * n_, g)
    for it in range(1, max_iter + 1):
        chi = b.interpolate(chi_nodes)
        Kspr = asm(spring, b, chi=(chi.value > 0.5).astype(float))
        K = bmat([[Sww + Kspr, Swx, Swy],
                  [Swx.T, Kxx + Stt, Kxy],
                  [Swy.T, Kxy.T, Kyy + Stt]], format="csr")
        U = solve(K, F)
        wv = U[:n]
        new = (wv > 0).astype(float)
        if np.array_equal(new, chi_nodes):
            break
        chi_nodes = new
    tx, ty = U[n:2 * n], U[2 * n:]

    # contact area and soil pressure (on a fine evaluation grid inside the disc)
    p_nodes = k * np.clip(wv, 0, None)
    area_c = asm(LinearForm(lambda v, w: (w.chi > 0) * v), b, chi=b.interpolate(chi_nodes).value).sum()
    contact = float(area_c / (math.pi * g.R**2))
    return PlateResult(mesh, b, wv, tx, ty, k, contact, float(p_nodes.max()), it, time.time() - t0, 3 * n, g)


def stress_resultants(res: PlateResult, c: Concrete, E=None):
    """Bending moments [kNm/m] and shear forces [kN/m] at the quadrature points (x, mx, my, mxy, qx, qy)."""
    b, g = res.basis, res.g
    E = E if E is not None else 22e3 * ((c.fck + 8) / 10) ** 0.3 * 1e3
    nu_c = 0.2
    x = b.global_coordinates().value
    h = g.thickness(np.hypot(x[0], x[1]))
    D = E * h**3 / (12 * (1 - nu_c**2))
    Gh = KAPPA * E / (2 * (1 + nu_c)) * h
    tx, ty, w = b.interpolate(res.tx), b.interpolate(res.ty), b.interpolate(res.w)
    mx = D * (tx.grad[0] + nu_c * ty.grad[1])
    my = D * (ty.grad[1] + nu_c * tx.grad[0])
    mxy = D * 0.5 * (1 - nu_c) * (tx.grad[1] + ty.grad[0])
    qx = Gh * (w.grad[0] - tx.value)
    qy = Gh * (w.grad[1] - ty.value)
    return x, mx, my, mxy, qx, qy


def radial_profiles(res: PlateResult, c: Concrete, nbins=60):
    """Envelope over the angle of the radial and tangential moments and of the radial shear, as functions
    of the radius. Negative moment = tension at the bottom face (soil pushing up)."""
    x, mx, my, mxy, qx, qy = stress_resultants(res, c)
    r = np.hypot(x[0], x[1]).ravel()
    cs, sn = (x[0] / np.maximum(np.hypot(x[0], x[1]), 1e-9)).ravel(), (x[1] / np.maximum(np.hypot(x[0], x[1]), 1e-9)).ravel()
    mx, my, mxy, qx, qy = (a.ravel() for a in (mx, my, mxy, qx, qy))
    mrr = mx * cs**2 + my * sn**2 + 2 * mxy * sn * cs
    mtt = mx * sn**2 + my * cs**2 - 2 * mxy * sn * cs
    qr = qx * cs + qy * sn
    edges = np.linspace(0, res.g.R, nbins + 1)
    idx = np.clip(np.digitize(r, edges) - 1, 0, nbins - 1)
    out = {"r": 0.5 * (edges[1:] + edges[:-1])}
    for name, v in [("mrr", mrr), ("mtt", mtt)]:
        out[name + "_bottom"] = np.array([-min(v[idx == i].min(), 0) if np.any(idx == i) else 0 for i in range(nbins)])
        out[name + "_top"] = np.array([max(v[idx == i].max(), 0) if np.any(idx == i) else 0 for i in range(nbins)])
    out["qr_abs"] = np.array([np.abs(qr[idx == i]).max() if np.any(idx == i) else 0 for i in range(nbins)])
    return out
