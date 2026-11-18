"""Onshore wind turbine gravity foundation: geometry, loads, soil and the classical design checks.

Units: m, kN, kPa (kN/m2), MPa where stated. Moments in kNm.

Geometry (circular footing, conical top):

        |<- r_p ->|            pedestal (anchor cage / "foundation basket" inside)
        +---------+
       /           \\           conical top, backfilled with soil
      /             \\
  +--+               +--+  h_e (edge thickness)
  |                     |
  +---------------------+  h_c = thickness at the pedestal (incl. pedestal)
  |<-------- R -------->|
"""
from dataclasses import dataclass, field, replace
import math
import numpy as np


@dataclass
class Geometry:
    R: float = 10.0            # footing radius [m]
    h_e: float = 1.0           # edge thickness [m]
    h_c: float = 3.2           # thickness at the pedestal, slab + pedestal [m]
    r_p: float = 3.0           # pedestal radius [m]
    r_a: float = 2.2           # anchor cage (bolt circle) radius [m]
    D: float = 3.0             # embedment depth of the underside below ground [m]

    def thickness(self, r):
        r = np.asarray(r, float)
        return np.where(r <= self.r_p, self.h_c,
                        self.h_c + (self.h_e - self.h_c) * (r - self.r_p) / (self.R - self.r_p))

    def concrete_volume(self, n=2000):
        r = np.linspace(0, self.R, n)
        return float(np.trapezoid(2 * np.pi * r * self.thickness(r), r))

    def backfill_volume(self, n=2000):
        r = np.linspace(0, self.R, n)
        return float(np.trapezoid(2 * np.pi * r * np.clip(self.D - self.thickness(r), 0, None), r))


@dataclass
class Loads:
    """Characteristic tower-base loads (assumed, representative of a 4-5 MW onshore turbine, ~120 m hub)."""
    N_turbine: float = 5000.0      # vertical load from tower + nacelle + rotor [kN]
    M_extreme: float = 90000.0     # resultant overturning moment, extreme load case [kNm]
    H_extreme: float = 1200.0      # horizontal shear, extreme [kN]
    M_normal: float = 55000.0      # overturning moment, normal production (frequent) [kNm]
    H_normal: float = 800.0


@dataclass
class Soil:
    G: float = 50e3            # shear modulus [kPa] (medium-dense sand; used for stiffness and subgrade)
    nu: float = 0.30
    gamma: float = 18.0        # backfill unit weight [kN/m3]
    q_allow: float = 300.0     # allowable (characteristic) bearing pressure on the effective area [kPa]
    mu: float = 0.5            # base friction coefficient (tan of interface friction angle)
    H_layer: float = 30.0      # depth to a stiff layer (for the stiffness formula) [m]


@dataclass
class Concrete:
    fck: float = 35.0          # C35/45 [MPa]
    gamma_c: float = 1.5
    unit_weight: float = 25.0  # reinforced concrete [kN/m3]

    @property
    def fctm(self):            # EC2 Table 3.1
        return 0.30 * self.fck ** (2 / 3)

    @property
    def fcd(self):
        return 0.85 * self.fck / self.gamma_c


@dataclass
class Steel:
    fyk: float = 500.0         # B500B [MPa]
    gamma_s: float = 1.15
    density: float = 7.85      # t/m3

    @property
    def fyd(self):
        return self.fyk / self.gamma_s


@dataclass
class Criteria:
    K_rot_required: float = 50e6   # minimum rotational stiffness of the foundation [kNm/rad] = 50 GNm/rad (assumed)
    gamma_F: float = 1.35          # partial factor on unfavourable (wind) loads, ULS
    gamma_G_fav: float = 0.9       # partial factor on favourable permanent loads (overturning)
    contact_extreme_min: float = 0.5   # at least half of the base in contact under extreme loads
    fs_overturning: float = 1.0    # design stabilising / destabilising moment (with partial factors) >= 1
    fs_sliding: float = 1.0


# ---------------------------------------------------------------------- weights and analytical checks

def vertical_load(g: Geometry, L: Loads, s: Soil, c: Concrete):
    """Total characteristic vertical load at the founding level [kN]: turbine + concrete + backfill."""
    return L.N_turbine + c.unit_weight * g.concrete_volume() + s.gamma * g.backfill_volume()


def rigid_contact(N, M, R, n=400):
    """Linear soil pressure under a rigid circular footing on tensionless soil.

    Returns max pressure [kPa], contact area fraction, and the eccentricity e = M/N.
    Solves p(x) = max(0, p0 + k x) so that the integral of p equals N and of p*x equals M."""
    e = M / N
    A = math.pi * R**2
    if e <= R / 4:                                   # full contact (inside the kern of a circle)
        W = math.pi * R**3 / 4
        return N / A + M / W, 1.0, e
    x = np.linspace(-R, R, n)
    chord = 2 * np.sqrt(np.clip(R**2 - x**2, 0, None))
    dx = x[1] - x[0]

    def resid(c):                                    # neutral axis at x = c, pressure p = k (x - c) for x > c
        m = x > c
        F = np.sum(chord[m] * (x[m] - c)) * dx
        Mx = np.sum(chord[m] * (x[m] - c) * x[m]) * dx
        return Mx / F - e, F
    lo, hi = -R, R * 0.999
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        r_, _ = resid(mid)
        lo, hi = (mid, hi) if r_ < 0 else (lo, mid)
    c0 = 0.5 * (lo + hi)
    _, F = resid(c0)
    k = N / F
    pmax = k * (R - c0)
    contact = np.sum(chord[x > c0]) * dx / A
    return pmax, contact, e


def effective_area_pressure(N, M, R):
    """Bearing pressure on the effective (compressed) area of a circular footing, as used in DNV/Risoe
    guidelines and EC7: the area of the circle segment centred on the load point (eccentricity e)."""
    e = M / N
    if e >= R:
        return float("inf"), 0.0
    A_eff = 2 * (R**2 * math.acos(e / R) - e * math.sqrt(R**2 - e**2))
    b_e = 2 * (R - e)
    l_e = 2 * R * math.sqrt(1 - (1 - b_e / (2 * R)) ** 2)
    return N / A_eff, A_eff


def rotational_stiffness(g: Geometry, s: Soil):
    """Rotational stiffness of a circular footing on elastic soil [kNm/rad] (DNV/Risoe 2002 guideline):
    K = 8 G R^3 / (3 (1 - nu)) * (1 + R/(6H)) * (1 + 2D/R) * (1 + 0.7 D/H)."""
    R, D, H = g.R, g.D, s.H_layer
    return 8 * s.G * R**3 / (3 * (1 - s.nu)) * (1 + R / (6 * H)) * (1 + 2 * D / R) * (1 + 0.7 * D / H)


def subgrade_modulus(g: Geometry, s: Soil):
    """Winkler spring stiffness k [kPa/m] calibrated so that a rigid disc on Winkler springs has the same
    rotational stiffness as the elastic half-space formula (k * pi R^4 / 4 = K_rot)."""
    return rotational_stiffness(g, s) / (math.pi * g.R**4 / 4)


def analytical_checks(g: Geometry, L: Loads, s: Soil, c: Concrete, cr: Criteria = Criteria()):
    N = vertical_load(g, L, s, c)
    # moment at the founding level includes the horizontal force times the foundation depth
    M_n = L.M_normal + L.H_normal * g.D
    M_x = L.M_extreme + L.H_extreme * g.D
    _, contact_normal, e_n = rigid_contact(N, M_n, g.R)
    p_ext, contact_ext, e_x = rigid_contact(N, M_x, g.R)
    q_eff, A_eff = effective_area_pressure(N, M_x, g.R)
    K = rotational_stiffness(g, s)
    # overturning about the edge (design values): stabilising N*R (favourable) vs destabilising M
    fs_ot = cr.gamma_G_fav * N * g.R / (cr.gamma_F * M_x)
    fs_sl = cr.gamma_G_fav * N * s.mu / (cr.gamma_F * L.H_extreme)
    res = {
        "N_total_kN": N, "concrete_m3": g.concrete_volume(), "backfill_m3": g.backfill_volume(),
        "e_normal_m": e_n, "contact_normal": float(contact_normal), "e_extreme_m": e_x,
        "contact_extreme": float(contact_ext), "p_max_rigid_kPa": p_ext, "q_effective_kPa": q_eff,
        "K_rot_GNm_rad": K / 1e6, "overturning_ratio": fs_ot, "sliding_ratio": fs_sl,
    }
    res["ok_no_gapping_normal"] = bool(contact_normal >= 0.999)
    res["ok_contact_extreme"] = bool(contact_ext >= cr.contact_extreme_min)
    res["ok_bearing"] = bool(q_eff <= s.q_allow)
    res["ok_stiffness"] = bool(K >= cr.K_rot_required)
    res["ok_overturning"] = bool(fs_ot >= cr.fs_overturning)
    res["ok_sliding"] = bool(fs_sl >= cr.fs_sliding)
    return res


# ---------------------------------------------------------------------- Eurocode 2 reinforcement

def ec2_bending_reinforcement(m_Ed, h, c: Concrete, st: Steel, cover=0.075, bar=0.032):
    """Required reinforcement [mm2/m] for a design moment m_Ed [kNm/m] in a slab of thickness h [m],
    EC2 rectangular stress block, with the EC2 9.2.1.1 minimum."""
    d = h - cover - bar
    mu = m_Ed * 1e-3 / (1.0 * d**2 * c.fcd)                 # MNm/m / (m * m2 * MPa)
    mu = min(mu, 0.295)
    z = d * min(0.5 * (1 + math.sqrt(1 - 2 * mu)), 0.95)
    As = m_Ed * 1e-3 / (z * st.fyd) * 1e6                    # mm2/m
    As_min = max(0.26 * c.fctm / st.fyk, 0.0013) * d * 1e6   # mm2/m (b = 1 m)
    return max(As, As_min), As_min, d


def ec2_shear_resistance(h, As_mm2_per_m, c: Concrete, cover=0.075, bar=0.032):
    """Shear resistance without shear reinforcement v_Rd,c [kN/m] (EC2 6.2.2)."""
    d = h - cover - bar
    k = min(1 + math.sqrt(0.2 / d), 2.0)
    rho = min(As_mm2_per_m * 1e-6 / d, 0.02)
    v = max(0.18 / c.gamma_c * k * (100 * rho * c.fck) ** (1 / 3), 0.035 * k**1.5 * c.fck**0.5)   # MPa
    return v * d * 1e3                                       # kN/m
