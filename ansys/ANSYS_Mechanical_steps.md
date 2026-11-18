# Repeating the foundation model in ANSYS Mechanical

This is the 3D version of the Python model, as it would be built in industry. Use the reference geometry from `src/foundation.py`. ANSYS Student limits the model to 128 000 nodes, so use a coarse soil mesh, or model half the foundation with a symmetry plane through the load direction.

## 1. Geometry (SpaceClaim)
1. **Foundation**: revolve the cross-section. The radius is 10 m, the edge is 1.0 m thick, the slab is 3.2 m thick at the pedestal, and the pedestal radius is 3.0 m with a conical top in between.
2. **Anchor cage ("foundation basket")**: two rings of radius 2.2 m (bottom anchor flange 0.6 m wide at 0.3 m above the base, top flange at the pedestal top), connected by line bodies for the bolts (e.g. 2 × 80 M42 bolts).
3. **Soil block**: a cylinder of 3–4 × R in radius and depth. Subtract the foundation from it and add the backfill on top of the slab as a separate body.
4. Use **Share Topology** for the backfill/soil, but **not** between the foundation and the soil, because they are connected by contact.

## 2. Materials
| Body | Model |
|---|---|
| Concrete C35/45 | Linear elastic: E = 34 GPa, ν = 0.2. For cracking, use the *Concrete* model of Menetrey–Willam, or check the stresses against f_ctm |
| Bolts, flanges | Structural steel, E = 210 GPa. Bolt pretension through **Bolt Pretension** |
| Soil | Elastic (G = 50 MPa, ν = 0.3) or **Mohr–Coulomb** (φ' = 35°, c' = 0) for non-linear analyses |
| Reinforcement | **Reinforcement** element (REINF264) as smeared layers, with the bar areas from `results/02_reference_design.csv` |

## 3. Contacts
- Foundation base to soil: **Frictional** contact, μ = 0.5, Augmented Lagrange. This captures gapping (lift-off) automatically.
- Foundation sides to backfill: frictional or frictionless.
- Anchor flanges to concrete: **Bonded**. Bolts to flanges: bonded line-to-surface (beam connection).

## 4. Loads and supports
1. **Standard Earth Gravity** (self-weight of the concrete and backfill).
2. **Remote Force/Moment** on the top flange of the anchor cage: Fz = −5 000 kN, Fx = 1 200 kN, My = 90 000 kNm (extreme). Repeat with My = 55 000 kNm and Fx = 800 kN (normal production).
3. Supports: soil block bottom **Fixed**, sides **Frictionless** (or a larger block with roller sides).
4. Turn on **Large Deflection** only if Mohr–Coulomb soil is used. Run gravity in step 1 and the turbine loads in step 2.

## 5. Results to compare with the Python model
| Result | Python model (reference design) |
|---|---|
| Contact status / gap at the base, extreme load | 87 % in contact |
| Max soil pressure, extreme load | ≈ 192 kPa |
| Tilt of the foundation | ≈ 0.3 mrad |
| Radial bending moment at the pedestal face (ULS) | ≈ 5 300 kNm/m (from a path result of the stresses: integrate σ_r·z over the thickness) |
| Bolt force on the tension side | ≈ 7 950 kN/m of bolt circle |

## 6. Parametric study and Python automation
1. Promote R, h_e and h_c to **Parameters** (tick the P box in SpaceClaim) and solve in a **Design Points** table.
2. Use **Response Surface Optimization** (minimise concrete volume, subject to contact ≥ 50 % and bearing ≤ q_allow), or drive Workbench from Python with **PyAnsys** (`ansys-mechanical-core`) and reuse `scripts/03_parametric_optimisation.py` as the outer loop.
