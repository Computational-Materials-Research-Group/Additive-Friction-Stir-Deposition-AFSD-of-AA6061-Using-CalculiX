#!/usr/bin/env python3
"""
Additive Friction Stir Deposition (AFSD) - CalculiX model generator
====================================================================

Builds a fully coupled thermo-mechanical CalculiX deck for AFSD of
AA6061 onto an AA6061 substrate:

  * hollow rotating shoulder (outer radius R) with a square solid feedstock
  * deposited track = stadium (swept shoulder footprint, rounded ends)
  * deposit elements are pre-meshed, removed in step 1 and switched on
    progressively under the moving shoulder (*MODEL CHANGE, ADD)
  * new material enters at the feedstock exit temperature T_DEP and is
    stress free at that temperature (it shrinks as it cools -> residual stress)
  * frictional / plastic heat from the shoulder applied as a surface flux
    q(r) = 3 P r / (2 pi R^3) on the material under the shoulder
  * convection on every exposed face (recomputed each step), anvil
    contact conductance on the substrate bottom
  * clamped edges during deposition, cooling, then unclamping to show
    the final distortion and residual stress

Units: mm, s, N, tonne, MPa, mJ, degC
Run:   python3 afsd_ccx.py            -> afsd.inp     (coupled thermo-mechanical)
       python3 afsd_ccx.py thermal    -> afsd_th.inp  (temperature only, fast)
"""
import json
import numpy as np

# ===================== PARAMETERS =====================================
JOB = "afsd"
ANALYSIS = "coupled"          # "coupled" (thermal + stress) or "thermal"

# --- substrate plate ---
PLATE_L, PLATE_W, PLATE_T = 200.0, 100.0, 10.0   # x, y, z
NZ_SUB = [4.0, 3.0, 2.0, 1.0]                     # element heights, bottom -> top

# --- tool / feedstock ---
R_SH = 19.05              # shoulder outer radius (38.1 mm dia)
FEED_A = 12.7             # square feedstock side
RPM = 300.0               # spindle speed (used for the tool animation note)

# --- process ---
N_LAYERS = 1              # deposited layers (alternate direction each layer)
LAYER_H = 2.0             # layer thickness
X_START, X_END = 45.0, 155.0   # tool-centre travel path along x (y = 0)
V_TRAVEL = 2.0            # traverse speed  mm/s  (120 mm/min)
DX_STEP = 2.0             # tool advance per analysis step (= frame)
T_PLUNGE = 4.0            # dwell at start of each layer (s)
T_INTERLAYER = 20.0       # tool lifted between layers (s)
P_TOOL = 2400.0           # heat into the work from the shoulder (W)
T_DEP = 380.0             # temperature of material leaving the tool (degC)
T_AMB = 20.0

# --- cooling / unclamp ---
COOL_TIMES = [2, 3, 5, 8, 12, 20, 30, 45, 60, 90, 125, 200]   # step lengths (s)
UNCLAMP = True

# --- heat losses (W/m^2K -> mW/mm^2K: x 1e-3) ---
H_AIR = 15e-3             # free convection (+ lumped radiation)
H_ANVIL = 300e-3          # substrate bottom to steel anvil
CLAMP_Y = 40.0            # |y| >= CLAMP_Y on top surface is clamped

# --- mesh ---
SYMMETRY = True           # half model (y >= 0), mirrored back for ParaView
H_FINE = 2.0              # in-plane element size around the track
H_COARSE_X, H_COARSE_Y = 8.0, 6.5

# --- AA6061 (approximate, temperature dependent) ---
MAT_T = [20, 100, 200, 300, 400, 500, 580]
E_MOD = [68.9e3, 66.2e3, 59.2e3, 47.8e3, 31.7e3, 15.0e3, 5.0e3]
YIELD = {20: 276, 100: 262, 150: 240, 200: 200, 250: 140, 300: 90,
         350: 50, 400: 30, 450: 18, 500: 12, 580: 5}
COND = {20: 167, 200: 180, 400: 195, 580: 205}               # W/mK = N/(s K)
CP = {20: 896e6, 200: 960e6, 400: 1050e6, 580: 1120e6}       # mm^2/(s^2 K)
ALPHA = {20: 23.4e-6, 200: 24.9e-6, 400: 27.5e-6, 580: 29.0e-6}
RHO = 2.7e-9
NU = 0.33
# ======================================================================


def axis_nodes(lo, hi, f_lo, f_hi, h_fine, h_coarse):
    """1D node coordinates: fine spacing in [f_lo, f_hi], coarse outside."""
    def seg(a, b, h):
        n = max(1, int(round(abs(b - a) / h)))
        return np.linspace(a, b, n + 1)
    pts = [np.array([lo])]
    for a, b, h in ((lo, f_lo, h_coarse), (f_lo, f_hi, h_fine), (f_hi, hi, h_coarse)):
        if b - a > 1e-9:
            pts.append(seg(a, b, h)[1:])
    return np.concatenate(pts)


def in_stadium(x, y, xa, xb, r):
    """Points inside the swept shoulder footprint between tool x = xa..xb."""
    lo, hi = min(xa, xb), max(xa, xb)
    xc = np.clip(x, lo, hi)
    return (x - xc) ** 2 + y ** 2 <= r ** 2 + 1e-9


def build():
    # ------------------------------------------------------------ grid
    f_lo = round(min(X_START, X_END) - R_SH - 2.0)
    f_hi = round(max(X_START, X_END) + R_SH + 2.0)
    xs = axis_nodes(0.0, PLATE_L, f_lo, f_hi, H_FINE, H_COARSE_X)
    yf = round(R_SH + 5.0)
    ys = (axis_nodes(0.0, PLATE_W / 2, 0.0, yf, H_FINE, H_COARSE_Y) if SYMMETRY else
          axis_nodes(-PLATE_W / 2, PLATE_W / 2, -yf, yf, H_FINE, H_COARSE_Y))
    zs = np.concatenate([[0.0], np.cumsum(NZ_SUB)])
    zs = zs * PLATE_T / zs[-1]
    nsub = len(zs) - 1
    zs = np.concatenate([zs, PLATE_T + LAYER_H * np.arange(1, N_LAYERS + 1)])
    nx, ny, nz = len(xs), len(ys), len(zs)

    nid = lambda i, j, k: 1 + i + nx * (j + ny * k)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    coords = np.zeros((nx * ny * nz, 3))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                coords[nid(i, j, k) - 1] = (xs[i], ys[j], zs[k])

    # ------------------------------------------------------------ elements
    xc = 0.5 * (xs[:-1] + xs[1:])
    yc = 0.5 * (ys[:-1] + ys[1:])
    track = in_stadium(xc[:, None], yc[None, :], X_START, X_END, R_SH)

    elems, eijk, region = [], [], []
    for k in range(nz - 1):
        layer = 0 if k < nsub else k - nsub + 1
        for j in range(ny - 1):
            for i in range(nx - 1):
                if layer > 0 and not track[i, j]:
                    continue
                n = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                     nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                elems.append(n)
                eijk.append((i, j, k))
                region.append(layer)
    elems = np.array(elems)
    eijk = np.array(eijk)
    region = np.array(region)
    ne = len(elems)
    eid = np.arange(1, ne + 1)
    lookup = {tuple(e): n for e, n in zip(map(tuple, eijk), eid)}

    used = np.unique(elems)
    sub_nodes = np.unique(elems[region == 0])
    dep_only = np.setdiff1d(used, sub_nodes)

    # ------------------------------------------------------------ step schedule
    steps = []   # dicts: name, dt, tool(x,y,z or None), add(list), heat(bool)
    t_dir = []
    active = set(eid[region == 0])
    for L in range(1, N_LAYERS + 1):
        a, b = (X_START, X_END) if L % 2 == 1 else (X_END, X_START)
        ztop = PLATE_T + L * LAYER_H
        if L > 1:
            steps.append(dict(name=f"Layer {L-1} interpass dwell", dt=T_INTERLAYER,
                              tool=None, add=[], layer=L, xt=None))
        nstep = int(round(abs(b - a) / DX_STEP))
        path = np.linspace(a, b, nstep + 1)
        lay = (region == L)
        cen_x = xc[eijk[:, 0]]
        cen_y = yc[eijk[:, 1]]
        for s, xt in enumerate(path):
            cover = lay & in_stadium(cen_x, cen_y, a, xt, R_SH)
            new = [e for e in eid[cover] if e not in active]
            active.update(new)
            steps.append(dict(name=(f"Layer {L} plunge" if s == 0 else f"Layer {L} pass"),
                              dt=(T_PLUNGE if s == 0 else DX_STEP / V_TRAVEL),
                              tool=(xt, 0.0, ztop), add=new, layer=L, xt=xt))
    for i, dt in enumerate(COOL_TIMES):
        steps.append(dict(name="Cooling", dt=float(dt), tool=None, add=[], layer=N_LAYERS, xt=None))
    if UNCLAMP and ANALYSIS == "coupled":
        steps.append(dict(name="Unclamped", dt=1.0, tool=None, add=[], layer=N_LAYERS,
                          xt=None, unclamp=True))
    return dict(xs=xs, ys=ys, zs=zs, nsub=nsub, coords=coords, elems=elems, eijk=eijk,
                region=region, lookup=lookup, used=used, dep_only=dep_only, steps=steps,
                xc=xc, yc=yc)


# C3D8 faces: S1 bottom(-z) S2 top(+z) S3 -y S4 +x S5 +y S6 -x
FACE_DIRS = {1: (0, 0, -1), 2: (0, 0, 1), 3: (0, -1, 0), 4: (1, 0, 0), 5: (0, 1, 0), 6: (-1, 0, 0)}


def table(d):
    return "".join(f"{v:.6g},{t}\n" for t, v in d.items())


def write_deck(m):
    xs, ys, coords, elems, eijk, region = m["xs"], m["ys"], m["coords"], m["elems"], m["eijk"], m["region"]
    lookup, steps = m["lookup"], m["steps"]
    ne = len(elems)
    eid = np.arange(1, ne + 1)
    xc, yc = m["xc"], m["yc"]
    dxs, dys = np.diff(xs), np.diff(ys)
    coupled = ANALYSIS == "coupled"

    out = []
    w = out.append
    w("** Additive Friction Stir Deposition (AFSD) - generated by afsd_ccx.py\n")
    w("** units: mm, s, N, tonne, MPa, mJ, degC\n")
    w("*NODE, NSET=NALL\n")
    for n in m["used"]:
        x, y, z = coords[n - 1]
        w(f"{n},{x:.6f},{y:.6f},{z:.6f}\n")
    w("*ELEMENT, TYPE=C3D8, ELSET=EALL\n")
    for e, conn in zip(eid, elems):
        w(f"{e}," + ",".join(map(str, conn)) + "\n")

    def elset(name, ids):
        w(f"*ELSET, ELSET={name}\n")
        ids = list(ids)
        for i in range(0, len(ids), 16):
            w(",".join(map(str, ids[i:i + 16])) + ",\n")

    def nset(name, ids):
        w(f"*NSET, NSET={name}\n")
        ids = list(ids)
        for i in range(0, len(ids), 16):
            w(",".join(map(str, ids[i:i + 16])) + ",\n")

    elset("SUB", eid[region == 0])
    for L in range(1, N_LAYERS + 1):
        elset(f"DEP{L}", eid[region == L])
    elset("DEPALL", eid[region > 0])

    # boundary node sets
    z = coords[:, 2]
    y = coords[:, 1]
    used = m["used"]
    bottom = [n for n in used if abs(z[n - 1]) < 1e-9]
    clamp = [n for n in used if abs(z[n - 1] - PLATE_T) < 1e-9 and abs(y[n - 1]) >= CLAMP_Y - 1e-9]
    nset("NBOT", bottom)
    nset("NCLAMP", clamp)
    if SYMMETRY:
        nset("NSYM", [n for n in used if abs(y[n - 1]) < 1e-9])
    # 3-2-1 support for the unclamped state (three bottom corners)
    x = coords[:, 0]
    def corner(xv, yv):
        return int(min(bottom, key=lambda n: (x[n-1]-xv)**2 + (y[n-1]-yv)**2))
    y0 = 0.0 if SYMMETRY else -PLATE_W / 2
    c1, c2, c3 = corner(0, y0), corner(PLATE_L, y0), corner(0, PLATE_W/2)

    # ---------------------------------------------------- material
    w("*MATERIAL, NAME=AA6061\n")
    if coupled:
        w("*ELASTIC\n")
        for t, e in zip(MAT_T, E_MOD):
            w(f"{e},{NU},{t}\n")
        w("*PLASTIC\n")
        for t, sy in YIELD.items():
            w(f"{sy},0.0,{t}\n{sy*1.15:.4g},0.2,{t}\n{sy*1.2:.4g},2.0,{t}\n")
        w(f"*EXPANSION, ZERO={T_AMB}\n" + table(ALPHA))
    w("*CONDUCTIVITY\n" + table(COND))
    w("*SPECIFIC HEAT\n" + table(CP))
    w(f"*DENSITY\n{RHO}\n")
    w("*SOLID SECTION, ELSET=EALL, MATERIAL=AA6061\n")

    # ---------------------------------------------------- initial temperature
    w("*INITIAL CONDITIONS, TYPE=TEMPERATURE\n")
    w(f"NALL,{T_AMB}\n")
    for n in m["dep_only"]:
        w(f"{n},{T_DEP}\n")   # inactive nodes keep this until activated

    # ---------------------------------------------------- per-step loads
    active = np.zeros(ne + 1, bool)
    active[eid[region == 0]] = True
    ijk_active = lambda i, j, k: (lookup.get((i, j, k), 0))
    t_now = 0.0
    sched = []
    for s_idx, st in enumerate(steps):
        for e in st["add"]:
            active[e] = True
        tool = st["tool"]
        # faces under the shoulder (top faces of topmost active element of the current layer)
        heated, flux_faces = set(), []
        if tool is not None:
            xt, yt, ztop = tool
            L = st["layer"]
            k_top = m["nsub"] + L - 1
            for e in eid[(region == L)]:
                if not active[e]:
                    continue
                i, j, k = eijk[e - 1]
                r = np.hypot(xc[i] - xt, yc[j] - yt)
                if r <= R_SH:
                    flux_faces.append((e, r, dxs[i] * dys[j]))
                    heated.add(e)
            if flux_faces:
                r = np.array([f[1] for f in flux_faces]); A = np.array([f[2] for f in flux_faces])
                p_model = P_TOOL * 1e3 * (0.5 if SYMMETRY else 1.0)     # mW in the model
                q = 3 * P_TOOL * 1e3 * r / (2 * np.pi * R_SH ** 3)       # mW/mm^2, q ~ r
                q *= p_model / max(np.sum(q * A), 1e-12)                # exact total power
                flux_faces = [(f[0], qq) for f, qq in zip(flux_faces, q)]

        # exposed faces for convection
        film = []
        for e in eid[active[1:]]:
            i, j, k = eijk[e - 1]
            for f, (di, dj, dk) in FACE_DIRS.items():
                nb = ijk_active(i + di, j + dj, k + dk)
                if nb and active[nb]:
                    continue
                if SYMMETRY and f == 3 and j == 0:
                    continue                      # symmetry plane is adiabatic
                if f == 1:
                    film.append((e, f, H_ANVIL))
                elif f == 2 and e in heated:
                    continue
                else:
                    film.append((e, f, H_AIR))

        dt = st["dt"]
        unclamp = st.get("unclamp", False)
        w(f"** ---------------- step {s_idx+1}: {st['name']}\n")
        w("*STEP, INC=1000\n")
        if coupled:
            w("*COUPLED TEMPERATURE-DISPLACEMENT\n")
        else:
            w("*HEAT TRANSFER\n")
        w(f"{dt/2:.6g},{dt:.6g},{dt*1e-5:.6g},{dt:.6g}\n" if dt > 1.5 else f"{dt:.6g},{dt:.6g},{dt*1e-5:.6g},{dt:.6g}\n")
        if s_idx == 0:
            w("*MODEL CHANGE, TYPE=ELEMENT, REMOVE\nDEPALL\n")
        if st["add"]:
            w("*MODEL CHANGE, TYPE=ELEMENT, ADD\n")
            for i in range(0, len(st["add"]), 16):
                w(",".join(map(str, st["add"][i:i + 16])) + ",\n")
        if coupled and (s_idx == 0 or unclamp):
            w("*BOUNDARY, OP=NEW\n")
            if unclamp:
                w(f"{c1},1,1\n{c1},3,3\n{c2},3,3\n{c3},3,3\n" if SYMMETRY else
                  f"{c1},1,3\n{c2},2,3\n{c3},3,3\n")
            else:
                w("NBOT,3,3\nNCLAMP,1,3\n")
            if SYMMETRY:
                w("NSYM,2,2\n")
        w("*DFLUX, OP=NEW\n")
        for e, qq in flux_faces:
            w(f"{e},S2,{qq:.6g}\n")
        w("*FILM, OP=NEW\n")
        for e, f, h in film:
            w(f"{e},F{f},{T_AMB},{h}\n")
        w("*NODE FILE, FREQUENCY=10000\nNT" + (",U\n" if coupled else "\n"))
        if coupled:
            w("*EL FILE, FREQUENCY=10000\nS,PEEQ\n")
        w("*END STEP\n")
        t_now += dt
        sched.append(dict(step=s_idx + 1, name=st["name"], t_end=t_now,
                          tool=None if tool is None else list(map(float, tool)),
                          layer=st["layer"], unclamp=unclamp,
                          n_active=int(active[1:].sum())))

    with open(JOB + ".inp", "w") as fh:
        fh.write("".join(out))

    # activation step for each element (0 = substrate, always active)
    act_step = np.zeros(ne, int)
    for s_idx, st in enumerate(steps):
        for e in st["add"]:
            act_step[e - 1] = s_idx + 1
    np.savez_compressed(JOB + "_model.npz", coords=coords, elems=elems, region=region,
                        act_step=act_step)
    meta = dict(job=JOB, analysis=ANALYSIS, steps=sched, R_SH=R_SH, FEED_A=FEED_A,
                LAYER_H=LAYER_H, PLATE_T=PLATE_T, SYMMETRY=SYMMETRY, N_LAYERS=N_LAYERS, RPM=RPM,
                V_TRAVEL=V_TRAVEL, X_START=X_START, X_END=X_END, T_AMB=T_AMB)
    with open(JOB + "_schedule.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    return ne, len(m["used"]), len(steps)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "thermal":     # python3 afsd_ccx.py thermal
        ANALYSIS, JOB = "thermal", "afsd_th"
    model = build()
    ne, nn, ns = write_deck(model)
    print(f"{JOB}.inp written: {nn} nodes, {ne} elements, {ns} steps "
          f"({ANALYSIS}, {N_LAYERS} layer(s), P = {P_TOOL:.0f} W, v = {V_TRAVEL} mm/s)")
