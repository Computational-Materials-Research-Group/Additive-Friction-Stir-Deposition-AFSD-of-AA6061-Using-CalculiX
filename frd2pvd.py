#!/usr/bin/env python3
"""
frd2pvd.py - convert the AFSD CalculiX result (.frd) into a ParaView animation

Writes, inside ./pvd/ :
  afsd_part.pvd + afsd_part_XXXX.vtu   substrate + deposit, only the elements
                                       active at that moment (the layer grows)
  afsd_tool.pvd + afsd_tool_XXXX.vtu   hollow rotating shoulder + square
                                       feedstock moving along the path
and in the working directory:
  afsd_probes.csv / afsd_probes.png    thermocouple-style temperature histories
  afsd_frames.csv                      time, tool position, peak T for each frame

Usage:  python3 frd2pvd.py [job]        (default job = afsd)
Only numpy is needed (no vtk / ccx2paraview).
"""
import base64
import csv
import json
import os
import sys
import zlib

import numpy as np

JOB = sys.argv[1] if len(sys.argv) > 1 else "afsd"
OUT = "pvd"
VIS_SPIN_DEG = 24.0     # visual shoulder rotation per frame (true 300 rpm would alias)


# ------------------------------------------------------------------ FRD reader
def read_frd(path):
    """Yield (time, {field: (node_ids, values)}) for every result frame."""
    frames = []
    cur_time, cur = None, {}
    field, ids, vals, ncomp = None, [], [], 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            key = line[:6]
            if key == "  100C":
                parts = line.split()
                t = float(parts[2])
                if cur_time is not None and abs(t - cur_time) > 1e-9:
                    frames.append((cur_time, cur)); cur = {}
                cur_time = t
            elif line.startswith(" -4"):
                field = line.split()[1]; ids, vals = [], []
                ncomp = int(line.split()[2])
            elif line.startswith(" -1") and field is not None:
                ids.append(int(line[3:13]))
                s = line[13:].rstrip("\n")
                nv = min(len(s.rstrip()) // 12, 6)
                vals.append([float(s[12 * c:12 * c + 12]) for c in range(nv)])
            elif line.startswith(" -3") and field is not None:
                if field in ("NDTEMP", "DISP", "STRESS", "PE"):
                    cur[field] = (np.array(ids), np.array(vals))
                field = None
    if cur_time is not None:
        frames.append((cur_time, cur))
    return frames


# ------------------------------------------------------------------ VTU writer
def _b64(arr):
    raw = np.ascontiguousarray(arr).tobytes()
    comp = zlib.compress(raw, 6)
    header = np.array([1, len(raw), len(raw), len(comp)], dtype=np.uint32).tobytes()
    return base64.b64encode(header).decode() + base64.b64encode(comp).decode()


def _vtype(a):
    return {np.dtype("float32"): "Float32", np.dtype("float64"): "Float64",
            np.dtype("int32"): "Int32", np.dtype("int64"): "Int64",
            np.dtype("uint8"): "UInt8"}[a.dtype]


def write_vtu(fname, pts, conn, ctype, pdata, cdata):
    """Unstructured grid with one cell type (12 = hexahedron)."""
    nc, npc = conn.shape
    offs = np.arange(1, nc + 1, dtype=np.int64) * npc
    types = np.full(nc, ctype, np.uint8)

    def arr(name, a, ncomp=None):
        a = np.asarray(a)
        nc_ = ncomp or (1 if a.ndim == 1 else a.shape[1])
        extra = ""
        if nc_ == 6:
            extra = (' ComponentName0="XX" ComponentName1="YY" ComponentName2="ZZ"'
                     ' ComponentName3="XY" ComponentName4="YZ" ComponentName5="XZ"')
        return (f'<DataArray type="{_vtype(a)}" Name="{name}" NumberOfComponents="{nc_}"{extra} '
                f'format="binary">{_b64(a)}</DataArray>\n')

    with open(fname, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" '
                'byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n'
                f'<UnstructuredGrid><Piece NumberOfPoints="{len(pts)}" NumberOfCells="{nc}">\n')
        f.write("<Points>" + arr("Points", pts.astype(np.float32), 3) + "</Points>\n")
        f.write("<Cells>" + arr("connectivity", conn.astype(np.int64).ravel(), 1)
                + arr("offsets", offs, 1) + arr("types", types, 1) + "</Cells>\n")
        f.write("<PointData>")
        for k, v in pdata.items():
            f.write(arr(k, v))
        f.write("</PointData>\n<CellData>")
        for k, v in cdata.items():
            f.write(arr(k, v))
        f.write("</CellData>\n</Piece></UnstructuredGrid></VTKFile>\n")


def write_pvd(fname, entries):
    with open(fname, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n<Collection>\n')
        for t, file in entries:
            f.write(f'<DataSet timestep="{t:.6g}" group="" part="0" file="{file}"/>\n')
        f.write("</Collection>\n</VTKFile>\n")


# ------------------------------------------------------------------ tool geometry
def tool_mesh(R, a, h_sh=22.0, feed_len=100.0, nseg=48):
    """Hex mesh of a hollow shoulder (square bore of side a) + square feedstock bar.
    Local frame: tool axis = z, shoulder face at z = 0."""
    pts, cells, part, stripe = [], [], [], []

    def add_pt(p):
        pts.append(p); return len(pts) - 1
    th = np.linspace(0, 2 * np.pi, nseg, endpoint=False)
    zs = [0.0, h_sh]
    rr = [0.0, 0.5, 1.0]                                   # bore -> outer
    idx = {}
    for it, t in enumerate(th):
        c, s = np.cos(t), np.sin(t)
        r_in = (a / 2) / max(abs(c), abs(s))
        for ir, f in enumerate(rr):
            r = r_in + f * (R - r_in)
            for iz, z in enumerate(zs):
                idx[it, ir, iz] = add_pt((r * c, r * s, z))
    for it in range(nseg):
        jt = (it + 1) % nseg
        for ir in range(len(rr) - 1):
            cells.append([idx[it, ir, 0], idx[jt, ir, 0], idx[jt, ir + 1, 0], idx[it, ir + 1, 0],
                          idx[it, ir, 1], idx[jt, ir, 1], idx[jt, ir + 1, 1], idx[it, ir + 1, 1]])
            part.append(0); stripe.append((it // 4) % 2)
    # feedstock bar (2 hexes stacked: inside the bore and above it)
    h = a / 2
    zb = [0.0, h_sh, feed_len]
    ring = []
    for z in zb:
        ring.append([add_pt((-h, -h, z)), add_pt((h, -h, z)), add_pt((h, h, z)), add_pt((-h, h, z))])
    for k in range(2):
        cells.append(ring[k] + ring[k + 1]); part.append(1); stripe.append(2)
    return np.array(pts), np.array(cells), np.array(part), np.array(stripe), ring


def main():
    meta = json.load(open(f"{JOB}_schedule.json"))
    mod = np.load(f"{JOB}_model.npz")
    coords, elems, region, act_step = mod["coords"], mod["elems"], mod["region"], mod["act_step"]
    steps = meta["steps"]
    t_ends = np.array([s["t_end"] for s in steps])
    T_amb = meta["T_AMB"]

    print(f"reading {JOB}.frd ...")
    frames = read_frd(f"{JOB}.frd")
    print(f"  {len(frames)} result frames")
    os.makedirs(OUT, exist_ok=True)

    # tool template
    tp, tc, tpart, tstripe, ring = tool_mesh(meta["R_SH"], meta["FEED_A"])
    feed_top_ids = ring[2]
    dep_area_per_mm = 2 * meta["R_SH"] * meta["LAYER_H"]
    feed_area = meta["FEED_A"] ** 2
    feed_per_mm = dep_area_per_mm / feed_area             # feedstock consumed per mm travelled

    # probes (nearest nodes): substrate under track mid-length at several depths + outside
    xm = 0.5 * (meta["X_START"] + meta["X_END"])
    zt = meta["PLATE_T"]
    probe_xyz = {"P1_interface_centre": (xm, 0, zt), "P2_substrate_3mm_deep": (xm, 0, zt - 3.0),
                 "P3_substrate_bottom": (xm, 0, 0.0), "P4_beside_track_y30": (xm, 30.0, zt),
                 "P5_deposit_top": (xm, 0, zt + meta["N_LAYERS"] * meta["LAYER_H"])}
    probes = {k: int(np.argmin(np.sum((coords - np.array(v)) ** 2, axis=1))) + 1
              for k, v in probe_xyz.items()}
    probe_rows = []

    part_entries, tool_entries, frame_rows = [], [], []
    last_tool, last_tool_fi = None, 0
    travelled = 0.0
    prev_x = None
    nn = len(coords)
    for fi, (t, data) in enumerate(frames):
        s_idx = int(np.searchsorted(t_ends, t - 1e-6))     # step index (0-based)
        s_idx = min(s_idx, len(steps) - 1)
        st = steps[s_idx]
        alive = act_step <= st["step"]
        con = elems[alive]
        used = np.unique(con)
        remap = np.zeros(nn + 1, np.int64); remap[used] = np.arange(len(used))
        conn = remap[con]

        def nodal(field, ncomp, fill=0.0):
            a = np.full((nn + 1, ncomp), fill, dtype=np.float64)
            if field in data:
                ids, v = data[field]
                a[ids, :] = v[:, :ncomp]
            return a[used]

        T = nodal("NDTEMP", 1, T_amb)[:, 0]
        pdata = {"T": T.astype(np.float32)}
        if "DISP" in data:
            U = nodal("DISP", 3)
            S = nodal("STRESS", 6)
            sx, sy, sz, txy, tyz, txz = S.T
            mises = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
                            + 3 * (txy ** 2 + tyz ** 2 + txz ** 2))
            pdata.update(U=U.astype(np.float32), S=S.astype(np.float32),
                         S_Mises=mises.astype(np.float32))
            if "PE" in data:
                pdata["PEEQ"] = nodal("PE", 1)[:, 0].astype(np.float32)
        cdata = {"Region": region[alive].astype(np.int32),
                 "ActivationStep": act_step[alive].astype(np.int32)}
        pts = coords[used - 1]
        if meta.get("SYMMETRY", False):              # mirror the half model about y = 0
            npt = len(pts)
            mp = pts.copy(); mp[:, 1] *= -1
            pts = np.vstack([pts, mp])
            conn = np.vstack([conn, conn[:, [0, 3, 2, 1, 4, 7, 6, 5]] + npt])
            for k, v in list(pdata.items()):
                m = v.copy()
                if k == "U":
                    m[:, 1] *= -1
                elif k == "S":
                    m[:, 3] *= -1; m[:, 4] *= -1     # XY and YZ change sign
                pdata[k] = np.concatenate([v, m])
            cdata = {k: np.concatenate([v, v]) for k, v in cdata.items()}
        fname = f"{JOB}_part_{fi:04d}.vtu"
        write_vtu(os.path.join(OUT, fname), pts, conn, 12, pdata, cdata)
        part_entries.append((t, fname))

        # ---- tool position
        if st["tool"] is not None:
            xt, yt, ztop = st["tool"]
            if prev_x is not None:
                travelled += abs(xt - prev_x)
            prev_x = xt
            last_tool = (xt, yt, ztop)
            last_tool_fi = fi
            lift = 0.0
        else:
            prev_x = None
            if last_tool is None:
                continue
            xt, yt, ztop = last_tool
            lift = min(60.0, 15.0 * (fi - last_tool_fi))
        ang = np.radians(VIS_SPIN_DEG * fi) * (0 if lift > 0 else 1)
        p = tp.copy()
        p[feed_top_ids, 2] = max(30.0, 100.0 - feed_per_mm * travelled)
        ca, sa = np.cos(ang), np.sin(ang)
        p[:, 0], p[:, 1] = ca * tp[:, 0] - sa * tp[:, 1], sa * tp[:, 0] + ca * tp[:, 1]
        p += np.array([xt, yt, ztop + lift])
        tname = f"{JOB}_tool_{fi:04d}.vtu"
        write_vtu(os.path.join(OUT, tname), p, tc, 12, {},
                  {"ToolPart": tpart.astype(np.int32), "Stripe": tstripe.astype(np.int32)})
        tool_entries.append((t, tname))

        Tall = np.full(nn + 1, np.nan)
        if "NDTEMP" in data:
            ids, v = data["NDTEMP"]; Tall[ids] = v[:, 0]
        probe_rows.append([t] + [Tall[n] for n in probes.values()])
        frame_rows.append([fi, t, st["name"], "" if st["tool"] is None else st["tool"][0],
                           int(alive.sum()), float(np.nanmax(T)),
                           float(pdata["S_Mises"].max()) if "S_Mises" in pdata else ""])

    write_pvd(os.path.join(OUT, f"{JOB}_part.pvd"), part_entries)
    # tool exists only while frames have a tool; pad with first tool pose before plunge
    write_pvd(os.path.join(OUT, f"{JOB}_tool.pvd"), tool_entries)

    with open(f"{JOB}_frames.csv", "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["frame", "time_s", "stage", "tool_x_mm", "active_elements", "Tmax_C", "Mises_max_MPa"])
        wcsv.writerows(frame_rows)
    with open(f"{JOB}_probes.csv", "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["time_s"] + list(probes.keys()))
        wcsv.writerows(probe_rows)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        a = np.array(probe_rows, float)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for k, name in enumerate(probes):
            ax.plot(a[:, 0], a[:, k + 1], marker=".", label=name.replace("_", " "))
        ax.set_xlabel("time (s)"); ax.set_ylabel("temperature (°C)")
        ax.set_title("AFSD thermal history at probe points (mid-track)")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(f"{JOB}_probes.png", dpi=150)
    except ImportError:
        pass
    print(f"  wrote {len(part_entries)} part frames, {len(tool_entries)} tool frames in ./{OUT}")
    print(f"  peak temperature {max(r[5] for r in frame_rows):.0f} °C")


if __name__ == "__main__":
    main()
