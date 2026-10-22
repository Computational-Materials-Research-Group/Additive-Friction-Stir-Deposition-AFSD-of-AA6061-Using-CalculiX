# Additive Friction Stir Deposition (AFSD) of AA6061 Using CalculiX

<p align="center">
  <img src="https://img.shields.io/badge/CalculiX-2.21-blue?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/AFSD-Solid%20State%20AM-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Thermo--Mechanical-Coupled-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/AA6061-Aluminium-silver?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Element%20Activation-MODEL%20CHANGE-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Residual%20Stress-Distortion-darkgreen?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/ParaView-PVD%20Animation-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.23558007"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23558007.svg" alt="DOI"/></a>
</p>

<p align="center">
  A <b>fully scripted, open-source process model</b> of
  <b>additive friction stir deposition (AFSD)</b>: a square
  <b>AA6061 feedstock</b> is fed through a <b>hollow rotating shoulder</b>
  and deposited as a layer on an <b>AA6061 substrate</b>. The model
  captures the <b>moving frictional heat source</b>, the
  <b>progressive growth of the deposited layer</b>, the
  <b>thermal history</b> of the build and substrate, and the
  <b>residual stress and distortion</b> after cooling and unclamping, and
  exports a <b>69-frame ParaView animation</b> with the rotating tool and
  shrinking feedstock visible.
</p>

<img width="1008" height="526" alt="afsd calculix" src="https://github.com/user-attachments/assets/bcedd0eb-ca38-4bdd-b61c-64ab548176d1" />

---

## Concept

Additive friction stir deposition is a **solid-state** metal additive
manufacturing process. A solid feedstock bar is pushed through a hollow,
non-consumable rotating shoulder. Friction and severe plastic
deformation heat the material to roughly 60 to 90 % of its melting
temperature, so it softens and flows without melting. The shoulder
forges the plasticised material onto the substrate while it travels,
leaving a deposited layer as wide as the shoulder. Repeating this layer
by layer builds a near-net-shape part with a fine, recrystallised,
fully dense microstructure and none of the porosity or hot cracking of
fusion-based processes.

The process is governed by heat. The temperature under the shoulder
controls how easily the material flows; the thermal cycle in the
deposit and substrate controls precipitate dissolution and overageing in
heat-treatable alloys such as AA6061; and the uneven heating and
cooling, restrained by the clamps, leave **residual stress** and
**distortion** in the finished part.

This repository is a compact testbed for that problem. One Python
script generates the complete CalculiX model (mesh, deposition path,
element activation schedule, heat input, heat losses, clamping and
cooling); the solver writes one frame per step; and a numpy-only
converter builds a ParaView animation that shows only the material
deposited so far, together with the moving tool.

| Aspect | Description |
|---|---|
| Process | Single-layer AFSD, 110 mm straight track at 2 mm/s |
| Tool | Hollow shoulder 38.1 mm diameter, square feedstock 12.7 mm |
| Substrate | AA6061 plate, 200 x 100 x 10 mm, edges clamped |
| Deposit | 2 mm layer, stadium footprint (swept shoulder) |
| Material | AA6061, temperature-dependent thermal, elastic and plastic properties |
| Physics | Transient heat conduction coupled with quasi-static elastoplasticity |
| Deposition | Pre-meshed layer, switched on under the tool (`*MODEL CHANGE`) |
| Heat source | Shoulder surface flux, q(r) proportional to r, 2400 W |
| Elements | C3D8 (8-node hexahedron), half model with symmetry |
| Solver | CalculiX `ccx` 2.21, `*COUPLED TEMPERATURE-DISPLACEMENT` |
| Output | `.frd` converted to `.vtu` / `.pvd`; CSV probe and frame histories |

---

## Process Sequence

```
  Stage A  MESH        structured hexahedral plate, fine 2 mm band around the
                       track; deposit layer pre-meshed in the stadium footprint
  Stage B  REMOVE      step 1: all deposit elements switched off,
                       plate clamped, at 20 degC
  Stage C  PLUNGE      tool dwells 4 s at the start (x = 45 mm), heating the
                       plate; the first disc of deposit is switched on
  Stage D  DEPOSIT     55 steps of 2 mm (1 s each): tool advances, the leading
                       crescent of deposit is switched on at 380 degC,
                       shoulder heat moves with the tool
  Stage E  COOL        tool lifted, 12 cooling steps of growing length (600 s)
  Stage F  UNCLAMP     clamps released, plate supported on three points:
                       residual stress redistributes, distortion appears
  Stage G  OUTPUT      one result frame per step (69 frames);
                       converter writes part + tool animation
```

---

## Geometry

```
      z
      |              feedstock 12.7 x 12.7 mm (rotates with the shoulder)
      |                     | |
      |                 +---+ +---+   hollow shoulder, D = 38.1 mm
      |                 |  shoulder |  ---> deposition direction (+x)
      |   ,-------------+-----------+--.
      |  (   deposited layer, 2 mm        )   z = 10 to 12 mm
      |   `----------------------------'
      |  +--------------------------------------------------+  z = 10 mm
      |  |              AA6061 substrate                     |
      |  +--------------------------------------------------+  z = 0
      +--------------------------------------------------------- x
      0        x = 45 mm  (start)          x = 155 mm (end)        200 mm

  Plate:       200 mm (x) x 100 mm (y) x 10 mm (z)
  Tool path:   y = 0, from x = 45 to x = 155 mm (110 mm)
  Deposit:     stadium = swept shoulder footprint, 38.1 mm wide,
               148.1 mm long including the rounded ends
  Clamps:      top surface strips |y| >= 40 mm
  Symmetry:    y = 0 plane (half model solved, mirrored for display)
```

A deposit element belongs to the layer if its centre lies inside the
stadium

```
  (x - clip(x, x_start, x_end))^2 + y^2 <= R^2,      R = 19.05 mm
```

and it is switched on in the first step in which the shoulder disc has
swept over its centre.

---

## Physics

### Heat conduction (transient)

```
  rho c(T) dT/dt = div( k(T) grad T )          in substrate and active deposit
```

### Frictional heat from the shoulder

The shoulder heat is applied as a surface flux on the top faces of the
deposit under the shoulder. For a rotating disc with a uniform contact
shear stress, the local heat generation grows with the sliding speed,
which is proportional to the radius:

```
  q(r) = 3 P r / (2 pi R^3),       0 <= r <= R

  integral over the disc:   int q 2 pi r dr = P
```

On the discrete mesh the face fluxes are rescaled so the applied total
is exactly P (P/2 in the half model).

### Material entering the layer

New material leaves the tool hot. Nodes that belong only to deposit
elements are given the initial temperature `T_DEP` and keep it while
inactive; when their elements are switched on they enter the analysis
at that temperature. Thermal strain is measured from the initial
temperature, so the deposit is **stress-free at 380 degC** and shrinks
as it cools, which is what generates the residual stress.

### Heat losses

```
  Exposed faces (air):         q = h_air   (T - T_amb),    h_air   = 15  W/m^2K
  Substrate bottom (anvil):    q = h_anvil (T - T_amb),    h_anvil = 300 W/m^2K
  Faces under the shoulder:    no loss (covered by the tool, receive the flux)
  Symmetry plane y = 0:        adiabatic
```

The set of exposed faces is recomputed every step from the currently
active elements, so a face stops losing heat as soon as new material
covers it. Radiation is lumped into `h_air` (aluminium has a low
emissivity).

### Elastoplasticity

```
  div(sigma) = 0                                      quasi-static
  eps = eps_el + eps_pl + eps_th,   eps_th = alpha(T) (T - T_0)
  von Mises yield,  sigma_y(T, eps_p),  isotropic hardening
```

Small-strain kinematics. Elements added by `*MODEL CHANGE, ADD` are
activated strain-free on the deformed substrate (verified in a
single-element test).

### Boundary conditions

| Node set | Region | Constraint | Steps |
|---|---|---|---|
| `NBOT` | Substrate bottom face | u_z = 0 (anvil) | Deposition, cooling |
| `NCLAMP` | Top-surface strips, \|y\| >= 40 mm | u_x = u_y = u_z = 0 | Deposition, cooling |
| `NSYM` | Symmetry plane y = 0 | u_y = 0 | All |
| Corner (0, 0, 0) | Bottom corner | u_x = u_z = 0 | Unclamped |
| Corner (200, 0, 0) | Bottom corner | u_z = 0 | Unclamped |
| Corner (0, 50, 0) | Bottom corner | u_z = 0 | Unclamped |

The three-point support in the unclamped step only removes rigid-body
motion, so the plate is free to bow.

---

## Material Properties

AA6061-T6, temperature dependent (representative literature values;
linear interpolation between points).

| T (degC) | E (MPa) | Yield (MPa) |
|---|---|---|
| 20 | 68,900 | 276 |
| 100 | 66,200 | 262 |
| 150 | | 240 |
| 200 | 59,200 | 200 |
| 250 | | 140 |
| 300 | 47,800 | 90 |
| 350 | | 50 |
| 400 | 31,700 | 30 |
| 450 | | 18 |
| 500 | 15,000 | 12 |
| 580 | 5,000 | 5 |

| T (degC) | k (W/mK) | c_p (J/kgK) | alpha (1/K) |
|---|---|---|---|
| 20 | 167 | 896 | 23.4e-6 |
| 200 | 180 | 960 | 24.9e-6 |
| 400 | 195 | 1050 | 27.5e-6 |
| 580 | 205 | 1120 | 29.0e-6 |

| Property | Value | Unit |
|---|---|---|
| Density | 2700 | kg/m^3 |
| Poisson's ratio | 0.33 | n/a |
| Hardening (each temperature) | (sigma_y, 0), (1.15 sigma_y, 0.2), (1.2 sigma_y, 2.0) | MPa, n/a |
| Expansion reference temperature | 20 | degC |

The T6 values do not account for the overageing that happens during
AFSD, so the post-process strength is overestimated. Replace the tables
with your own data for quantitative work.

---

## Model Parameters

### Substrate and tool

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Plate length / width / thickness | `PLATE_L`, `PLATE_W`, `PLATE_T` | 200, 100, 10 | mm |
| Substrate element heights (bottom to top) | `NZ_SUB` | 4, 3, 2, 1 | mm |
| Shoulder radius | `R_SH` | 19.05 | mm |
| Feedstock side | `FEED_A` | 12.7 | mm |
| Spindle speed (animation only) | `RPM` | 300 | rpm |

### Process

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Number of layers | `N_LAYERS` | 1 | n/a |
| Layer thickness | `LAYER_H` | 2 | mm |
| Path start / end | `X_START`, `X_END` | 45, 155 | mm |
| Traverse speed | `V_TRAVEL` | 2 | mm/s |
| Tool advance per step (frame) | `DX_STEP` | 2 | mm |
| Plunge dwell | `T_PLUNGE` | 4 | s |
| Interlayer dwell (multi-layer) | `T_INTERLAYER` | 20 | s |
| Shoulder heat input | `P_TOOL` | 2400 | W |
| Deposit entry temperature | `T_DEP` | 380 | degC |
| Ambient temperature | `T_AMB` | 20 | degC |

### Cooling, losses and clamping

| Parameter | Symbol | Default | Unit |
|---|---|---|---|
| Cooling step lengths | `COOL_TIMES` | 2, 3, 5, 8, 12, 20, 30, 45, 60, 90, 125, 200 | s |
| Unclamp after cooling | `UNCLAMP` | True | n/a |
| Air convection | `H_AIR` | 15 | W/m^2K |
| Anvil conductance | `H_ANVIL` | 300 | W/m^2K |
| Clamp strips | `CLAMP_Y` | \|y\| >= 40 | mm |

### Mesh and analysis

| Parameter | Symbol | Default |
|---|---|---|
| Analysis type | `ANALYSIS` | `coupled` (or `thermal`) |
| Half model with symmetry | `SYMMETRY` | True |
| In-plane size near the track | `H_FINE` | 2 mm |
| In-plane size away from the track | `H_COARSE_X`, `H_COARSE_Y` | 8, 6.5 mm |

All parameters are set at the top of `afsd_ccx.py`. Setting
`N_LAYERS > 1` builds a multi-layer wall with the tool direction
reversed every layer and an interlayer dwell.

---

## Derived Design Values (computed from the default inputs)

These follow directly from the parameters above. They are **not**
simulation results.

| Quantity | Value |
|---|---|
| Deposit cross-section (38.1 x 2 mm) | 76.2 mm^2 |
| Feedstock cross-section | 161.3 mm^2 |
| Feedstock consumed per mm of travel | 0.47 mm |
| Feedstock feed speed | 0.94 mm/s |
| Deposition rate | 152 mm^3/s (9.1 cm^3/min) |
| Deposited footprint area | 5,331 mm^2 |
| Heat input per unit length (P / v) | 1,200 J/mm |
| Peak shoulder flux at r = R, 3P / (2 pi R^2) | 3.2 W/mm^2 |
| Shoulder rim speed at 300 rpm | 0.60 m/s |
| Deposition time (plunge + 55 steps) | 59 s |
| Total simulated time | 660 s |

---

## Mesh

| Quantity | Value (half model) |
|---|---|
| Nodes | 7,830 |
| Substrate elements | 5,248 C3D8 |
| Deposit elements | 690 C3D8 |
| Total elements | 5,938 |
| Element size near the track | 2 x 2 mm in plane |
| Steps (frames) | 69 coupled / 68 thermal |

Treat the console line printed by `afsd_ccx.py` as authoritative; it
changes with any geometry, mesh or process parameter.

---

## Results

### Thermal history (thermal-only run, default parameters)

Peak temperatures at probe points in the middle of the track
(x = 100 mm):

| Probe | Location | Peak T (degC) |
|---|---|---|
| P5 | Deposit top surface (z = 12 mm) | 401 |
| P1 | Deposit / substrate interface (z = 10 mm) | 380 |
| P2 | Substrate, 3 mm below the interface | 358 |
| P3 | Substrate bottom (anvil side) | 337 |
| P4 | Substrate top, beside the track (y = 30 mm) | 247 |
| | Highest temperature anywhere in the model | 435 |

The peak temperature rises along the track as the plate accumulates
heat, so the end of the deposit runs hotter than the start. These
values lie in the range reported for AA6061 AFSD (roughly 350 to
480 degC). `P_TOOL` should be calibrated against thermocouple data
(mechanical power is approximately spindle torque x angular speed).



### Coupled run

<!-- Fill in from afsd_frames.csv and ParaView after your run -->

| Quantity | Value |
|---|---|
| Peak von Mises stress during deposition | from `afsd_frames.csv` |
| Peak equivalent plastic strain | from ParaView (`PEEQ`) |
| Residual von Mises stress after unclamping | from the last frame |
| Maximum distortion after unclamping | from the last frame (`U`) |

### What to look for

| Field | Expected behaviour |
|---|---|
| `T` | Hot zone travels with the tool; trail of cooling deposit behind it; plate heats up along the track |
| `S_Mises` | Compressive yielding ahead of and under the tool (hot material restrained by cold surroundings); tension in the deposit and near-interface substrate after cooling |
| `S` (xx) | Longitudinal tension along the cooled deposit, balanced by compression in the rest of the plate |
| `PEEQ` | Plastic strain concentrated in the deposit and the substrate just below it |
| `U` | After unclamping the plate bows along the track (use `WARP_SCALE` to exaggerate) |
| Tool | Shoulder and feedstock move along the path; the feedstock shortens; the tool lifts away during cooling |

### Known limitations

* **No material flow.** The deposit geometry is prescribed, not
  predicted. Use a CFD (OpenFOAM), SPH or MPM model for material flow,
  layer shape and flow-related defects.
* **Heat input is prescribed.** Friction and plastic work are not
  computed; `P_TOOL` and `T_DEP` must be calibrated.
* **No tool downforce.** The forging pressure of the shoulder is not
  applied.
* **Slow plasticity convergence.** Hot aluminium is almost perfectly
  plastic, so CalculiX cuts the increment several times in many steps.
  The unclamping step alone can need around 100 increments.
* **No microstructure.** Precipitate evolution and overageing are not
  modelled; the thermal history is the input for that.

---

## Output Files

| File | Description |
|---|---|
| `afsd.inp` | CalculiX input deck written by `afsd_ccx.py` |
| `afsd_model.npz` | Mesh and element activation step (used by the converter) |
| `afsd_schedule.json` | Step list: end time, tool position, layer, active elements |
| `afsd.frd` | Field results, one frame per step |
| `afsd.sta`, `.cvg` | Solver status and convergence history |
| `afsd_ccx.log` | Solver console output |
| `afsd_frames.csv` | Time, stage, tool position, active elements, peak T and peak stress per frame |
| `afsd_probes.csv` / `.png` | Temperature history at the five probe points |
| `pvd/afsd_part.pvd` | ParaView collection: substrate + deposit |
| `pvd/afsd_part_0000.vtu` ... | Part frames (only active elements, mirrored to the full plate) |
| `pvd/afsd_tool.pvd` | ParaView collection: shoulder + feedstock |
| `pvd/afsd_tool_0000.vtu` ... | Tool frames |

For the thermal-only run every name starts with `afsd_th` instead of
`afsd`.

### VTU data fields (part)

| Field | Type | Description |
|---|---|---|
| `T` | point | Temperature (degC) |
| `U` | point | Displacement (mm), coupled run only |
| `S` | point | Stress tensor, 6 components (MPa), coupled run only |
| `S_Mises` | point | Von Mises stress (MPa), coupled run only |
| `PEEQ` | point | Equivalent plastic strain, coupled run only |
| `Region` | cell | 0 = substrate, 1..N = deposit layer number |
| `ActivationStep` | cell | Step in which the element was deposited |

### VTU data fields (tool)

| Field | Type | Description |
|---|---|---|
| `ToolPart` | cell | 0 = shoulder, 1 = feedstock |
| `Stripe` | cell | Alternating shoulder sectors (makes rotation visible), 2 = feedstock |

### CSV columns

```
afsd_frames.csv : frame, time_s, stage, tool_x_mm, active_elements, Tmax_C, Mises_max_MPa
afsd_probes.csv : time_s, P1_interface_centre, P2_substrate_3mm_deep,
                  P3_substrate_bottom, P4_beside_track_y30, P5_deposit_top
```

---

## Repository Structure

```
.
├── afsd_ccx.py           # Generates the mesh, activation schedule and CalculiX deck
├── run_afsd.sh           # Runs the full pipeline: generate, solve, convert
├── frd2pvd.py            # .frd -> part + tool .vtu/.pvd animation, probe curves
├── afsd_view.py          # ParaView script: part + tool, colour map, camera, movie
├── afsd_probes.png
└── README.md             # This file
```

---

## How to Run

### Requirements
* CalculiX `ccx` 2.21 or later: http://www.calculix.de
* Python 3 with numpy
* matplotlib (optional, for the probe plot)
* ParaView 5.10 or later: https://www.paraview.org

No ccx2paraview or vtk is needed: `frd2pvd.py` writes the `.vtu` files
directly.

On Ubuntu or WSL:

```bash
sudo apt install calculix-ccx python3-pip
pip install numpy matplotlib
```

If pip reports an externally managed environment, use
`sudo apt install python3-numpy python3-matplotlib` instead.

### Step 1: Run the full pipeline

```bash
bash run_afsd.sh thermal     # temperature only, a few minutes
bash run_afsd.sh             # coupled temperature + stress + distortion
```

Set `NP` in `run_afsd.sh` to the number of CPU cores (`nproc`).
Follow solver progress from a second terminal with `tail -f afsd.sta`;
the run is complete when the total time reaches 660 s (659 s for the
thermal run).

### Step 2: Or run the stages separately

```bash
python3 afsd_ccx.py            # writes afsd.inp  (python3 afsd_ccx.py thermal -> afsd_th.inp)
ccx afsd                       # solve
python3 frd2pvd.py afsd        # pvd/ frames, probe and frame CSVs
```

The job name is given to `ccx` without the `.inp` extension. If the
solve has finished, only the last command needs to be repeated after a
change to the converter.

---

## Visualisation in ParaView

### Option A: Script (recommended)

```
Edit > Reset Session
View > Python Shell > Run Script > afsd_view.py
Press Play
```

The script loads both series, warps the part by `U` (coupled run),
colours it by `T`, shows the shoulder in grey and the feedstock in
orange, adds a time label and sets an isometric view. At the top of the
script:

| Setting | Meaning |
|---|---|
| `PVD_DIR` | Folder with the `.pvd` files; on Windows use the full path with `/` |
| `FIELD` | `T`, `S_Mises`, `PEEQ` or `U` |
| `RANGES` | Colour range per field |
| `WARP_SCALE` | Displacement magnification (0 = off) |
| `SAVE_MOVIE` | `True` writes PNG frames to `./movie/` |

Combine saved frames into a video with

```bash
ffmpeg -framerate 12 -i movie/afsd.%04d.png -pix_fmt yuv420p afsd.mp4
```

### Option B: By hand

```
File > Open > pvd/afsd_part.pvd and pvd/afsd_tool.pvd (Ctrl-click), Apply
Select afsd_part.pvd: colour by T or S_Mises
Optional: Filters > Warp By Vector, Vectors = U, Scale Factor = 10
Select afsd_tool.pvd: colour by ToolPart or Solid Color
Press Play
```

### Option C: See under the shoulder

```
Select afsd_tool.pvd, set Opacity = 0.4
```

The shoulder covers the hottest region; a semi-transparent tool shows
the temperature field underneath while keeping the tool visible.

---

## Key Numerical Insights

### Element activation instead of remeshing
The whole deposit is meshed in advance. `*MODEL CHANGE, TYPE=ELEMENT,
REMOVE` switches it off in step 1 and `ADD` switches on the crescent of
elements swept by the shoulder in each step. No remeshing or mapping is
needed, and the same mesh carries the results from deposition to
unclamping.

### Hot, stress-free activation
Nodes of inactive elements keep their initial temperature until they
are activated, so setting the deposit-only nodes to `T_DEP` makes the
new material enter hot. Because added elements are activated
strain-free and thermal strain is referenced to the initial
temperature, the deposit carries no stress at the moment it is placed
and builds up stress only as it cools.

### Moving surfaces
Free surfaces change every step as the layer grows. The generator
rebuilds the `*FILM` list from the current active elements in every
step, and leaves out faces under the shoulder, which receive the flux
instead.

### Exact heat input on a coarse mesh
The flux q(r) is evaluated at each face centre and rescaled so the
faces under the shoulder receive exactly `P_TOOL`, independent of mesh
size or how the disc cuts the element grid.

### One frame per step
`*NODE FILE, FREQUENCY=10000` writes results only at the end of each
step, so the animation has exactly one frame per tool position no
matter how many increments the solver needs.

### Half model
The track, heat input and clamps are symmetric about y = 0, so only
half the plate is solved (about 4x faster per iteration). The converter
mirrors coordinates, connectivity and fields (u_y, S_xy and S_yz change
sign) to display the full plate.

### Numpy-only converter
CalculiX writes only the active nodes in each frame. `frd2pvd.py` reads
the fixed-width `.frd` records, keeps the elements active at that time
and writes compressed binary `.vtu` files, so the layer visibly grows
in ParaView without any thresholding.

### Tool animation
The tool series is built from the step schedule: the shoulder (hollow,
square bore) and feedstock follow the tool position, the feedstock
shortens by 0.47 mm per mm of travel, the shoulder turns 24 degrees per
frame (the true 300 rpm would alias) and the tool lifts away during
cooling.

---

## Common Errors and Fixes

| Error / symptom | Cause | Fix |
|---|---|---|
| `ccx: command not found` | CalculiX not installed or not on PATH | `sudo apt install calculix-ccx` |
| `run_afsd.sh: No such file or directory` | Zip extracted into a folder of the same name | `cd` into the inner `afsd_calculix` folder |
| `$'\r': command not found` | Windows line endings | `sed -i 's/\r$//' run_afsd.sh` |
| `error: externally-managed-environment` | System Python protected | `sudo apt install python3-numpy python3-matplotlib` |
| `ValueError: could not convert string to float: ''` in `read_frd` | Old converter read 4 displacement components | Use the updated `frd2pvd.py`, then rerun only `python3 frd2pvd.py afsd` |
| Coupled run very slow, many cut-backs in `afsd.sta` | Near-perfect plasticity of hot aluminium | Expected; use more threads, run `thermal` first, or add more hardening in `*PLASTIC` |
| No `pvd` folder | Conversion not run or failed | `python3 frd2pvd.py afsd` and read its output |
| ParaView cannot find the files | Script run from another working directory | Set `PVD_DIR` to the full path, e.g. `C:/Users/<you>/.../pvd` |
| Only temperature available in ParaView | Thermal-only result (`afsd_th`) loaded | Run `bash run_afsd.sh` for stress, then reload |
| Hot spot hidden by the tool | Shoulder covers the heated area | Set tool Opacity to 0.4 |
| Plate looks undeformed | Displacements are small | Increase `WARP_SCALE` |
| Overlapping duplicate objects | Script run several times | Edit > Reset Session, then run once |
| SyntaxError when pasting the script | Python Shell cannot take pasted indented blocks | Use **Run Script**, not paste |
| Temperatures too low / too high | Heat input not calibrated | Adjust `P_TOOL` (and `T_DEP`) against thermocouples |

---

## Extending the Model

| Extension | What to change |
|---|---|
| Multi-layer wall | `N_LAYERS` (direction reverses each layer), `T_INTERLAYER` |
| Process window study | Sweep `V_TRAVEL`, `P_TOOL` (from torque x speed), `T_DEP` |
| Different alloys | Replace the property tables (AA2024, AA7075, Cu, Ti-6Al-4V, steels) |
| Dissimilar deposition | Separate material for substrate and deposit |
| Tool downforce | Add a `*DLOAD` pressure under the shoulder in each deposition step |
| Calibrated heat input | Fit `P_TOOL` and `T_DEP` to measured thermocouple curves (`afsd_probes.csv`) |
| Faster thermo-mechanics | Sequential coupling: thermal run first, then a static run driven by the temperature field |
| Microstructure | Feed the probe histories into a precipitation / grain-growth model |
| Material flow | Couple with an OpenFOAM, SPH or MPM model of the region under the tool |
| Non-straight paths | Replace the straight path with a list of tool positions |
| Parametric study / ML dataset | Loop the generator over process parameters |

---

## Citation

Mishra, A. (2026). *Additive Friction Stir Deposition (AFSD) of AA6061
Using CalculiX* [Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.23558007

```bibtex
@software{mishra2026afsd,
  author       = {Mishra, Akshansh},
  title        = {Additive Friction Stir Deposition (AFSD) of AA6061
                   Using CalculiX},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23558007},
  url          = {https://doi.org/10.5281/zenodo.23558007}
}
```

---

## License

<p>
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
<img alt="Creative Commons Licence" style="border-width:0"
  src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"/>
</a>
<br/>
This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
Creative Commons Attribution NonCommercial 4.0 International License</a>.
</p>

You are free to:
* **Share**: copy and redistribute in any medium or format
* **Adapt**: remix, transform, and build upon the material

Under the following terms:
* **Attribution**: give appropriate credit and link to this repository
* **NonCommercial**: not for commercial use without permission

Copyright 2026 Akshansh Mishra.
