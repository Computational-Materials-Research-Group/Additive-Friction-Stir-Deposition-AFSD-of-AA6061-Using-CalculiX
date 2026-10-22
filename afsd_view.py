# ParaView script for the AFSD animation
#   GUI:      Edit > Reset Session, then View > Python Shell > Run Script > afsd_view.py
#   batch:    pvpython afsd_view.py        (set SAVE_MOVIE = True to write frames)
#
# Loads pvd/afsd_part.pvd (substrate + growing deposit) and pvd/afsd_tool.pvd
# (rotating hollow shoulder + feedstock), colours the part by temperature or
# stress, and sets an isometric camera. Press Play.
from paraview.simple import *
import os

PVD_DIR = "pvd"              # folder holding the .pvd files (relative to this script's working dir)
JOB = "afsd"                 # "afsd" (coupled) or "afsd_th" (thermal only)
FIELD = "T"                  # "T" (degC), "S_Mises" (MPa), "PEEQ", or "U"
RANGES = {"T": (20, 450), "S_Mises": (0, 250), "PEEQ": (0, 0.05), "U": (0, 0.3)}
WARP_SCALE = 10.0            # displacement magnification (0 = no warp)
SHOW_EDGES = False
SAVE_MOVIE = False           # True -> writes PNG frames to ./movie/
MOVIE_RES = [1600, 900]

base = os.path.join(os.getcwd(), PVD_DIR)
if not os.path.exists(os.path.join(base, f"{JOB}_part.pvd")) and os.path.exists(os.path.join(base, "afsd_th_part.pvd")):
    JOB = "afsd_th"            # fall back to the thermal-only animation
part = PVDReader(registrationName="AFSD part", FileName=os.path.join(base, f"{JOB}_part.pvd"))
tool = PVDReader(registrationName="AFSD tool", FileName=os.path.join(base, f"{JOB}_tool.pvd"))

view = GetActiveViewOrCreate("RenderView")
view.OrientationAxesVisibility = 1
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass

src = part
part.UpdatePipeline()
has_u = "U" in part.PointData.keys()
if not has_u:
    FIELD = "T"
if has_u and WARP_SCALE > 0:
    src = WarpByVector(registrationName="Warp", Input=part)
    src.Vectors = ["POINTS", "U"]
    src.ScaleFactor = WARP_SCALE

disp = Show(src, view)
disp.Representation = "Surface With Edges" if SHOW_EDGES else "Surface"
ColorBy(disp, ("POINTS", FIELD, "Magnitude") if FIELD == "U" else ("POINTS", FIELD))
lut = GetColorTransferFunction(FIELD)
lut.ApplyPreset("Inferno (matplotlib)" if FIELD == "T" else "Turbo", True)
lut.RescaleTransferFunction(*RANGES.get(FIELD, (0, 1)))
lut.AutomaticRescaleRangeMode = "Never"
disp.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = {"T": "Temperature (°C)", "S_Mises": "von Mises (MPa)",
             "PEEQ": "Equivalent plastic strain", "U": "Displacement (mm)"}.get(FIELD, FIELD)
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0.1, 0.1, 0.1]

tdisp = Show(tool, view)
ColorBy(tdisp, ("CELLS", "Stripe"))
tlut = GetColorTransferFunction("Stripe")
tlut.InterpretValuesAsCategories = 1
tlut.Annotations = ["0", "shoulder", "1", "shoulder", "2", "feedstock"]
tlut.IndexedColors = [0.62, 0.64, 0.68, 0.45, 0.47, 0.52, 0.85, 0.55, 0.25]
tdisp.SetScalarBarVisibility(view, False)
tdisp.Specular = 0.4

txt = AnnotateTimeFilter(registrationName="Time", Input=part)
try:
    txt.Format = "t = {time:.1f} s"
except Exception:
    txt.Format = "t = %.1f s"
tshow = Show(txt, view)
tshow.FontSize = 14
tshow.Color = [0.1, 0.1, 0.1]
tshow.WindowLocation = "Upper Left Corner"

view.CameraPosition = [-60, -210, 170]
view.CameraFocalPoint = [100, 0, 8]
view.CameraViewUp = [0, 0, 1]
view.CameraParallelScale = 90

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.PlayMode = "Snap To TimeSteps"
Render()

if SAVE_MOVIE:
    os.makedirs("movie", exist_ok=True)
    SaveAnimation(os.path.join("movie", "afsd.png"), view, ImageResolution=MOVIE_RES)
    print("PNG frames written to ./movie (combine with ffmpeg, see README)")
