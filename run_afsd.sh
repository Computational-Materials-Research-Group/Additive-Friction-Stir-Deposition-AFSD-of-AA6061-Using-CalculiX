#!/bin/bash
# ============================================================
#  AFSD (additive friction stir deposition) with CalculiX -> ParaView
#
#  One-time install (Ubuntu / WSL):
#     sudo apt install calculix-ccx python3-pip
#     pip install numpy matplotlib
#
#  Run:   bash run_afsd.sh            coupled thermo-mechanical (temperature + stress)
#         bash run_afsd.sh thermal    temperature only (a few minutes)
# ============================================================
set -e
MODE=${1:-coupled}
NP=14                                   # CPU threads
export OMP_NUM_THREADS=$NP
export CCX_NPROC_STIFFNESS=$NP
export CCX_NPROC_EQUATION_SOLVER=$NP
export CCX_NPROC_RESULTS=$NP

if [ "$MODE" = "thermal" ]; then JOB=afsd_th; python3 afsd_ccx.py thermal
else JOB=afsd; python3 afsd_ccx.py; fi

echo ">>> Solving $JOB with CalculiX ($NP threads). Progress: tail -f ${JOB}.sta"
ccx $JOB > ${JOB}_ccx.log 2>&1
tail -n 1 $JOB.sta

echo ">>> Converting to ParaView"
python3 frd2pvd.py $JOB

echo "Done. In ParaView: View > Python Shell > Run Script > afsd_view.py, then press Play."
