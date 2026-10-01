#!/bin/bash
# ADD-ON: stricter within-assay dimer limit. Leaves the standard scripts, results and
# figures untouched; everything is written under primer_design/strict/.
#   bash run_strict_addon.sh        (about 5 min on one core; needs steps 00-04 done)
# The standard steps 06 (verification) and 07 (figures) are reused unmodified through
# symlinks: they locate work/, results/ and figures/ relative to their own path.
set -euo pipefail
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PIPE_DIR}"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate dpcr-design
export THREADS="${LSB_DJOB_NUMPROC:-1}"

mkdir -p strict/scripts strict/work strict/results strict/figures
for f in 06_verify_report.py 07_figures.py insilico_pcr.py; do
    ln -sfn "../../scripts/${f}" "strict/scripts/${f}"
done
# share the standard intermediate files read-only; files written by 06 land in strict/work
for f in work/*; do
    b=$(basename "$f")
    case "$b" in final.*|pool*|ko_*) continue ;; esac
    [ -e "strict/work/$b" ] || ln -s "../../work/$b" "strict/work/$b"
done

python scripts/08_strict_dimer_addon.py
python strict/scripts/06_verify_report.py
python strict/scripts/07_figures.py
rm -f strict/figures/fig1_workflow.*   # the workflow figure is the same as the standard one
