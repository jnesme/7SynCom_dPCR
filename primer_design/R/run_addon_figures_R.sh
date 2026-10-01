#!/bin/bash
# R figures for the low-dimer add-on set (Figures 2-4).
# The add-on selection itself exists in Python only (../scripts/08_strict_dimer_addon.py,
# run with ../run_strict_addon.sh first). This script only DRAWS its results in R:
# it copies the add-on results next to symlinks of the R intermediate files and runs the
# unmodified R figure script (scripts/07_figures.R) from that folder.
#   bash run_addon_figures_R.sh        (2-3 min on one core)
# Output: primer_design/R/strict/figures/*_R.png / .pdf
set -euo pipefail
R_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${R_DIR}/../strict/results"          # results of the Python add-on
cd "${R_DIR}"
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate dpcr-design-r
export THREADS="${LSB_DJOB_NUMPROC:-1}"

[ -f "${SRC}/combination_scores.tsv" ] || { echo "ERROR: run ../run_strict_addon.sh first"; exit 1; }
mkdir -p strict/work strict/results strict/figures
ln -sfn ../scripts strict/scripts
for f in work/*; do
    b=$(basename "$f")
    case "$b" in final.*|pool*|fig3.*|combination_scores.tsv) continue ;; esac
    [ -e "strict/work/$b" ] || ln -s "../../work/$b" "strict/work/$b"
done
cp "${SRC}/final_multiplex.tsv" "${SRC}/candidates_pool.tsv" "${SRC}/candidates_pool_wide.tsv" \
   "${SRC}/specificity_report.tsv" strict/results/
cp "${SRC}/combination_scores.tsv" strict/work/
# the Python table has an empty first header cell; the R script expects "oligo"
sed '1s/^/oligo/' "${SRC}/dimer_matrix.tsv" > strict/results/dimer_matrix.tsv

cd strict
Rscript scripts/07_figures.R
rm -f figures/fig1_workflow_R.*      # the workflow figure does not apply to the add-on
