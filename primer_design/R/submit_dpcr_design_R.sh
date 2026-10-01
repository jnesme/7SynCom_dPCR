#!/bin/bash
### General options
### -- specify queue --
#BSUB -q hpc
### -- set the job Name --
#BSUB -J dpcr_design_R
### -- ask for number of cores --
#BSUB -n 8
### -- all cores on one host, 4GB per core (8 x 4GB = 32GB) --
#BSUB -R "span[hosts=1] rusage[mem=4GB]"
### -- specify that we want the job to get killed if it exceeds 4.5GB per core/slot --
#BSUB -M 4500MB
### -- set walltime limit: hh:mm --
#BSUB -W 04:00
### -- set the email address --
#BSUB -u josne@dtu.dk
### -- send notification at completion --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o logs/dpcr_design_R_%J.out
#BSUB -e logs/dpcr_design_R_%J.err

#==========================================================================
# EDIT THESE BEFORE SUBMITTING
#==========================================================================
# Steps to run (space-separated). Each step skips work whose output already exists.
STEPS="00 01 01a 02 03 04 05 06 07"
#==========================================================================

# Pipeline directory (submit from here: bsub < submit_dpcr_design_R.sh)
PIPE_DIR="/work3/josne/Projects/7SynCom_Kolter/primer_design/R"
CONDA_ENV="dpcr-design-r"

cd "${PIPE_DIR}" || exit 1
mkdir -p logs work results

# Load environment
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}" || { echo "ERROR: conda env ${CONDA_ENV} not found"; exit 1; }

# Print job information
echo "=========================================="
echo "dPCR multiplex design (R version) - 7 SynCom genomes"
echo "Job started on $(date)"
echo "Job ID: $LSB_JOBID"
echo "Running on node: $(hostname)"
echo "Cores: ${LSB_DJOB_NUMPROC:-1}"
echo "Steps: ${STEPS}"
echo "=========================================="

export THREADS="${LSB_DJOB_NUMPROC:-1}"
EXIT_CODE=0
for STEP in ${STEPS}; do
    SCRIPT=$(ls scripts/${STEP}_*.R 2>/dev/null | head -1)
    if [ -z "${SCRIPT}" ]; then
        echo "WARNING: no script for step ${STEP}, skipping"
        continue
    fi
    echo "---- $(date +%T) running ${SCRIPT}"
    Rscript "${SCRIPT}"
    EXIT_CODE=$?
    if [ ${EXIT_CODE} -ne 0 ]; then
        echo "ERROR: ${SCRIPT} failed (exit ${EXIT_CODE})"
        break
    fi
done

# Print completion information
echo "=========================================="
echo "Job finished on $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "=========================================="

exit ${EXIT_CODE}
