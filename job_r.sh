#!/bin/bash
# =============================================================================
# job_r.sh  -  runs any R script of the project on a SLURM cluster
# Usage: sbatch job_r.sh <script.R> [arguments...]
#   sbatch job_r.sh 08_sdm_model.R epilosa tgb
#   sbatch job_r.sh 09_variable_importance.R epilosa
# Edit the module name, memory and time to match your cluster.
# =============================================================================
#SBATCH --job-name=syn_r
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
module load R/R-4.4.0
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
echo "Job $SLURM_JOB_ID | $* | start: $(date)"
Rscript "$@"
STATUS=$?
echo "End: $(date) | exit code: $STATUS"
exit $STATUS
