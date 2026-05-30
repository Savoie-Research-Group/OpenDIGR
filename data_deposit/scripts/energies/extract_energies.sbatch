#!/bin/bash
#SBATCH --job-name=extract_energies
#SBATCH --output=extract_energies_%j.out
#SBATCH --error=extract_energies_%j.err
#SBATCH -A bsavoie
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00

set -euo pipefail

# --- CONFIGURE FOR YOUR SYSTEM ---
# SCRIPT  : path to extract_energies.py
# KEYS    : newline-separated reaction keys to process
# ZIP_DIR : root of the per-reaction .zip archives
# OUTPUT  : output CSV (job-id suffix optional)
SCRIPT=${SCRIPT:-/path/to/extract_energies.py}
KEYS=${KEYS:-/path/to/reaction_keys.txt}
ZIP_DIR=${ZIP_DIR:-/path/to/zip_archives}
OUTPUT=${OUTPUT:-energies_${SLURM_JOB_ID}.csv}

echo "Job:       $SLURM_JOB_ID on $(hostname)"
echo "CPUs:      $SLURM_CPUS_PER_TASK"
echo "Keys:      $(wc -l < "$KEYS") rows"
echo "Output:    $OUTPUT"
echo "Start:     $(date)"
echo

python3 "$SCRIPT" \
    --keys "$KEYS" \
    --zip-dir "$ZIP_DIR" \
    -o "$OUTPUT" \
    -j "$SLURM_CPUS_PER_TASK"

echo
echo "Rows written: $(($(wc -l < "$OUTPUT") - 1))"
echo "End:       $(date)"
