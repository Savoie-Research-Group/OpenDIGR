#!/bin/bash
#SBATCH --job-name=patch_xtb_keep
#SBATCH --output=patch_logs/chunk_%A_%a.out
#SBATCH --error=patch_logs/chunk_%A_%a.err
#SBATCH -A bsavoie
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --array=1-200%20

set -euo pipefail

# --- CONFIGURE FOR YOUR SYSTEM ---
# WORKDIR  : directory containing patch_xtb_keep_worker.py and chunk lists.
# SRC_DIR  : source archives with the xTB-scan content to copy over.
# DST_DIR  : pre-existing slim archives to patch (must already exist).
WORKDIR=${WORKDIR:-/path/to/workdir}
SCRIPT="$WORKDIR/patch_xtb_keep_worker.py"
SRC_DIR=${SRC_DIR:-/path/to/source/archives}
DST_DIR=${DST_DIR:-/path/to/slim/output}
CHUNKDIR="$WORKDIR/trim_chunks"   # reuse same chunks from the trim job

IDX=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
KEYS="$CHUNKDIR/chunk_${IDX}.txt"

echo "Array task: $SLURM_ARRAY_TASK_ID  Job: $SLURM_JOB_ID  Node: $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Keys file:  $KEYS  ($(wc -l < "$KEYS") zips)"
echo "Src (RO):   $SRC_DIR"
echo "Dst (RW):   $DST_DIR"
echo "Start:      $(date)"

# Inode-quota guard (no NEW inodes are created by this patch -- we only
# rewrite existing slim zips in place via tmp + rename -- but keep the
# guard for symmetry with the trim job).
INODE_GUARD=9700000
FILES_USED=$(lfs quota -u "$USER" "${LUSTRE_MOUNT:-/scratch}" 2>/dev/null | tail -1 | awk '{print $5}')
FILES_LIMIT=$(lfs quota -u "$USER" "${LUSTRE_MOUNT:-/scratch}" 2>/dev/null | tail -1 | awk '{print $7}')
echo "Quota:      ${FILES_USED} / ${FILES_LIMIT} inodes  (guard at ${INODE_GUARD})"
if [ -n "$FILES_USED" ] && [ "$FILES_USED" -ge "$INODE_GUARD" ]; then
    echo "ABORT: inode usage $FILES_USED has reached guard threshold $INODE_GUARD."
    exit 0
fi
echo

python3 "$SCRIPT" \
    --keys "$KEYS" \
    --src-dir "$SRC_DIR" \
    --dst-dir "$DST_DIR" \
    -j "$SLURM_CPUS_PER_TASK"

echo
echo "End:        $(date)"
