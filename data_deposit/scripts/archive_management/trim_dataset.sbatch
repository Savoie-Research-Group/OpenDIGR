#!/bin/bash
#SBATCH --job-name=trim_dataset
#SBATCH --output=trim_logs/chunk_%A_%a.out
#SBATCH --error=trim_logs/chunk_%A_%a.err
#SBATCH -A bsavoie
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --array=1-200%20

set -euo pipefail

# --- CONFIGURE FOR YOUR SYSTEM ---
# WORKDIR  : directory containing trim_dataset_worker.py and the chunk lists.
# SRC_DIR  : source archives (read-only).
# DST_DIR  : destination for slimmed archives.
WORKDIR=${WORKDIR:-/path/to/workdir}
SCRIPT="$WORKDIR/trim_dataset_worker.py"
SRC_DIR=${SRC_DIR:-/path/to/source/archives}
DST_DIR=${DST_DIR:-/path/to/slim/output}
CHUNKDIR="$WORKDIR/trim_chunks"

IDX=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
KEYS="$CHUNKDIR/chunk_${IDX}.txt"

echo "Array task: $SLURM_ARRAY_TASK_ID  Job: $SLURM_JOB_ID  Node: $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Keys file:  $KEYS  ($(wc -l < "$KEYS") zips)"
echo "Src (RO):   $SRC_DIR"
echo "Dst:        $DST_DIR"
echo "Start:      $(date)"

# --- Pre-task inode-quota guard ---
# Abort cleanly (exit 0) if scratch inode usage already exceeds the
# threshold, so SLURM does not retry and we leave headroom for
# whatever else is consuming inodes on the account.
INODE_GUARD=9700000
FILES_USED=$(lfs quota -u "$USER" "${LUSTRE_MOUNT:-/scratch}" 2>/dev/null | tail -1 | awk '{print $5}')
FILES_LIMIT=$(lfs quota -u "$USER" "${LUSTRE_MOUNT:-/scratch}" 2>/dev/null | tail -1 | awk '{print $7}')
echo "Quota:      ${FILES_USED} / ${FILES_LIMIT} inodes  (guard at ${INODE_GUARD})"
if [ -n "$FILES_USED" ] && [ "$FILES_USED" -ge "$INODE_GUARD" ]; then
    echo "ABORT: inode usage $FILES_USED has reached guard threshold $INODE_GUARD."
    echo "       Task $SLURM_ARRAY_TASK_ID exiting cleanly without processing."
    echo "       Re-submit (or use --array=<skipped task IDs>) to resume once quota frees up."
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
